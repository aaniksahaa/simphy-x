#include "internal.hpp"

#include <algorithm>
#include <array>
#include <charconv>
#include <cctype>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <map>
#include <mutex>
#include <optional>
#include <sstream>

namespace wqfm::internal {
namespace {

struct BenchmarkPhaseTotal {
  std::string name;
  long long ns = 0;
};

struct BenchmarkSubproblemRecord {
  uint64_t path = 0;
  int depth = 0;
  int active_identities = 0;
  int active_taxa = 0;
  BenchmarkIdentitySummary identities;

  bool has_initial = false;
  std::string initial_source;
  int initial_candidate_count = 0;
  double initial_score = 0.0;
  int initial_representative_order = -1;
  int initial_left_size = 0;
  int initial_right_size = 0;
  BenchmarkIdentitySummary initial_left;
  BenchmarkIdentitySummary initial_right;

  bool has_final = false;
  int final_left_size = 0;
  int final_right_size = 0;
  long long final_balance_ratio_ppm = 0;
  BenchmarkIdentitySummary final_left;
  BenchmarkIdentitySummary final_right;

  std::vector<BenchmarkPhaseTotal> phases;
};

struct BenchmarkSubproblemStore {
  std::mutex mutex;
  std::vector<BenchmarkSubproblemRecord> records;
  std::map<uint64_t, std::size_t> index_by_path;
};

BenchmarkSubproblemStore& benchmark_subproblem_store() {
  static BenchmarkSubproblemStore store;
  return store;
}

std::atomic<bool>& benchmark_telemetry_gate() {
  static std::atomic<bool> enabled{true};
  return enabled;
}

std::string path_hex(uint64_t path) {
  std::ostringstream out;
  out << "0x" << std::hex << std::setw(16) << std::setfill('0') << path;
  return out.str();
}

std::string hash_hex(uint64_t hash) {
  std::ostringstream out;
  out << "0x" << std::hex << std::setw(16) << std::setfill('0') << hash;
  return out.str();
}

class StableIdentityHasher {
 public:
  void add_u64(uint64_t value) {
    for (int byte = 0; byte < 8; ++byte) {
      add_byte(static_cast<unsigned char>(value & 0xffU));
      value >>= 8;
    }
  }

  void add_bytes(std::string_view bytes) {
    for (unsigned char byte : bytes) add_byte(byte);
  }

  uint64_t value() const { return value_; }

 private:
  void add_byte(unsigned char byte) {
    value_ ^= static_cast<uint64_t>(byte);
    value_ *= 1099511628211ULL;
  }

  uint64_t value_ = 14695981039346656037ULL;
};

std::string_view identity_key_view(const Identity& identity,
                                   std::array<char, 32>& id_buffer) {
  if (!identity.stable_key.empty()) return identity.stable_key;
  const auto result =
      std::to_chars(id_buffer.data(), id_buffer.data() + id_buffer.size(),
                    identity.id);
  return std::string_view(id_buffer.data(), result.ptr);
}

template <typename Included>
BenchmarkIdentitySummary summarize_identities(
    const std::vector<Identity>& identities,
    Included&& included) {
  BenchmarkIdentitySummary summary;
  for (std::size_t i = 0; i < identities.size(); ++i) {
    if (included(i)) ++summary.count;
  }

  // The version and fixed-width length frames make the hash independent of
  // host endianness and unambiguous for arbitrary identity key bytes.
  StableIdentityHasher hasher;
  hasher.add_u64(1);
  hasher.add_u64(static_cast<uint64_t>(summary.count));
  std::array<char, 32> id_buffer{};
  for (std::size_t i = 0; i < identities.size(); ++i) {
    if (!included(i)) continue;
    const std::string_view key = identity_key_view(identities[i], id_buffer);
    summary.key_bytes += static_cast<uint64_t>(key.size());
    hasher.add_u64(static_cast<uint64_t>(key.size()));
    hasher.add_bytes(key);
  }
  summary.ordered_hash64 = hasher.value();
  return summary;
}

BenchmarkIdentitySummary summarize_identities(
    const std::vector<Identity>& identities) {
  return summarize_identities(identities, [](std::size_t) { return true; });
}

BenchmarkIdentitySummary summarize_assignment(
    const std::vector<Side>& assignment,
    const std::vector<Identity>& identities,
    Side side) {
  return summarize_identities(identities, [&](std::size_t i) {
    return i < assignment.size() && assignment[i] == side;
  });
}

long long scaled_balance_ratio_ppm(int left_size, int right_size) {
  const int smaller = std::min(left_size, right_size);
  const int larger = std::max(left_size, right_size);
  if (smaller <= 0) return 0;
  return static_cast<long long>(
      std::llround(1000000.0 * static_cast<double>(larger) /
                   static_cast<double>(smaller)));
}

int assignment_side_count(const std::vector<Side>& assignment, Side side) {
  return static_cast<int>(
      std::count(assignment.begin(), assignment.end(), side));
}

void write_identity_summary_json(std::ostream& out,
                                 const BenchmarkIdentitySummary& summary) {
  out << "{"
      << "\"count\":" << summary.count << ","
      << "\"key_bytes\":" << summary.key_bytes << ","
      << "\"ordered_hash64\":\"" << hash_hex(summary.ordered_hash64)
      << "\"}";
}

BenchmarkSubproblemRecord* find_record_locked(BenchmarkSubproblemStore& store,
                                              uint64_t path) {
  const auto it = store.index_by_path.find(path);
  if (it == store.index_by_path.end()) return nullptr;
  return &store.records[it->second];
}

void add_phase_locked(BenchmarkSubproblemRecord& record,
                      std::string phase,
                      long long ns) {
  for (BenchmarkPhaseTotal& total : record.phases) {
    if (total.name == phase) {
      total.ns += ns;
      return;
    }
  }
  record.phases.push_back(BenchmarkPhaseTotal{std::move(phase), ns});
}

void write_optional_bipartition_json(std::ostream& out,
                                     const char* key,
                                     bool present,
                                     const std::string& source,
                                     int candidate_count,
                                     double score,
                                     int representative_order,
                                     int left_size,
                                     int right_size,
                                     const BenchmarkIdentitySummary& left,
                                     const BenchmarkIdentitySummary& right) {
  out << "      \"" << key << "\": ";
  if (!present) {
    out << "null";
    return;
  }
  out << "{"
      << "\"source\":\"" << json_escape(source) << "\","
      << "\"candidate_count\":" << candidate_count << ","
      << "\"score\":" << std::setprecision(17) << score << ","
      << "\"representative_order\":" << representative_order << ","
      << "\"left_size\":" << left_size << ","
      << "\"right_size\":" << right_size << ","
      << "\"left\":";
  write_identity_summary_json(out, left);
  out << ",\"right\":";
  write_identity_summary_json(out, right);
  out << "}";
}

void append_identity_json_string(std::string& out, const Identity& identity) {
  out.push_back('"');
  if (!identity.stable_key.empty()) {
    out += json_escape(identity.stable_key);
  } else {
    out += std::to_string(identity.id);
  }
  out.push_back('"');
}

}  // namespace

/// Converts a split bitset into side assignments.
std::vector<Side> assignment_from_split(const Bitset& split, int m) {
  std::vector<Side> assignment(m, Side::Right);
  for (int i = 0; i < m; ++i) {
    if (split.test(i)) assignment[i] = Side::Left;
  }
  return assignment;
}

/// Returns the active compressed gene bank when one is present.
const std::vector<GeneIndex>& scoring_genes(const SolverContext& ctx) {
  if (ctx.owned_active_genes != nullptr) {
    return *ctx.owned_active_genes;
  }
  return ctx.active_genes == nullptr ? base_genes(ctx) : *ctx.active_genes;
}

/// Returns the immutable base gene bank for a solver context.
const std::vector<GeneIndex>& base_genes(const SolverContext& ctx) {
  return ctx.genes_ref == nullptr ? ctx.genes : *ctx.genes_ref;
}

/// Returns the taxon registry for a solver context.
const TaxonRegistry& context_taxa(const SolverContext& ctx) {
  return ctx.taxa_ref == nullptr ? ctx.taxa : *ctx.taxa_ref;
}

/// Returns the consensus tree for a solver context.
const ConsensusTree& context_consensus(const SolverContext& ctx) {
  return ctx.consensus_ref == nullptr ? ctx.consensus : *ctx.consensus_ref;
}

/// Returns process-wide counters owned by the diagnostic inference run gate.
BenchmarkStats& benchmark_stats() {
  static BenchmarkStats stats;
  return stats;
}

bool benchmark_telemetry_enabled() noexcept {
  return benchmark_telemetry_gate().load(std::memory_order_acquire);
}

void set_benchmark_telemetry_enabled(bool enabled) noexcept {
  benchmark_telemetry_gate().store(enabled, std::memory_order_release);
}

/// Clears all benchmark counters before an inference run.
void reset_benchmark_stats() {
  BenchmarkStats& stats = benchmark_stats();
  stats.subproblems.store(0, std::memory_order_relaxed);
  stats.completed_subproblems.store(0, std::memory_order_relaxed);
  stats.max_depth.store(0, std::memory_order_relaxed);
  stats.max_active_identities.store(0, std::memory_order_relaxed);
  stats.max_active_taxa.store(0, std::memory_order_relaxed);
  stats.worst_balance_ratio_ppm.store(0, std::memory_order_relaxed);
  stats.candidate_sets.store(0, std::memory_order_relaxed);
  stats.candidates.store(0, std::memory_order_relaxed);
  stats.max_candidates.store(0, std::memory_order_relaxed);
  stats.active_consensus_candidate_sets.store(0, std::memory_order_relaxed);
  stats.active_consensus_candidates.store(0, std::memory_order_relaxed);
  stats.active_gene_split_candidates_seen.store(0, std::memory_order_relaxed);
  stats.active_gene_split_candidates_added.store(0, std::memory_order_relaxed);
  stats.active_gene_split_candidate_duplicates.store(0,
                                                     std::memory_order_relaxed);
  stats.active_total_candidates_after_expansion.store(0,
                                                      std::memory_order_relaxed);
  stats.attachment_calls.store(0, std::memory_order_relaxed);
  stats.attachment_internal_nodes.store(0, std::memory_order_relaxed);
  stats.attachment_four_color_passes.store(0, std::memory_order_relaxed);
  stats.attachment_edge_ties.store(0, std::memory_order_relaxed);
  stats.attachment_normalized_one_x_calls.store(0, std::memory_order_relaxed);
  stats.hat_subproblems.store(0, std::memory_order_relaxed);
  stats.hat_max_active_taxa.store(0, std::memory_order_relaxed);
  stats.hat_consensus_index_builds.store(0, std::memory_order_relaxed);
  stats.hat_consensus_index_nodes.store(0, std::memory_order_relaxed);
  stats.hat_active_local_taxon_entries.store(0, std::memory_order_relaxed);
  stats.hat_active_local_index_entries.store(0, std::memory_order_relaxed);
  stats.hat_virtual_query_selected_leaves.store(0,
                                                std::memory_order_relaxed);
  stats.hat_virtual_query_lca_queries.store(0, std::memory_order_relaxed);
  stats.hat_virtual_query_retained_nodes.store(0,
                                               std::memory_order_relaxed);
  stats.hat_localized_scoring_gene_nodes.store(0,
                                               std::memory_order_relaxed);
  stats.hat_localized_scoring_calls.store(0, std::memory_order_relaxed);
  stats.hat_outer_recursive_calls.store(0, std::memory_order_relaxed);
  stats.hat_outer_nonbase_calls.store(0, std::memory_order_relaxed);
  stats.hat_outer_max_depth.store(0, std::memory_order_relaxed);
  stats.hat_outer_active_identity_mass.store(0, std::memory_order_relaxed);
  stats.hat_outer_active_support_word_mass.store(0,
                                                 std::memory_order_relaxed);
  stats.hat_outer_active_taxon_mass.store(0, std::memory_order_relaxed);
  stats.hat_outer_split_identity_tests.store(0, std::memory_order_relaxed);
  stats.hat_outer_partition_identity_moves.store(0,
                                                 std::memory_order_relaxed);
  stats.hat_outer_partition_support_word_moves.store(
      0, std::memory_order_relaxed);
  stats.hat_outer_partition_taxon_moves.store(0,
                                              std::memory_order_relaxed);
  stats.hat_outer_balanced_splits.store(0, std::memory_order_relaxed);
  stats.hat_outer_maximally_skewed_splits.store(0,
                                               std::memory_order_relaxed);
  stats.hat_outer_catalog_identity_entries.store(0,
                                                 std::memory_order_relaxed);
  stats.hat_outer_catalog_taxon_entries.store(0,
                                              std::memory_order_relaxed);
  stats.hat_outer_catalog_lookup_entries.store(0,
                                               std::memory_order_relaxed);
  stats.hat_outer_tree_node_visits.store(0, std::memory_order_relaxed);
  stats.hat_outer_tree_identity_lookups.store(0,
                                             std::memory_order_relaxed);
  stats.hat_outer_tree_taxon_entries.store(0, std::memory_order_relaxed);
  stats.hat_outer_taxon_radix_entries.store(0, std::memory_order_relaxed);
  stats.hat_outer_taxon_radix_bucket_ops.store(0,
                                              std::memory_order_relaxed);
  stats.hat_outer_taxon_unique_entries.store(0,
                                             std::memory_order_relaxed);
  stats.hat_outer_root_bank_nodes.store(0, std::memory_order_relaxed);
  stats.hat_outer_root_bank_other_entries.store(0,
                                                std::memory_order_relaxed);
  stats.hat_outer_live_payload_units.store(0, std::memory_order_relaxed);
  stats.hat_outer_peak_payload_units.store(0, std::memory_order_relaxed);
  stats.hat_outer_live_identity_entries.store(0, std::memory_order_relaxed);
  stats.hat_outer_peak_identity_entries.store(0, std::memory_order_relaxed);
  stats.hat_outer_live_support_words.store(0, std::memory_order_relaxed);
  stats.hat_outer_peak_support_words.store(0, std::memory_order_relaxed);
  stats.hat_outer_live_identity_taxon_entries.store(
      0, std::memory_order_relaxed);
  stats.hat_outer_peak_identity_taxon_entries.store(
      0, std::memory_order_relaxed);
  stats.hat_outer_live_identity_key_bytes.store(0,
                                                std::memory_order_relaxed);
  stats.hat_outer_peak_identity_key_bytes.store(0,
                                                std::memory_order_relaxed);
  stats.hat_outer_live_catalog_entries.store(0, std::memory_order_relaxed);
  stats.hat_outer_peak_catalog_entries.store(0, std::memory_order_relaxed);
  stats.hat_outer_live_bank_nodes.store(0, std::memory_order_relaxed);
  stats.hat_outer_peak_bank_nodes.store(0, std::memory_order_relaxed);
  stats.hat_outer_live_bank_other_entries.store(0,
                                                std::memory_order_relaxed);
  stats.hat_outer_peak_bank_other_entries.store(0,
                                                std::memory_order_relaxed);
  stats.hat_outer_live_bank_identity_entries.store(
      0, std::memory_order_relaxed);
  stats.hat_outer_peak_bank_identity_entries.store(
      0, std::memory_order_relaxed);
  stats.hat_outer_live_continuation_identity_entries.store(
      0, std::memory_order_relaxed);
  stats.hat_outer_peak_continuation_identity_entries.store(
      0, std::memory_order_relaxed);
  stats.hat_outer_live_continuation_support_words.store(
      0, std::memory_order_relaxed);
  stats.hat_outer_peak_continuation_support_words.store(
      0, std::memory_order_relaxed);
  stats.hat_outer_live_continuation_taxon_entries.store(
      0, std::memory_order_relaxed);
  stats.hat_outer_peak_continuation_taxon_entries.store(
      0, std::memory_order_relaxed);
  stats.hat_outer_live_continuation_key_bytes.store(
      0, std::memory_order_relaxed);
  stats.hat_outer_peak_continuation_key_bytes.store(
      0, std::memory_order_relaxed);
  stats.hat_outer_live_completed_tree_nodes.store(
      0, std::memory_order_relaxed);
  stats.hat_outer_peak_completed_tree_nodes.store(
      0, std::memory_order_relaxed);
  stats.hat_outer_live_completed_tree_adjacencies.store(
      0, std::memory_order_relaxed);
  stats.hat_outer_peak_completed_tree_adjacencies.store(
      0, std::memory_order_relaxed);
  stats.hat_outer_live_taxon_scratch_entries.store(
      0, std::memory_order_relaxed);
  stats.hat_outer_peak_taxon_scratch_entries.store(
      0, std::memory_order_relaxed);
  stats.hat_outer_live_accounting_violations.store(
      0, std::memory_order_relaxed);
  stats.hat_consensus_edges_scored.store(0, std::memory_order_relaxed);
  stats.hat_hdt_banks_built.store(0, std::memory_order_relaxed);
  stats.hat_hdt_bank_build_ns.store(0, std::memory_order_relaxed);
  stats.hat_hdt_recolor_ops.store(0, std::memory_order_relaxed);
  stats.hat_hdt_nodes_recomputed.store(0, std::memory_order_relaxed);
  stats.hat_hdt_queries.store(0, std::memory_order_relaxed);
  stats.hat_reduced_hdt_nodes.store(0, std::memory_order_relaxed);
  stats.hat_contract_calls.store(0, std::memory_order_relaxed);
  stats.hat_extract_calls.store(0, std::memory_order_relaxed);
  stats.hat_contract_input_nodes.store(0, std::memory_order_relaxed);
  stats.hat_contract_output_nodes.store(0, std::memory_order_relaxed);
  stats.hat_contract_live_leaves.store(0, std::memory_order_relaxed);
  stats.hat_contract_max_output_per_live_ppm.store(0,
                                                   std::memory_order_relaxed);
  stats.hat_extract_input_nodes.store(0, std::memory_order_relaxed);
  stats.hat_extract_output_nodes.store(0, std::memory_order_relaxed);
  stats.hat_extract_live_leaves.store(0, std::memory_order_relaxed);
  stats.hat_extract_max_output_per_live_ppm.store(0,
                                                  std::memory_order_relaxed);
  stats.hat_extract_source_hdt_nodes.store(0, std::memory_order_relaxed);
  stats.hat_extract_marked_hdt_nodes.store(0, std::memory_order_relaxed);
  stats.hat_extract_hanging_components.store(0, std::memory_order_relaxed);
  stats.hat_extract_intermediate_units.store(0, std::memory_order_relaxed);
  stats.hat_extract_live_leaves_exact.store(0, std::memory_order_relaxed);
  stats.hat_contract_light_calls.store(0, std::memory_order_relaxed);
  stats.hat_contract_heavy_calls.store(0, std::memory_order_relaxed);
  stats.hat_contract_input_units.store(0, std::memory_order_relaxed);
  stats.hat_contract_output_units.store(0, std::memory_order_relaxed);
  stats.hat_contract_live_leaves_exact.store(0, std::memory_order_relaxed);
  stats.hat_contract_corridors_emitted.store(0, std::memory_order_relaxed);
  stats.hat_contract_background_units_absorbed.store(
      0, std::memory_order_relaxed);
  stats.hat_contract_polytomy_arms_folded.store(0,
                                                std::memory_order_relaxed);
  stats.hat_contract_cc_ops.store(0, std::memory_order_relaxed);
  stats.hat_contract_gg_ops.store(0, std::memory_order_relaxed);
  stats.hat_contract_ig_ops.store(0, std::memory_order_relaxed);
  stats.hat_contract_c_to_g_ops.store(0, std::memory_order_relaxed);
  stats.hat_contract_bound_violations.store(0, std::memory_order_relaxed);
  stats.hat_heavy_input_units.store(0, std::memory_order_relaxed);
  stats.hat_heavy_output_units.store(0, std::memory_order_relaxed);
  stats.hat_heavy_calls.store(0, std::memory_order_relaxed);
  stats.hat_heavy_shrink_violations.store(0, std::memory_order_relaxed);
  stats.hat_hdt_initial_atoms.store(0, std::memory_order_relaxed);
  stats.hat_hdt_rebuild_nodes.store(0, std::memory_order_relaxed);
  stats.hat_hdt_local_balance_violations.store(0,
                                               std::memory_order_relaxed);
  stats.hat_hdt_max_update_height.store(0, std::memory_order_relaxed);
  stats.hat_recolor_changed_leaves.store(0, std::memory_order_relaxed);
  stats.hat_recolor_marked_nodes.store(0, std::memory_order_relaxed);
  stats.hat_dp_transition_evals.store(0, std::memory_order_relaxed);
  stats.hat_dp_scalar_fmas.store(0, std::memory_order_relaxed);
  stats.hat_extract_scratch_peak_bytes.store(0, std::memory_order_relaxed);
  stats.hat_contract_scratch_peak_bytes.store(0, std::memory_order_relaxed);
  stats.hat_state_arena_peak_bytes.store(0, std::memory_order_relaxed);
  stats.hat_worker_total_peak_bytes.store(0, std::memory_order_relaxed);
  stats.hat_score_partial_peak_bytes.store(0, std::memory_order_relaxed);
  stats.hat_contract_mark_ns.store(0, std::memory_order_relaxed);
  stats.hat_contract_background_fold_ns.store(0,
                                              std::memory_order_relaxed);
  stats.hat_contract_corridor_fold_ns.store(0,
                                            std::memory_order_relaxed);
  stats.hat_contract_emit_ns.store(0, std::memory_order_relaxed);
  stats.hat_hdt_rebuild_ns.store(0, std::memory_order_relaxed);
  stats.complexity_contract_satisfied.store(1, std::memory_order_relaxed);
  stats.hat_frozen_segments.store(0, std::memory_order_relaxed);
  stats.hat_frozen_segment_source_degree_sum.store(0,
                                                   std::memory_order_relaxed);
  stats.hat_recompute_child_folds.store(0, std::memory_order_relaxed);
  stats.hat_max_recompute_degree.store(0, std::memory_order_relaxed);
  stats.hat_split_score_ns.store(0, std::memory_order_relaxed);
  stats.hat_certificate_index_ns.store(0, std::memory_order_relaxed);
  stats.hat_certificate_ns.store(0, std::memory_order_relaxed);
  stats.hat_certificate_tests.store(0, std::memory_order_relaxed);
  stats.hat_certificate_accepts.store(0, std::memory_order_relaxed);
  stats.hat_certificate_rejects.store(0, std::memory_order_relaxed);
  stats.hat_certificate_negative_accepts.store(0, std::memory_order_relaxed);
  stats.hat_certificate_max_negative_gap_twice_k.store(0, std::memory_order_relaxed);
  stats.hat_certificate_singletons.store(0, std::memory_order_relaxed);
  stats.hat_certificate_distance_queries.store(0, std::memory_order_relaxed);
  stats.hat_one_x_calls.store(0, std::memory_order_relaxed);
  stats.hat_one_x_medians.store(0, std::memory_order_relaxed);
  stats.hat_one_x_ns.store(0, std::memory_order_relaxed);
  stats.hat_root_source_scoring_calls.store(0, std::memory_order_relaxed);
  stats.hat_child_gene_bank_reuses.store(0, std::memory_order_relaxed);
  stats.hat_child_bank_state_attempts.store(0, std::memory_order_relaxed);
  stats.hat_child_bank_state_accepts.store(0, std::memory_order_relaxed);
  stats.hat_child_bank_state_heavy_reuses.store(0,
                                                std::memory_order_relaxed);
  stats.hat_child_bank_state_pair_violations.store(
      0, std::memory_order_relaxed);
  stats.hat_child_bank_state_input_nodes.store(0,
                                               std::memory_order_relaxed);
  stats.hat_child_bank_state_output_nodes.store(0,
                                                std::memory_order_relaxed);
  stats.hat_child_bank_state_max_accepted_ratio_ppm.store(
      0, std::memory_order_relaxed);
  stats.hat_max_scoring_source_identities.store(0, std::memory_order_relaxed);
  stats.hat_max_scoring_source_inflation_ppm.store(
      0, std::memory_order_relaxed);
  stats.hat_scoring_source_inflation_violations.store(
      0, std::memory_order_relaxed);
  stats.hat_structural_serial_child_subproblems.store(
      0, std::memory_order_relaxed);
  stats.hat_structural_parallel_child_subproblems.store(
      0, std::memory_order_relaxed);
  stats.hat_structural_parallel_rejections.store(
      0, std::memory_order_relaxed);
  stats.hat_structural_parallel_physical_rejections.store(
      0, std::memory_order_relaxed);
  stats.hat_active_structural_parallel_slots.store(
      0, std::memory_order_relaxed);
  stats.hat_max_structural_parallel_slots.store(
      0, std::memory_order_relaxed);
  stats.hat_nni_passes_executed.store(0, std::memory_order_relaxed);
  stats.hat_nni_edges_scored.store(0, std::memory_order_relaxed);
  stats.hat_nni_full_score_passes.store(0, std::memory_order_relaxed);
  stats.hat_nni_local_score_passes.store(0, std::memory_order_relaxed);
  stats.hat_nni_local_edges_scored.store(0, std::memory_order_relaxed);
  stats.hat_nni_moves_accepted.store(0, std::memory_order_relaxed);
  stats.hat_nni_moves_rejected_by_threshold.store(
      0, std::memory_order_relaxed);
  stats.hat_nni_peak_proposals.store(0, std::memory_order_relaxed);
  stats.hat_nni_max_qsplit_leaves.store(0, std::memory_order_relaxed);
  stats.hat_nni_qsplit_owner_bytes_avoided.store(0,
                                                 std::memory_order_relaxed);
  stats.hat_nni_gain_threshold.store(0.0, std::memory_order_relaxed);
  stats.hat_nni_last_accepted_gain.store(0.0, std::memory_order_relaxed);
  stats.hat_nni_max_accepted_gain.store(0.0, std::memory_order_relaxed);
  stats.hat_nni_verified_gain.store(0.0, std::memory_order_relaxed);
  stats.hat_nni_prefix_gate_stops.store(0, std::memory_order_relaxed);
  stats.hat_nni_prefix_gate_quartet_stops.store(0,
                                                std::memory_order_relaxed);
  stats.hat_nni_prefix_gate_split_stops.store(0, std::memory_order_relaxed);
  stats.hat_nni_split_support_delta.store(0, std::memory_order_relaxed);
  stats.hat_nni_split_support_unique_splits.store(0,
                                                  std::memory_order_relaxed);
  stats.hat_nni_split_support_ns.store(0, std::memory_order_relaxed);
  stats.hat_nni_ns.store(0, std::memory_order_relaxed);
  stats.hat_peak_bank_bytes.store(0, std::memory_order_relaxed);
  stats.active_child_compressed_bank_identity_slots.store(
      0, std::memory_order_relaxed);
  stats.max_active_child_compressed_bank_identity_slots.store(
      0, std::memory_order_relaxed);
  stats.consensus_distinct_splits.store(0, std::memory_order_relaxed);
  stats.consensus_accepted_splits.store(0, std::memory_order_relaxed);
  stats.consensus_hash_low.store(0, std::memory_order_relaxed);
  stats.mre_candidates_checked.store(0, std::memory_order_relaxed);
  stats.mre_candidates_rejected.store(0, std::memory_order_relaxed);
  stats.mre_case_a_accepts.store(0, std::memory_order_relaxed);
  stats.mre_case_b_accepts.store(0, std::memory_order_relaxed);
  stats.mre_moved_children.store(0, std::memory_order_relaxed);
  stats.mre_max_child_degree.store(0, std::memory_order_relaxed);
  stats.mre_radix_passes.store(0, std::memory_order_relaxed);
  stats.mre_radix_items.store(0, std::memory_order_relaxed);
  stats.mre_child_index_operations.store(0, std::memory_order_relaxed);
  stats.mre_child_index_rotations.store(0, std::memory_order_relaxed);
  stats.mre_child_index_max_height.store(0, std::memory_order_relaxed);
  stats.mre_child_index_height_violations.store(
      0, std::memory_order_relaxed);
  stats.mre_child_index_enumeration_visits.store(
      0, std::memory_order_relaxed);
  stats.mre_topology_slots_allocated.store(0, std::memory_order_relaxed);
  stats.mre_topology_peak_live_slots.store(0, std::memory_order_relaxed);
  stats.compression_materializations.store(0, std::memory_order_relaxed);
  stats.active_compressed_banks.store(0, std::memory_order_relaxed);
  stats.max_active_compressed_banks.store(0, std::memory_order_relaxed);
  stats.active_child_compressed_banks.store(0, std::memory_order_relaxed);
  stats.max_active_child_compressed_banks.store(0, std::memory_order_relaxed);
  stats.max_child_compressed_bank_bytes.store(0, std::memory_order_relaxed);
  stats.choose_initial_ns.store(0, std::memory_order_relaxed);
  stats.attachment_range_add_ns.store(0, std::memory_order_relaxed);
  stats.attachment_score_ns.store(0, std::memory_order_relaxed);
  stats.attachment_total_ns.store(0, std::memory_order_relaxed);
  stats.parse_input_ns.store(0, std::memory_order_relaxed);
  stats.taxon_registry_ns.store(0, std::memory_order_relaxed);
  stats.taxon_registry_physical_labels.store(0,
                                             std::memory_order_relaxed);
  stats.taxon_registry_unique_labels.store(0,
                                           std::memory_order_relaxed);
  stats.taxon_registry_sort_items.store(0, std::memory_order_relaxed);
  stats.gene_index_ns.store(0, std::memory_order_relaxed);
  stats.consensus_count_ns.store(0, std::memory_order_relaxed);
  stats.consensus_assembly_ns.store(0, std::memory_order_relaxed);
  stats.startup_ns.store(0, std::memory_order_relaxed);
  stats.compression_materialization_ns.store(0, std::memory_order_relaxed);
  stats.cpu_work_units.store(0, std::memory_order_relaxed);
  BenchmarkSubproblemStore& store = benchmark_subproblem_store();
  std::lock_guard<std::mutex> lock(store.mutex);
  store.records.clear();
  store.index_by_path.clear();
}

/// Atomically increments a benchmark counter.
void benchmark_counter_add(std::atomic<long long>& counter, long long value) {
  if (!benchmark_telemetry_enabled()) return;
  counter.fetch_add(value, std::memory_order_relaxed);
}

/// Atomically updates a benchmark maximum.
void benchmark_counter_max(std::atomic<long long>& counter, long long value) {
  if (!benchmark_telemetry_enabled()) return;
  long long current = counter.load(std::memory_order_relaxed);
  while (current < value &&
         !counter.compare_exchange_weak(current, value,
                                        std::memory_order_relaxed,
                                        std::memory_order_relaxed)) {
  }
}

namespace {

void benchmark_live_underflow(long long previous, long long value) {
  if (value < 0 || previous < value) {
    benchmark_counter_add(
        benchmark_stats().hat_outer_live_accounting_violations);
  }
}

}  // namespace

void benchmark_live_payload_acquire(std::atomic<long long>& current,
                                    std::atomic<long long>& peak,
                                    long long value) {
  if (value <= 0) return;
  const long long category =
      current.fetch_add(value, std::memory_order_relaxed) + value;
  benchmark_counter_max(peak, category);
  BenchmarkStats& stats = benchmark_stats();
  const long long total = stats.hat_outer_live_payload_units.fetch_add(
                              value, std::memory_order_relaxed) +
                          value;
  benchmark_counter_max(stats.hat_outer_peak_payload_units, total);
}

void benchmark_live_payload_release(std::atomic<long long>& current,
                                    long long value) {
  if (value <= 0) return;
  const long long previous =
      current.fetch_sub(value, std::memory_order_relaxed);
  benchmark_live_underflow(previous, value);
  BenchmarkStats& stats = benchmark_stats();
  const long long total_previous = stats.hat_outer_live_payload_units.fetch_sub(
      value, std::memory_order_relaxed);
  benchmark_live_underflow(total_previous, value);
}

void benchmark_live_gauge_acquire(std::atomic<long long>& current,
                                  std::atomic<long long>& peak,
                                  long long value) {
  if (value <= 0) return;
  const long long observed =
      current.fetch_add(value, std::memory_order_relaxed) + value;
  benchmark_counter_max(peak, observed);
}

void benchmark_live_gauge_release(std::atomic<long long>& current,
                                  long long value) {
  if (value <= 0) return;
  const long long previous =
      current.fetch_sub(value, std::memory_order_relaxed);
  benchmark_live_underflow(previous, value);
}

/// Adds elapsed time to a benchmark timing counter.
void benchmark_timer_add(std::atomic<long long>& counter,
                         std::chrono::nanoseconds elapsed) {
  if (!benchmark_telemetry_enabled()) return;
  counter.fetch_add(elapsed.count(), std::memory_order_relaxed);
}

void benchmark_subproblem_start(uint64_t path,
                                int depth,
                                int active_identities,
                                int active_taxa,
                                const std::vector<Identity>& identities) {
  if (!benchmark_telemetry_enabled()) return;
  const BenchmarkIdentitySummary identity_summary =
      summarize_identities(identities);
  BenchmarkSubproblemStore& store = benchmark_subproblem_store();
  std::lock_guard<std::mutex> lock(store.mutex);
  if (store.index_by_path.find(path) != store.index_by_path.end()) return;
  BenchmarkSubproblemRecord record;
  record.path = path;
  record.depth = depth;
  record.active_identities = active_identities;
  record.active_taxa = active_taxa;
  record.identities = identity_summary;
  store.index_by_path.emplace(path, store.records.size());
  store.records.push_back(std::move(record));
}

void benchmark_subproblem_add_phase(uint64_t path,
                                    std::string phase,
                                    std::chrono::nanoseconds elapsed) {
  if (!benchmark_telemetry_enabled()) return;
  BenchmarkSubproblemStore& store = benchmark_subproblem_store();
  std::lock_guard<std::mutex> lock(store.mutex);
  BenchmarkSubproblemRecord* record = find_record_locked(store, path);
  if (record == nullptr) return;
  add_phase_locked(*record, std::move(phase), elapsed.count());
}

void benchmark_subproblem_set_initial(uint64_t path,
                                      std::string source,
                                      int candidate_count,
                                      double score,
                                      int representative_order,
                                      const std::vector<Side>& assignment,
                                      const std::vector<Identity>& identities) {
  if (!benchmark_telemetry_enabled()) return;
  const int left_size = assignment_side_count(assignment, Side::Left);
  const int right_size = assignment_side_count(assignment, Side::Right);
  const BenchmarkIdentitySummary left_summary =
      summarize_assignment(assignment, identities, Side::Left);
  const BenchmarkIdentitySummary right_summary =
      summarize_assignment(assignment, identities, Side::Right);
  BenchmarkSubproblemStore& store = benchmark_subproblem_store();
  std::lock_guard<std::mutex> lock(store.mutex);
  BenchmarkSubproblemRecord* record = find_record_locked(store, path);
  if (record == nullptr) return;
  record->has_initial = true;
  record->initial_source = std::move(source);
  record->initial_candidate_count = candidate_count;
  record->initial_score = score;
  record->initial_representative_order = representative_order;
  record->initial_left_size = left_size;
  record->initial_right_size = right_size;
  record->initial_left = left_summary;
  record->initial_right = right_summary;
}

void benchmark_subproblem_set_final(uint64_t path,
                                    const std::vector<Identity>& left,
                                    const std::vector<Identity>& right) {
  if (!benchmark_telemetry_enabled()) return;
  const int left_size = static_cast<int>(left.size());
  const int right_size = static_cast<int>(right.size());
  const long long balance_ratio_ppm =
      scaled_balance_ratio_ppm(left_size, right_size);
  benchmark_counter_max(benchmark_stats().worst_balance_ratio_ppm,
                        balance_ratio_ppm);
  const BenchmarkIdentitySummary left_summary = summarize_identities(left);
  const BenchmarkIdentitySummary right_summary = summarize_identities(right);
  BenchmarkSubproblemStore& store = benchmark_subproblem_store();
  std::lock_guard<std::mutex> lock(store.mutex);
  BenchmarkSubproblemRecord* record = find_record_locked(store, path);
  if (record == nullptr) return;
  record->has_final = true;
  record->final_left_size = left_size;
  record->final_right_size = right_size;
  record->final_balance_ratio_ppm = balance_ratio_ppm;
  record->final_left = left_summary;
  record->final_right = right_summary;
}

void write_benchmark_subproblem_records_json(std::ostream& out) {
  BenchmarkSubproblemStore& store = benchmark_subproblem_store();
  std::lock_guard<std::mutex> lock(store.mutex);
  std::vector<const BenchmarkSubproblemRecord*> ordered_records;
  ordered_records.reserve(store.records.size());
  for (const BenchmarkSubproblemRecord& record : store.records) {
    ordered_records.push_back(&record);
  }
  std::sort(ordered_records.begin(), ordered_records.end(),
            [](const BenchmarkSubproblemRecord* lhs,
               const BenchmarkSubproblemRecord* rhs) {
              return lhs->path < rhs->path;
            });

  out << "  \"subproblem_record_schema\": 2,\n";
  out << "  \"subproblem_identity_hash\": "
         "\"fnv1a64-length-framed-ordered-v1\",\n";
  out << "  \"subproblem_records\": [\n";
  for (std::size_t i = 0; i < ordered_records.size(); ++i) {
    const BenchmarkSubproblemRecord& record = *ordered_records[i];
    out << "    {\n";
    out << "      \"path\": " << record.path << ",\n";
    out << "      \"path_hex\": \"" << path_hex(record.path) << "\",\n";
    out << "      \"depth\": " << record.depth << ",\n";
    out << "      \"identity_count\": " << record.active_identities << ",\n";
    out << "      \"active_taxa\": " << record.active_taxa << ",\n";
    out << "      \"identities\": ";
    write_identity_summary_json(out, record.identities);
    out << ",\n";
    write_optional_bipartition_json(
        out, "initial_bipartition", record.has_initial,
        record.initial_source, record.initial_candidate_count,
        record.initial_score, record.initial_representative_order,
        record.initial_left_size, record.initial_right_size,
        record.initial_left, record.initial_right);
    out << ",\n";
    out << "      \"final_bipartition\": ";
    if (!record.has_final) {
      out << "null";
    } else {
      out << "{"
          << "\"left_size\":" << record.final_left_size << ","
          << "\"right_size\":" << record.final_right_size << ","
          << "\"balance_ratio_ppm\":" << record.final_balance_ratio_ppm << ","
          << "\"left\":";
      write_identity_summary_json(out, record.final_left);
      out << ",\"right\":";
      write_identity_summary_json(out, record.final_right);
      out << "}";
    }
    out << ",\n";
    out << "      \"phase_ns\": {";
    std::vector<const BenchmarkPhaseTotal*> ordered_phases;
    ordered_phases.reserve(record.phases.size());
    for (const BenchmarkPhaseTotal& phase : record.phases) {
      ordered_phases.push_back(&phase);
    }
    std::sort(ordered_phases.begin(), ordered_phases.end(),
              [](const BenchmarkPhaseTotal* lhs,
                 const BenchmarkPhaseTotal* rhs) {
                return lhs->name < rhs->name;
              });
    for (std::size_t p = 0; p < ordered_phases.size(); ++p) {
      if (p != 0) out << ", ";
      out << "\"" << json_escape(ordered_phases[p]->name) << "\": "
          << ordered_phases[p]->ns;
    }
    out << "}\n";
    out << "    }" << (i + 1 == ordered_records.size() ? "\n" : ",\n");
  }
  out << "  ]";
}

/// Returns a JSON/debug name for a progress mode.
const char* progress_mode_name(ProgressMode mode) {
  switch (mode) {
    case ProgressMode::Off:
      return "off";
    case ProgressMode::Auto:
      return "auto";
    case ProgressMode::On:
      return "on";
  }
  return "unknown";
}

/// Tests whether tracing is enabled at a depth.
bool trace_enabled_at(const Options& options, int depth) {
  return !options.debug_trace.empty() &&
         (options.debug_max_depth < 0 || depth <= options.debug_max_depth);
}

/// Escapes a string for JSON output.
std::string json_escape(const std::string& value) {
  std::string out;
  out.reserve(value.size() + 8);
  for (char ch : value) {
    switch (ch) {
      case '\\':
        out += "\\\\";
        break;
      case '"':
        out += "\\\"";
        break;
      case '\n':
        out += "\\n";
        break;
      case '\r':
        out += "\\r";
        break;
      case '\t':
        out += "\\t";
        break;
      default:
        if (static_cast<unsigned char>(ch) < 0x20) {
          out += '?';
        } else {
          out += ch;
        }
        break;
    }
  }
  return out;
}

/// Emits one JSONL debug event if tracing is enabled for this recursion depth.
void trace_json_line(const Options& options, int depth, const std::string& body) {
  if (!trace_enabled_at(options, depth)) return;
  static std::mutex trace_mutex;
  std::lock_guard<std::mutex> lock(trace_mutex);
  std::ofstream out(options.debug_trace, std::ios::app);
  if (!out) return;
  out << body << '\n';
}

/// Returns the trace key for an identity.
std::string identity_trace_key(const Identity& identity) {
  if (!identity.stable_key.empty()) return identity.stable_key;
  return std::to_string(identity.id);
}

/// Serializes identities as a JSON array.
std::string identity_array_json(const std::vector<Identity>& identities) {
  std::string out;
  out.reserve(2 + identities.size() * 8);
  out.push_back('[');
  for (int i = 0; i < static_cast<int>(identities.size()); ++i) {
    if (i != 0) out.push_back(',');
    append_identity_json_string(out, identities[i]);
  }
  out.push_back(']');
  return out;
}

/// Serializes one side of a split.
std::string split_array_json(const Bitset& split,
                             const std::vector<Identity>& identities,
                             bool left_side) {
  std::string out;
  out.reserve(2 + identities.size() * 4);
  out.push_back('[');
  bool first = true;
  for (int i = 0; i < static_cast<int>(identities.size()); ++i) {
    if (split.test(i) != left_side) continue;
    if (!first) out.push_back(',');
    first = false;
    append_identity_json_string(out, identities[i]);
  }
  out.push_back(']');
  return out;
}

/// Serializes one side of an assignment.
std::string assignment_array_json(const std::vector<Side>& assignment,
                                  const std::vector<Identity>& identities,
                                  Side side) {
  std::string out;
  out.reserve(2 + identities.size() * 4);
  out.push_back('[');
  bool first = true;
  for (int i = 0; i < static_cast<int>(assignment.size()); ++i) {
    if (assignment[i] != side) continue;
    if (!first) out.push_back(',');
    first = false;
    append_identity_json_string(out, identities[i]);
  }
  out.push_back(']');
  return out;
}

/// Serializes all taxa in a set of identities.
std::string identity_taxa_array_json(const std::vector<Identity>& identities) {
  std::vector<std::string> taxa;
  for (const Identity& identity : identities) {
    for (int taxon : identity.support_taxa) {
      taxa.push_back(std::to_string(taxon));
    }
  }
  std::sort(taxa.begin(), taxa.end());
  taxa.erase(std::unique(taxa.begin(), taxa.end()), taxa.end());
  std::ostringstream out;
  out << "[";
  for (int i = 0; i < static_cast<int>(taxa.size()); ++i) {
    if (i != 0) out << ",";
    out << "\"" << taxa[i] << "\"";
  }
  out << "]";
  return out.str();
}

/// Serializes taxa selected by one side of an assignment.
std::string assignment_taxa_array_json(const std::vector<Side>& assignment,
                                       const std::vector<Identity>& identities,
                                       Side side) {
  std::vector<Identity> selected;
  for (int i = 0; i < static_cast<int>(assignment.size()); ++i) {
    if (assignment[i] == side) selected.push_back(identities[i]);
  }
  return identity_taxa_array_json(selected);
}

}  // namespace wqfm::internal
