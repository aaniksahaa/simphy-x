#include "internal.hpp"
#include "hat/hat_split_certificate.hpp"

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <condition_variable>
#include <fstream>
#include <iterator>
#include <memory>
#include <mutex>
#include <numeric>
#include <optional>
#include <set>
#include <sstream>
#include <string_view>
#include <thread>
#include <type_traits>

namespace wqfm::internal {
namespace {

struct InferenceRunState {
  std::mutex mutex;
  std::condition_variable changed;
  int active_shared_runs = 0;
  int waiting_writers = 0;
  int configured_workers = 0;
  bool exclusive_run = false;
};

InferenceRunState& inference_run_state() {
  static InferenceRunState state;
  return state;
}

int requested_worker_count(const Options& options, int worker_cap) {
  const int requested =
      options.threads > 0
          ? options.threads
          : static_cast<int>(std::thread::hardware_concurrency());
  return std::min(std::max(1, worker_cap), std::max(1, requested));
}

bool requires_exclusive_inference_run(const Options& options) {
  return options.progress != ProgressMode::Off ||
         !options.benchmark_json.empty() || !options.dump_consensus.empty() ||
         !options.debug_trace.empty() ||
         !options.hat_nni_pass_tree_dir.empty() ||
         !options.hat_debug_trace.empty() || options.debug_stop_after_root ||
         debug_timing_enabled();
}

/// Coordinates the shared executor and process-global diagnostics. Quiet runs
/// may overlap only while the executor configuration is stable.
class ScopedInferenceRun {
 public:
  ScopedInferenceRun(const Options& options, int worker_cap)
      : requested_workers_(requested_worker_count(options, worker_cap)),
        exclusive_(requires_exclusive_inference_run(options)) {
    if (exclusive_) {
      acquire_exclusive(options);
    } else {
      acquire_shared(options);
    }
  }

  ScopedInferenceRun(const ScopedInferenceRun&) = delete;
  ScopedInferenceRun& operator=(const ScopedInferenceRun&) = delete;

  ~ScopedInferenceRun() noexcept {
    InferenceRunState& state = inference_run_state();
    {
      std::lock_guard<std::mutex> lock(state.mutex);
      if (exclusive_) {
        state.exclusive_run = false;
      } else {
        --state.active_shared_runs;
        if (state.active_shared_runs == 0) {
          set_benchmark_telemetry_enabled(true);
        }
      }
    }
    state.changed.notify_all();
  }

  int workers() const noexcept { return workers_; }

  bool owns_process_telemetry() const noexcept { return exclusive_; }

  // Structural child reservations currently live in BenchmarkStats. Shared
  // runs keep recursive child scheduling serial so another run's resource
  // gauges cannot affect this run's scheduling decisions.
  bool permits_recursive_subflows() const noexcept { return exclusive_; }

 private:
  int configure_requested_workers(const Options& options) const {
    Options effective = options;
    effective.threads = requested_workers_;
    return configure_threads(effective);
  }

  void acquire_exclusive(const Options& options) {
    InferenceRunState& state = inference_run_state();
    std::unique_lock<std::mutex> lock(state.mutex);
    ++state.waiting_writers;
    state.changed.wait(lock, [&] {
      return !state.exclusive_run && state.active_shared_runs == 0;
    });
    --state.waiting_writers;
    state.exclusive_run = true;
    set_benchmark_telemetry_enabled(true);
    try {
      workers_ = configure_requested_workers(options);
      state.configured_workers = workers_;
    } catch (...) {
      state.exclusive_run = false;
      lock.unlock();
      state.changed.notify_all();
      throw;
    }
  }

  void acquire_shared(const Options& options) {
    InferenceRunState& state = inference_run_state();
    std::unique_lock<std::mutex> lock(state.mutex);
    bool waiting_for_reconfiguration = false;
    for (;;) {
      if (!state.exclusive_run && state.active_shared_runs > 0 &&
          state.configured_workers == requested_workers_ &&
          state.waiting_writers == 0) {
        ++state.active_shared_runs;
        workers_ = state.configured_workers;
        return;
      }

      if (!state.exclusive_run && state.active_shared_runs == 0 &&
          (state.waiting_writers == 0 || waiting_for_reconfiguration)) {
        if (waiting_for_reconfiguration) {
          --state.waiting_writers;
          waiting_for_reconfiguration = false;
        }

        if (state.configured_workers == requested_workers_ &&
            get_num_workers() == requested_workers_) {
          set_benchmark_telemetry_enabled(false);
          state.active_shared_runs = 1;
          workers_ = state.configured_workers;
          lock.unlock();
          state.changed.notify_all();
          return;
        }

        state.exclusive_run = true;
        try {
          workers_ = configure_requested_workers(options);
          state.configured_workers = workers_;
        } catch (...) {
          state.exclusive_run = false;
          lock.unlock();
          state.changed.notify_all();
          throw;
        }
        state.exclusive_run = false;
        set_benchmark_telemetry_enabled(false);
        state.active_shared_runs = 1;
        lock.unlock();
        state.changed.notify_all();
        return;
      }

      if (!waiting_for_reconfiguration &&
          state.configured_workers != requested_workers_) {
        ++state.waiting_writers;
        waiting_for_reconfiguration = true;
      }
      state.changed.wait(lock);
    }
  }

  int requested_workers_ = 1;
  int workers_ = 1;
  bool exclusive_ = false;
};

/// Ensures an exceptional startup/consensus path cannot release an exclusive
/// run while a progress heartbeat still accesses process-wide diagnostics.
class ScopedInferenceProgressFinish {
 public:
  explicit ScopedInferenceProgressFinish(bool active) : active_(active) {}
  ScopedInferenceProgressFinish(const ScopedInferenceProgressFinish&) = delete;
  ScopedInferenceProgressFinish& operator=(
      const ScopedInferenceProgressFinish&) = delete;

  ~ScopedInferenceProgressFinish() noexcept {
    if (!active_) return;
    try {
      progress_finish();
    } catch (...) {
    }
  }

  void finish() {
    if (!active_) return;
    progress_finish();
    active_ = false;
  }

 private:
  bool active_ = false;
};

const char* algorithm_name(AlgorithmMode mode) {
  switch (mode) {
    case AlgorithmMode::HatV1:
      return "hat-v1";
    case AlgorithmMode::HatConsensus:
      return "hat-consensus";
    case AlgorithmMode::HatConsensusNniV1:
      return "hat-consensus-nni-v1";
  }
  return "unknown";
}

const char* taxon_order_name(TaxonOrderMode mode) {
  switch (mode) {
    case TaxonOrderMode::Sorted:
      return "sorted";
    case TaxonOrderMode::JavaHashSet:
      return "java-hashset";
  }
  return "unknown";
}

const char* hat_split_name(HatSplitMode mode) {
  switch (mode) {
    case HatSplitMode::Raw:
      return "raw";
    case HatSplitMode::ConsensusBalanced:
      return "consensus-balanced";
    case HatSplitMode::ConsensusPolytomyBalanced:
      return "consensus-polytomy-balanced";
  }
  return "unknown";
}

const char* hat_one_x_name(HatOneXMode mode) {
  switch (mode) {
    case HatOneXMode::Raw:
      return "raw";
    case HatOneXMode::Normalized:
      return "normalized";
  }
  return "unknown";
}

int java_string_hash(const std::string& value) {
  int32_t hash = 0;
  for (unsigned char ch : value) {
    hash = static_cast<int32_t>(31 * hash + ch);
  }
  return hash;
}

std::vector<std::string> java_hashset_taxon_order(
    const std::vector<std::string>& labels) {
  std::vector<std::string> unique;
  unique.reserve(labels.size());
  std::set<std::string_view, std::less<>> seen;
  for (const std::string& label : labels) {
    if (seen.insert(label).second) {
      unique.push_back(label);
    }
  }

  int capacity = 16;
  while (static_cast<int>(unique.size()) > capacity * 3 / 4) {
    capacity *= 2;
  }
  std::stable_sort(unique.begin(), unique.end(), [&](const std::string& a,
                                                     const std::string& b) {
    const uint32_t ah = static_cast<uint32_t>(java_string_hash(a));
    const uint32_t bh = static_cast<uint32_t>(java_string_hash(b));
    const uint32_t aspread = ah ^ (ah >> 16);
    const uint32_t bspread = bh ^ (bh >> 16);
    return (aspread & static_cast<uint32_t>(capacity - 1)) <
           (bspread & static_cast<uint32_t>(capacity - 1));
  });
  return unique;
}

TaxonRegistry build_registry_with_order(std::vector<std::string> labels,
                                        TaxonOrderMode order) {
  if (order == TaxonOrderMode::JavaHashSet) {
    return TaxonRegistry::from_ordered_unique(
        java_hashset_taxon_order(labels));
  }
  return TaxonRegistry(std::move(labels));
}

std::string startup_settings_summary(const Options& options) {
  std::ostringstream out;
  out << "settings algorithm=" << algorithm_name(options.algorithm)
      << " taxon_order=" << taxon_order_name(options.taxon_order)
      << " hat_split=" << hat_split_name(options.hat_split_mode)
      << " hat_one_third_singleton="
      << (options.hat_one_third_singleton ? "true" : "false")
      << " hat_split_certificate="
      << (options.hat_split_certificate ? "true" : "false")
      << " hat_one_x=" << hat_one_x_name(options.hat_one_x_mode)
      << " structural_budget_factor="
      << options.hat_structural_budget_factor;
  return out.str();
}

/// Adds an undirected edge to a solver tree.
void add_edge(SolverTree& tree, int a, int b) {
  if (std::find(tree.nodes[a].neighbors.begin(), tree.nodes[a].neighbors.end(), b) ==
      tree.nodes[a].neighbors.end()) {
    tree.nodes[a].neighbors.push_back(b);
  }
  if (std::find(tree.nodes[b].neighbors.begin(), tree.nodes[b].neighbors.end(), a) ==
      tree.nodes[b].neighbors.end()) {
    tree.nodes[b].neighbors.push_back(a);
  }
}

/// Removes an undirected edge from a solver tree.
void remove_edge(SolverTree& tree, int a, int b) {
  auto& an = tree.nodes[a].neighbors;
  an.erase(std::remove(an.begin(), an.end(), b), an.end());
  auto& bn = tree.nodes[b].neighbors;
  bn.erase(std::remove(bn.begin(), bn.end(), a), bn.end());
}

/// Returns active neighbors of a solver-tree node.
std::vector<int> active_neighbors(const SolverTree& tree, int node) {
  std::vector<int> result;
  for (int nb : tree.nodes[node].neighbors) {
    if (tree.nodes[nb].active) result.push_back(nb);
  }
  return result;
}

/// Removes inactive degree-two internal nodes from the final solution tree.
void suppress_solver_degree_two(SolverTree& tree) {
  bool changed = true;
  while (changed) {
    changed = false;
    for (int i = 0; i < static_cast<int>(tree.nodes.size()); ++i) {
      if (!tree.nodes[i].active || tree.nodes[i].identity_id >= 0) continue;
      const std::vector<int> nbs = active_neighbors(tree, i);
      if (nbs.size() != 2) continue;
      remove_edge(tree, i, nbs[0]);
      remove_edge(tree, i, nbs[1]);
      tree.nodes[i].active = false;
      tree.nodes[i].neighbors.clear();
      add_edge(tree, nbs[0], nbs[1]);
      changed = true;
    }
  }
}

/// Converts an internal solver tree to the public species tree.
SpeciesTree to_species(const SolverTree& solved, const TaxonRegistry& taxa) {
  SolverTree compact = solved;
  suppress_solver_degree_two(compact);
  SpeciesTree species;
  std::vector<int> map(compact.nodes.size(), -1);
  for (int i = 0; i < static_cast<int>(compact.nodes.size()); ++i) {
    if (!compact.nodes[i].active) continue;
    if (compact.nodes[i].identity_id >= static_cast<int64_t>(taxa.size())) continue;
    map[i] = static_cast<int>(species.nodes.size());
    SpeciesTree::Node node;
    if (compact.nodes[i].identity_id >= 0) {
      node.label = taxa.label_for(static_cast<int>(compact.nodes[i].identity_id));
    }
    species.nodes.push_back(std::move(node));
  }
  for (int i = 0; i < static_cast<int>(compact.nodes.size()); ++i) {
    if (map[i] < 0) continue;
    for (int nb : compact.nodes[i].neighbors) {
      if (map[nb] >= 0) species.nodes[map[i]].neighbors.push_back(map[nb]);
    }
  }
  return species;
}

/// Writes benchmark counters and options as JSON.
void write_benchmark_json(const Options& options, const SolverContext& ctx) {
  if (options.benchmark_json.empty()) return;
  std::ofstream out(options.benchmark_json);
  if (!out) return;
  const BenchmarkStats& stats = benchmark_stats();
  out << "{\n";
  out << "  \"algorithm\": \"" << algorithm_name(options.algorithm) << "\",\n";
  out << "  \"taxon_order\": \"" << taxon_order_name(options.taxon_order) << "\",\n";
  out << "  \"progress\": \"" << progress_mode_name(options.progress) << "\",\n";
  out << "  \"genes\": " << base_genes(ctx).size() << ",\n";
  out << "  \"taxa\": " << context_taxa(ctx).size() << ",\n";
  out << "  \"workers\": " << ctx.workers << ",\n";
  out << "  \"hat_split_mode\": \""
      << hat_split_name(options.hat_split_mode) << "\",\n";
  out << "  \"hat_one_third_singleton\": "
      << (options.hat_one_third_singleton ? "true" : "false") << ",\n";
  out << "  \"hat_split_certificate\": "
      << (options.hat_split_certificate ? "true" : "false") << ",\n";
  out << "  \"hat_certificate_available\": "
      << (ctx.split_certificate_index &&
                  ctx.split_certificate_index->complete_single_copy
              ? "true" : "false") << ",\n";
  out << "  \"hat_certificate_coefficient\": "
      << options.hat_certificate_coefficient << ",\n";
  out << "  \"hat_certificate_crossover\": "
      << options.hat_certificate_crossover << ",\n";
  out << "  \"hat_certificate_range_scale\": "
      << 2.0 * std::max(0, internal::context_taxa(ctx).size() - 1) << ",\n";
  out << "  \"hat_certificate_threshold\": ";
  if (ctx.split_certificate_index) {
    out << internal::hat_split_certificate_threshold(
        static_cast<int>(ctx.split_certificate_index->genes.size()),
        ctx.split_certificate_index->taxon_count,
        options.hat_certificate_coefficient, options.hat_certificate_crossover);
  } else {
    out << "null";
  }
  out << ",\n";
  out << "  \"hat_certificate_candidates\": "
      << options.hat_certificate_candidates << ",\n";
  out << "  \"hat_certificate_anchors\": "
      << options.hat_certificate_anchors << ",\n";
  out << "  \"hat_one_x_mode\": \""
      << hat_one_x_name(options.hat_one_x_mode) << "\",\n";
  out << "  \"hat_norm_alpha\": " << options.hat_norm_alpha << ",\n";
  out << "  \"hat_nni_passes\": " << options.hat_nni_passes << ",\n";
  out << "  \"hat_nni_passes_explicit\": "
      << (options.hat_nni_passes_explicit ? "true" : "false") << ",\n";
  out << "  \"hat_nni_unbounded\": "
      << (options.hat_nni_unbounded ? "true" : "false") << ",\n";
  out << "  \"hat_nni_unbounded_explicit\": "
      << (options.hat_nni_unbounded_explicit ? "true" : "false") << ",\n";
  out << "  \"hat_nni_local_rescore\": "
      << (options.hat_nni_local_rescore ? "true" : "false") << ",\n";
  out << "  \"hat_nni_min_gain_per_gene\": "
      << options.hat_nni_min_gain_per_gene << ",\n";
  out << "  \"hat_nni_prefix_gate\": "
      << (options.hat_nni_prefix_gate ? "true" : "false") << ",\n";
  out << "  \"hat_nni_prefix_gate_explicit\": "
      << (options.hat_nni_prefix_gate_explicit ? "true" : "false") << ",\n";
  out << "  \"hat_nni_gate_low_first_gain_per_gene\": "
      << options.hat_nni_gate_low_first_gain_per_gene << ",\n";
  out << "  \"hat_nni_gate_high_first_gain_per_gene\": "
      << options.hat_nni_gate_high_first_gain_per_gene << ",\n";
  out << "  \"hat_nni_gate_high_gain_pass_cap\": "
      << options.hat_nni_gate_high_gain_pass_cap << ",\n";
  out << "  \"hat_nni_gate_middle_min_gain_per_gene\": "
      << options.hat_nni_gate_middle_min_gain_per_gene << ",\n";
  out << "  \"hat_nni_gate_middle_rel_max_gain\": "
      << options.hat_nni_gate_middle_rel_max_gain << ",\n";
  out << "  \"hat_nni_gate_middle_min_passes\": "
      << options.hat_nni_gate_middle_min_passes << ",\n";
  out << "  \"hat_nni_gate_max_passes\": "
      << options.hat_nni_gate_max_passes << ",\n";
  out << "  \"hat_nni_gate_split_drop_v1\": "
      << options.hat_nni_gate_split_drop_v1 << ",\n";
  out << "  \"hat_nni_gate_split_drop_per_gene_v1\": "
      << options.hat_nni_gate_split_drop_per_gene_v1 << ",\n";
  out << "  \"hat_nni_pass_tree_dir\": \""
      << json_escape(options.hat_nni_pass_tree_dir) << "\",\n";
  out << "  \"hat_structural_budget_factor\": "
      << options.hat_structural_budget_factor << ",\n";
  auto write_counter = [&](const char* name, const std::atomic<long long>& c) {
    out << "  \"" << name << "\": "
        << c.load(std::memory_order_relaxed) << ",\n";
  };
  auto write_double_counter = [&](const char* name,
                                  const std::atomic<double>& c) {
    out << "  \"" << name << "\": "
        << c.load(std::memory_order_relaxed) << ",\n";
  };

  write_counter("consensus_distinct_splits", stats.consensus_distinct_splits);
  write_counter("consensus_accepted_splits", stats.consensus_accepted_splits);
  write_counter("consensus_hash_low", stats.consensus_hash_low);
  write_counter("candidate_sets", stats.candidate_sets);
  write_counter("candidates", stats.candidates);
  write_counter("max_candidates", stats.max_candidates);
  write_counter("active_consensus_candidate_sets",
                stats.active_consensus_candidate_sets);
  write_counter("active_consensus_candidates",
                stats.active_consensus_candidates);
  write_counter("active_gene_split_candidates_seen",
                stats.active_gene_split_candidates_seen);
  write_counter("active_gene_split_candidates_added",
                stats.active_gene_split_candidates_added);
  write_counter("active_gene_split_candidate_duplicates",
                stats.active_gene_split_candidate_duplicates);
  write_counter("active_total_candidates_after_expansion",
                stats.active_total_candidates_after_expansion);
  write_counter("attachment_calls", stats.attachment_calls);
  write_counter("attachment_internal_nodes", stats.attachment_internal_nodes);
  write_counter("attachment_four_color_passes",
                stats.attachment_four_color_passes);
  write_counter("attachment_edge_ties", stats.attachment_edge_ties);
  write_counter("attachment_normalized_one_x_calls",
                stats.attachment_normalized_one_x_calls);
  write_counter("hat_subproblems", stats.hat_subproblems);
  write_counter("hat_max_active_taxa", stats.hat_max_active_taxa);
  write_counter("hat_consensus_index_builds",
                stats.hat_consensus_index_builds);
  write_counter("hat_consensus_index_nodes", stats.hat_consensus_index_nodes);
  write_counter("hat_active_local_taxon_entries",
                stats.hat_active_local_taxon_entries);
  write_counter("hat_active_local_index_entries",
                stats.hat_active_local_index_entries);
  write_counter("hat_virtual_query_selected_leaves",
                stats.hat_virtual_query_selected_leaves);
  write_counter("hat_virtual_query_lca_queries",
                stats.hat_virtual_query_lca_queries);
  write_counter("hat_virtual_query_retained_nodes",
                stats.hat_virtual_query_retained_nodes);
  write_counter("hat_localized_scoring_gene_nodes",
                stats.hat_localized_scoring_gene_nodes);
  write_counter("hat_localized_scoring_calls",
                stats.hat_localized_scoring_calls);
  write_counter("hat_outer_recursive_calls", stats.hat_outer_recursive_calls);
  write_counter("hat_outer_nonbase_calls", stats.hat_outer_nonbase_calls);
  write_counter("hat_outer_max_depth", stats.hat_outer_max_depth);
  write_counter("hat_outer_active_identity_mass",
                stats.hat_outer_active_identity_mass);
  write_counter("hat_outer_active_support_word_mass",
                stats.hat_outer_active_support_word_mass);
  write_counter("hat_outer_active_taxon_mass",
                stats.hat_outer_active_taxon_mass);
  write_counter("hat_outer_split_identity_tests",
                stats.hat_outer_split_identity_tests);
  write_counter("hat_outer_partition_identity_moves",
                stats.hat_outer_partition_identity_moves);
  write_counter("hat_outer_partition_support_word_moves",
                stats.hat_outer_partition_support_word_moves);
  write_counter("hat_outer_partition_taxon_moves",
                stats.hat_outer_partition_taxon_moves);
  write_counter("hat_outer_balanced_splits",
                stats.hat_outer_balanced_splits);
  write_counter("hat_outer_maximally_skewed_splits",
                stats.hat_outer_maximally_skewed_splits);
  write_counter("hat_outer_catalog_identity_entries",
                stats.hat_outer_catalog_identity_entries);
  write_counter("hat_outer_catalog_taxon_entries",
                stats.hat_outer_catalog_taxon_entries);
  write_counter("hat_outer_catalog_lookup_entries",
                stats.hat_outer_catalog_lookup_entries);
  write_counter("hat_outer_tree_node_visits",
                stats.hat_outer_tree_node_visits);
  write_counter("hat_outer_tree_identity_lookups",
                stats.hat_outer_tree_identity_lookups);
  write_counter("hat_outer_tree_taxon_entries",
                stats.hat_outer_tree_taxon_entries);
  write_counter("hat_outer_taxon_radix_entries",
                stats.hat_outer_taxon_radix_entries);
  write_counter("hat_outer_taxon_radix_bucket_ops",
                stats.hat_outer_taxon_radix_bucket_ops);
  write_counter("hat_outer_taxon_unique_entries",
                stats.hat_outer_taxon_unique_entries);
  write_counter("hat_outer_root_bank_nodes",
                stats.hat_outer_root_bank_nodes);
  write_counter("hat_outer_root_bank_other_entries",
                stats.hat_outer_root_bank_other_entries);
  write_counter("hat_outer_live_payload_units",
                stats.hat_outer_live_payload_units);
  write_counter("hat_outer_peak_payload_units",
                stats.hat_outer_peak_payload_units);
  write_counter("hat_outer_live_identity_entries",
                stats.hat_outer_live_identity_entries);
  write_counter("hat_outer_peak_identity_entries",
                stats.hat_outer_peak_identity_entries);
  write_counter("hat_outer_live_support_words",
                stats.hat_outer_live_support_words);
  write_counter("hat_outer_peak_support_words",
                stats.hat_outer_peak_support_words);
  write_counter("hat_outer_live_identity_taxon_entries",
                stats.hat_outer_live_identity_taxon_entries);
  write_counter("hat_outer_peak_identity_taxon_entries",
                stats.hat_outer_peak_identity_taxon_entries);
  write_counter("hat_outer_live_identity_key_bytes",
                stats.hat_outer_live_identity_key_bytes);
  write_counter("hat_outer_peak_identity_key_bytes",
                stats.hat_outer_peak_identity_key_bytes);
  write_counter("hat_outer_live_catalog_entries",
                stats.hat_outer_live_catalog_entries);
  write_counter("hat_outer_peak_catalog_entries",
                stats.hat_outer_peak_catalog_entries);
  write_counter("hat_outer_live_bank_nodes",
                stats.hat_outer_live_bank_nodes);
  write_counter("hat_outer_peak_bank_nodes",
                stats.hat_outer_peak_bank_nodes);
  write_counter("hat_outer_live_bank_other_entries",
                stats.hat_outer_live_bank_other_entries);
  write_counter("hat_outer_peak_bank_other_entries",
                stats.hat_outer_peak_bank_other_entries);
  write_counter("hat_outer_live_bank_identity_entries",
                stats.hat_outer_live_bank_identity_entries);
  write_counter("hat_outer_peak_bank_identity_entries",
                stats.hat_outer_peak_bank_identity_entries);
  write_counter("hat_outer_live_continuation_identity_entries",
                stats.hat_outer_live_continuation_identity_entries);
  write_counter("hat_outer_peak_continuation_identity_entries",
                stats.hat_outer_peak_continuation_identity_entries);
  write_counter("hat_outer_live_continuation_support_words",
                stats.hat_outer_live_continuation_support_words);
  write_counter("hat_outer_peak_continuation_support_words",
                stats.hat_outer_peak_continuation_support_words);
  write_counter("hat_outer_live_continuation_taxon_entries",
                stats.hat_outer_live_continuation_taxon_entries);
  write_counter("hat_outer_peak_continuation_taxon_entries",
                stats.hat_outer_peak_continuation_taxon_entries);
  write_counter("hat_outer_live_continuation_key_bytes",
                stats.hat_outer_live_continuation_key_bytes);
  write_counter("hat_outer_peak_continuation_key_bytes",
                stats.hat_outer_peak_continuation_key_bytes);
  write_counter("hat_outer_live_completed_tree_nodes",
                stats.hat_outer_live_completed_tree_nodes);
  write_counter("hat_outer_peak_completed_tree_nodes",
                stats.hat_outer_peak_completed_tree_nodes);
  write_counter("hat_outer_live_completed_tree_adjacencies",
                stats.hat_outer_live_completed_tree_adjacencies);
  write_counter("hat_outer_peak_completed_tree_adjacencies",
                stats.hat_outer_peak_completed_tree_adjacencies);
  write_counter("hat_outer_live_taxon_scratch_entries",
                stats.hat_outer_live_taxon_scratch_entries);
  write_counter("hat_outer_peak_taxon_scratch_entries",
                stats.hat_outer_peak_taxon_scratch_entries);
  write_counter("hat_outer_live_accounting_violations",
                stats.hat_outer_live_accounting_violations);
  write_counter("hat_consensus_edges_scored", stats.hat_consensus_edges_scored);
  write_counter("hat_hdt_banks_built", stats.hat_hdt_banks_built);
  write_counter("hat_hdt_bank_build_ns", stats.hat_hdt_bank_build_ns);
  write_counter("hat_hdt_recolor_ops", stats.hat_hdt_recolor_ops);
  write_counter("hat_hdt_nodes_recomputed", stats.hat_hdt_nodes_recomputed);
  write_counter("hat_hdt_queries", stats.hat_hdt_queries);
  write_counter("hat_reduced_hdt_nodes", stats.hat_reduced_hdt_nodes);
  write_counter("hat_contract_calls", stats.hat_contract_calls);
  write_counter("hat_extract_calls", stats.hat_extract_calls);
  write_counter("hat_contract_input_nodes", stats.hat_contract_input_nodes);
  write_counter("hat_contract_output_nodes", stats.hat_contract_output_nodes);
  write_counter("hat_contract_live_leaves", stats.hat_contract_live_leaves);
  write_counter("hat_contract_max_output_per_live_ppm",
                stats.hat_contract_max_output_per_live_ppm);
  write_counter("hat_extract_input_nodes", stats.hat_extract_input_nodes);
  write_counter("hat_extract_output_nodes", stats.hat_extract_output_nodes);
  write_counter("hat_extract_live_leaves", stats.hat_extract_live_leaves);
  write_counter("hat_extract_max_output_per_live_ppm",
                stats.hat_extract_max_output_per_live_ppm);
  write_counter("hat_extract_source_hdt_nodes",
                stats.hat_extract_source_hdt_nodes);
  write_counter("hat_extract_marked_hdt_nodes",
                stats.hat_extract_marked_hdt_nodes);
  write_counter("hat_extract_hanging_components",
                stats.hat_extract_hanging_components);
  write_counter("hat_extract_intermediate_units",
                stats.hat_extract_intermediate_units);
  write_counter("hat_extract_live_leaves_exact",
                stats.hat_extract_live_leaves_exact);
  write_counter("hat_contract_light_calls", stats.hat_contract_light_calls);
  write_counter("hat_contract_heavy_calls", stats.hat_contract_heavy_calls);
  write_counter("hat_contract_input_units", stats.hat_contract_input_units);
  write_counter("hat_contract_output_units", stats.hat_contract_output_units);
  write_counter("hat_contract_live_leaves_exact",
                stats.hat_contract_live_leaves_exact);
  write_counter("hat_contract_corridors_emitted",
                stats.hat_contract_corridors_emitted);
  write_counter("hat_contract_background_units_absorbed",
                stats.hat_contract_background_units_absorbed);
  write_counter("hat_contract_polytomy_arms_folded",
                stats.hat_contract_polytomy_arms_folded);
  write_counter("hat_contract_cc_ops", stats.hat_contract_cc_ops);
  write_counter("hat_contract_gg_ops", stats.hat_contract_gg_ops);
  write_counter("hat_contract_ig_ops", stats.hat_contract_ig_ops);
  write_counter("hat_contract_c_to_g_ops", stats.hat_contract_c_to_g_ops);
  write_counter("hat_contract_bound_violations",
                stats.hat_contract_bound_violations);
  write_counter("hat_heavy_input_units", stats.hat_heavy_input_units);
  write_counter("hat_heavy_output_units", stats.hat_heavy_output_units);
  write_counter("hat_heavy_calls", stats.hat_heavy_calls);
  write_counter("hat_heavy_shrink_violations",
                stats.hat_heavy_shrink_violations);
  write_counter("hat_hdt_initial_atoms", stats.hat_hdt_initial_atoms);
  write_counter("hat_hdt_rebuild_nodes", stats.hat_hdt_rebuild_nodes);
  write_counter("hat_hdt_local_balance_violations",
                stats.hat_hdt_local_balance_violations);
  write_counter("hat_hdt_max_update_height",
                stats.hat_hdt_max_update_height);
  write_counter("hat_recolor_changed_leaves",
                stats.hat_recolor_changed_leaves);
  write_counter("hat_recolor_marked_nodes", stats.hat_recolor_marked_nodes);
  write_counter("hat_dp_transition_evals", stats.hat_dp_transition_evals);
  write_counter("hat_dp_scalar_fmas", stats.hat_dp_scalar_fmas);
  write_counter("hat_extract_scratch_peak_bytes",
                stats.hat_extract_scratch_peak_bytes);
  write_counter("hat_contract_scratch_peak_bytes",
                stats.hat_contract_scratch_peak_bytes);
  write_counter("hat_state_arena_peak_bytes",
                stats.hat_state_arena_peak_bytes);
  write_counter("hat_worker_total_peak_bytes",
                stats.hat_worker_total_peak_bytes);
  write_counter("hat_score_partial_peak_bytes",
                stats.hat_score_partial_peak_bytes);
  write_counter("hat_contract_mark_ns", stats.hat_contract_mark_ns);
  write_counter("hat_contract_background_fold_ns",
                stats.hat_contract_background_fold_ns);
  write_counter("hat_contract_corridor_fold_ns",
                stats.hat_contract_corridor_fold_ns);
  write_counter("hat_contract_emit_ns", stats.hat_contract_emit_ns);
  write_counter("hat_hdt_rebuild_ns", stats.hat_hdt_rebuild_ns);
  write_counter("complexity_contract_satisfied",
                stats.complexity_contract_satisfied);
  write_counter("hat_frozen_segments", stats.hat_frozen_segments);
  write_counter("hat_frozen_segment_source_degree_sum",
                stats.hat_frozen_segment_source_degree_sum);
  write_counter("hat_recompute_child_folds", stats.hat_recompute_child_folds);
  write_counter("hat_max_recompute_degree", stats.hat_max_recompute_degree);
  write_counter("hat_split_score_ns", stats.hat_split_score_ns);
  write_counter("hat_certificate_index_ns", stats.hat_certificate_index_ns);
  write_counter("hat_certificate_ns", stats.hat_certificate_ns);
  write_counter("hat_certificate_tests", stats.hat_certificate_tests);
  write_counter("hat_certificate_accepts", stats.hat_certificate_accepts);
  write_counter("hat_certificate_rejects", stats.hat_certificate_rejects);
  write_counter("hat_certificate_negative_accepts", stats.hat_certificate_negative_accepts);
  write_counter("hat_certificate_max_negative_gap_twice_k",
                stats.hat_certificate_max_negative_gap_twice_k);
  write_counter("hat_certificate_singletons", stats.hat_certificate_singletons);
  write_counter("hat_certificate_distance_queries",
                stats.hat_certificate_distance_queries);
  write_counter("hat_one_x_calls", stats.hat_one_x_calls);
  write_counter("hat_one_x_medians", stats.hat_one_x_medians);
  write_counter("hat_one_x_ns", stats.hat_one_x_ns);
  write_counter("hat_root_source_scoring_calls",
                stats.hat_root_source_scoring_calls);
  write_counter("hat_child_gene_bank_reuses", stats.hat_child_gene_bank_reuses);
  write_counter("hat_child_bank_state_attempts",
                stats.hat_child_bank_state_attempts);
  write_counter("hat_child_bank_state_accepts",
                stats.hat_child_bank_state_accepts);
  write_counter("hat_child_bank_state_heavy_reuses",
                stats.hat_child_bank_state_heavy_reuses);
  write_counter("hat_child_bank_state_pair_violations",
                stats.hat_child_bank_state_pair_violations);
  write_counter("hat_child_bank_state_input_nodes",
                stats.hat_child_bank_state_input_nodes);
  write_counter("hat_child_bank_state_output_nodes",
                stats.hat_child_bank_state_output_nodes);
  write_counter("hat_child_bank_state_max_accepted_ratio_ppm",
                stats.hat_child_bank_state_max_accepted_ratio_ppm);
  write_counter("hat_max_scoring_source_identities",
                stats.hat_max_scoring_source_identities);
  write_counter("hat_max_scoring_source_inflation_ppm",
                stats.hat_max_scoring_source_inflation_ppm);
  write_counter("hat_scoring_source_inflation_violations",
                stats.hat_scoring_source_inflation_violations);
  write_counter("hat_structural_serial_child_subproblems",
                stats.hat_structural_serial_child_subproblems);
  write_counter("hat_structural_parallel_child_subproblems",
                stats.hat_structural_parallel_child_subproblems);
  write_counter("hat_structural_parallel_rejections",
                stats.hat_structural_parallel_rejections);
  write_counter("hat_structural_parallel_physical_rejections",
                stats.hat_structural_parallel_physical_rejections);
  write_counter("hat_active_structural_parallel_slots",
                stats.hat_active_structural_parallel_slots);
  write_counter("hat_max_structural_parallel_slots",
                stats.hat_max_structural_parallel_slots);
  write_counter("hat_nni_passes_executed", stats.hat_nni_passes_executed);
  write_counter("hat_nni_edges_scored", stats.hat_nni_edges_scored);
  write_counter("hat_nni_full_score_passes", stats.hat_nni_full_score_passes);
  write_counter("hat_nni_local_score_passes", stats.hat_nni_local_score_passes);
  write_counter("hat_nni_local_edges_scored", stats.hat_nni_local_edges_scored);
  write_counter("hat_nni_moves_accepted", stats.hat_nni_moves_accepted);
  write_counter("hat_nni_moves_rejected_by_threshold",
                stats.hat_nni_moves_rejected_by_threshold);
  write_counter("hat_nni_peak_proposals", stats.hat_nni_peak_proposals);
  write_counter("hat_nni_max_qsplit_leaves", stats.hat_nni_max_qsplit_leaves);
  write_counter("hat_nni_qsplit_owner_bytes_avoided",
                stats.hat_nni_qsplit_owner_bytes_avoided);
  write_double_counter("hat_nni_gain_threshold", stats.hat_nni_gain_threshold);
  write_double_counter("hat_nni_last_accepted_gain",
                       stats.hat_nni_last_accepted_gain);
  write_double_counter("hat_nni_max_accepted_gain",
                       stats.hat_nni_max_accepted_gain);
  write_double_counter("hat_nni_verified_gain", stats.hat_nni_verified_gain);
  write_counter("hat_nni_prefix_gate_stops", stats.hat_nni_prefix_gate_stops);
  write_counter("hat_nni_prefix_gate_quartet_stops",
                stats.hat_nni_prefix_gate_quartet_stops);
  write_counter("hat_nni_prefix_gate_split_stops",
                stats.hat_nni_prefix_gate_split_stops);
  write_counter("hat_nni_split_support_delta",
                stats.hat_nni_split_support_delta);
  write_counter("hat_nni_split_support_unique_splits",
                stats.hat_nni_split_support_unique_splits);
  write_counter("hat_nni_split_support_ns", stats.hat_nni_split_support_ns);
  write_counter("hat_nni_ns", stats.hat_nni_ns);
  write_counter("mre_candidates_checked", stats.mre_candidates_checked);
  write_counter("mre_candidates_rejected", stats.mre_candidates_rejected);
  write_counter("mre_case_a_accepts", stats.mre_case_a_accepts);
  write_counter("mre_case_b_accepts", stats.mre_case_b_accepts);
  write_counter("mre_moved_children", stats.mre_moved_children);
  write_counter("mre_max_child_degree", stats.mre_max_child_degree);
  write_counter("mre_radix_passes", stats.mre_radix_passes);
  write_counter("mre_radix_items", stats.mre_radix_items);
  write_counter("mre_child_index_operations",
                stats.mre_child_index_operations);
  write_counter("mre_child_index_rotations",
                stats.mre_child_index_rotations);
  write_counter("mre_child_index_max_height",
                stats.mre_child_index_max_height);
  write_counter("mre_child_index_height_violations",
                stats.mre_child_index_height_violations);
  write_counter("mre_child_index_enumeration_visits",
                stats.mre_child_index_enumeration_visits);
  write_counter("mre_topology_slots_allocated",
                stats.mre_topology_slots_allocated);
  write_counter("mre_topology_peak_live_slots",
                stats.mre_topology_peak_live_slots);
  write_counter("hat_peak_bank_bytes", stats.hat_peak_bank_bytes);
  write_counter("subproblems", stats.subproblems);
  write_counter("completed_subproblems", stats.completed_subproblems);
  write_counter("max_depth", stats.max_depth);
  write_counter("max_active_identities", stats.max_active_identities);
  write_counter("max_active_taxa", stats.max_active_taxa);
  write_counter("worst_balance_ratio_ppm", stats.worst_balance_ratio_ppm);
  write_counter("choose_initial_ns", stats.choose_initial_ns);
  write_counter("attachment_range_add_ns", stats.attachment_range_add_ns);
  write_counter("attachment_score_ns", stats.attachment_score_ns);
  write_counter("attachment_total_ns", stats.attachment_total_ns);
  write_counter("parse_input_ns", stats.parse_input_ns);
  write_counter("taxon_registry_ns", stats.taxon_registry_ns);
  write_counter("taxon_registry_physical_labels",
                stats.taxon_registry_physical_labels);
  write_counter("taxon_registry_unique_labels",
                stats.taxon_registry_unique_labels);
  write_counter("taxon_registry_sort_items",
                stats.taxon_registry_sort_items);
  write_counter("gene_index_ns", stats.gene_index_ns);
  write_counter("consensus_count_ns", stats.consensus_count_ns);
  write_counter("consensus_assembly_ns", stats.consensus_assembly_ns);
  write_counter("startup_ns", stats.startup_ns);
  write_counter("compression_materializations",
                stats.compression_materializations);
  write_counter("compression_materialization_ns",
                stats.compression_materialization_ns);
  write_counter("active_compressed_banks", stats.active_compressed_banks);
  write_counter("max_active_compressed_banks",
                stats.max_active_compressed_banks);
  write_counter("active_child_compressed_banks",
                stats.active_child_compressed_banks);
  write_counter("max_active_child_compressed_banks",
                stats.max_active_child_compressed_banks);
  write_counter("active_child_compressed_bank_identity_slots",
                stats.active_child_compressed_bank_identity_slots);
  write_counter("max_active_child_compressed_bank_identity_slots",
                stats.max_active_child_compressed_bank_identity_slots);
  write_counter("max_child_compressed_bank_bytes",
                stats.max_child_compressed_bank_bytes);
  write_counter("cpu_work_units", stats.cpu_work_units);

  write_benchmark_subproblem_records_json(out);
  out << "\n";
  out << "}\n";
}

#undef WQFM_ALWAYS_INLINE

}  // namespace

}  // namespace wqfm::internal

namespace wqfm {

/// Runs the full species tree inference pipeline.
SpeciesTree infer_species_tree(std::vector<GeneTree> genes, const Options& options) {
  if (genes.empty()) {
    throw std::invalid_argument("at least one gene tree is required");
  }
  if (genes.size() >
      static_cast<std::size_t>(std::numeric_limits<int>::max())) {
    throw std::length_error("gene count exceeds the integer-indexed domain");
  }
  if (!std::isfinite(options.hat_structural_budget_factor) ||
      options.hat_structural_budget_factor < 1.0) {
    throw std::invalid_argument(
        "HAT structural budget factor must be finite and at least one");
  }
  if (!std::isfinite(options.hat_norm_alpha)) {
    throw std::invalid_argument("HAT normalization alpha must be finite");
  }
  if (!std::isfinite(options.hat_certificate_coefficient) ||
      options.hat_certificate_coefficient <= 0.0 ||
      !std::isfinite(options.hat_certificate_crossover) ||
      options.hat_certificate_crossover < 0.0 ||
      options.hat_certificate_candidates < 1 ||
      options.hat_certificate_candidates > 16 ||
      options.hat_certificate_anchors < 1 ||
      options.hat_certificate_anchors > 4) {
    throw std::invalid_argument(
        "HAT certificate requires a finite positive coefficient, finite "
        "nonnegative crossover, 1..16 candidates, and 1..4 anchor pairs");
  }
  if (options.hat_one_third_singleton) {
    if (options.algorithm != AlgorithmMode::HatV1) {
      throw std::invalid_argument(
          "HAT one-third/singleton mode requires --algorithm hat-v1");
    }
    if (!options.consensus_path.empty()) {
      throw std::invalid_argument(
          "HAT one-third/singleton mode does not accept an external consensus");
    }
    if (options.hat_split_mode != HatSplitMode::Raw &&
        options.hat_split_mode != HatSplitMode::ConsensusBalanced) {
      throw std::invalid_argument(
          "HAT one-third/singleton mode requires --hat-split raw or "
          "consensus-balanced");
    }
  }
  std::size_t physical_nodes = 0;
  for (const GeneTree& gene : genes) {
    physical_nodes = std::min(
        static_cast<std::size_t>(std::numeric_limits<int>::max()),
        physical_nodes + std::min(
                             gene.nodes.size(),
                             static_cast<std::size_t>(
                                 std::numeric_limits<int>::max()) -
                                 physical_nodes));
  }
  const int worker_cap = std::max(1, static_cast<int>(physical_nodes));
  internal::ScopedInferenceRun inference_run(options, worker_cap);
  internal::ScopedInferenceProgressFinish progress_finish(
      options.progress == ProgressMode::On);
  internal::DebugScope scope("api.infer_species_tree");
  const auto startup_started = std::chrono::steady_clock::now();
  if (inference_run.owns_process_telemetry()) {
    internal::reset_benchmark_stats();
    internal::benchmark_stats().parse_input_ns.store(
        options.parse_input_ns, std::memory_order_relaxed);
  }
  if (!options.debug_trace.empty()) {
    std::ofstream clear_trace(options.debug_trace, std::ios::trunc);
  }
  internal::SolverContext ctx;
  ctx.options = options;
  if (ctx.options.algorithm == AlgorithmMode::HatV1 ||
      ctx.options.algorithm == AlgorithmMode::HatConsensus ||
      ctx.options.algorithm == AlgorithmMode::HatConsensusNniV1) {
    const bool runs_hat_nni =
        ctx.options.algorithm == AlgorithmMode::HatV1 ||
        ctx.options.algorithm == AlgorithmMode::HatConsensusNniV1;
    if (runs_hat_nni && !ctx.options.hat_nni_passes_explicit &&
        !ctx.options.hat_nni_unbounded_explicit) {
      if (!ctx.options.hat_nni_prefix_gate_explicit) {
        ctx.options.hat_nni_prefix_gate = true;
      }
      if (ctx.options.hat_nni_prefix_gate) {
        ctx.options.hat_nni_unbounded = false;
        ctx.options.hat_nni_passes =
            std::max(ctx.options.hat_nni_passes,
                     ctx.options.hat_nni_gate_max_passes);
      }
    }
  }
  const Options& run_options = ctx.options;
  ctx.workers = inference_run.workers();
  const std::string startup_settings =
      internal::startup_settings_summary(run_options);
  internal::progress_startup_status(
      run_options, "taxa", startup_settings + " building registry");
  {
    const auto stage_started = std::chrono::steady_clock::now();
    std::vector<std::string> labels;
    for (const GeneTree& gene : genes) {
      auto leaf_labels = gene.leaf_labels();
      labels.insert(labels.end(), leaf_labels.begin(), leaf_labels.end());
    }
    if (labels.size() > static_cast<std::size_t>(
                            std::numeric_limits<long long>::max())) {
      throw std::length_error("physical taxon-label count exceeds telemetry domain");
    }
    const long long physical_labels = static_cast<long long>(labels.size());
    ctx.taxa = internal::build_registry_with_order(std::move(labels),
                                                   run_options.taxon_order);
    internal::benchmark_stats().taxon_registry_physical_labels.store(
        physical_labels, std::memory_order_relaxed);
    internal::benchmark_stats().taxon_registry_unique_labels.store(
        ctx.taxa.size(), std::memory_order_relaxed);
    internal::benchmark_stats().taxon_registry_sort_items.store(
        ctx.taxa.size(), std::memory_order_relaxed);
    internal::benchmark_timer_add(
        internal::benchmark_stats().taxon_registry_ns,
        std::chrono::duration_cast<std::chrono::nanoseconds>(
            std::chrono::steady_clock::now() - stage_started));
  }
  internal::progress_startup_status(
      run_options, "index", startup_settings + " indexing gene trees");
  ctx.genes.resize(genes.size());
  const auto index_started = std::chrono::steady_clock::now();
  internal::parallel_for_inner(
      0, static_cast<int>(genes.size()),
      [&](int i) {
        ctx.genes[static_cast<std::size_t>(i)] =
            internal::index_gene(genes[static_cast<std::size_t>(i)], ctx.taxa);
      },
      1);
  std::vector<GeneTree>().swap(genes);
  internal::benchmark_timer_add(
      internal::benchmark_stats().gene_index_ns,
      std::chrono::duration_cast<std::chrono::nanoseconds>(
          std::chrono::steady_clock::now() - index_started));

  if (run_options.algorithm == AlgorithmMode::HatV1 &&
      run_options.hat_split_certificate && !run_options.hat_one_third_singleton) {
    ctx.split_certificate_index = internal::make_hat_split_certificate_index(
        ctx.genes, ctx.taxa.size());
    if (!ctx.split_certificate_index->complete_single_copy) {
      ctx.split_certificate_index.reset();
    }
  }

  internal::progress_startup_status(
      run_options, "consensus",
      startup_settings +
          (run_options.hat_one_third_singleton
               ? " building strict one-third greedy consensus"
               : " building MRE"));
  if (run_options.consensus_path.empty()) {
    ctx.consensus =
        internal::build_consensus(ctx.genes, ctx.taxa.size(), run_options);
  } else {
    std::ifstream consensus_in(run_options.consensus_path);
    if (!consensus_in) {
      throw std::runtime_error("failed to open consensus file: " +
                               run_options.consensus_path);
    }
    const std::string consensus_text(
        (std::istreambuf_iterator<char>(consensus_in)),
        std::istreambuf_iterator<char>());
    ctx.consensus =
        internal::consensus_from_tree(parse_newick(consensus_text), ctx.taxa);
  }
  internal::write_consensus_dump(run_options, ctx);
  std::vector<internal::Identity> identities = internal::make_initial_identities(ctx.taxa);
  internal::benchmark_timer_add(
      internal::benchmark_stats().startup_ns,
      std::chrono::duration_cast<std::chrono::nanoseconds>(
          std::chrono::steady_clock::now() - startup_started));
  internal::progress_start(run_options, std::max(0, ctx.taxa.size() - 3));

  internal::SolverTree solved;
  tf::Taskflow solve_flow;
  solve_flow.emplace([&](tf::Subflow& sf) {
    switch (run_options.algorithm) {
      case AlgorithmMode::HatV1:
        // Shared quiet runs keep process-global structural budget gauges out of
        // their scheduling decisions. Scoring loops still use the executor.
        solved = internal::solve_hat_v1(ctx, std::move(identities),
                                        0x6a09e667f3bcc909ULL, 0,
                                        inference_run.permits_recursive_subflows()
                                            ? &sf
                                            : nullptr);
        break;
      case AlgorithmMode::HatConsensus:
        solved = internal::solve_hat_consensus(ctx, false);
        break;
      case AlgorithmMode::HatConsensusNniV1:
        solved = internal::solve_hat_consensus(ctx, true);
        break;
    }
  });
  internal::get_executor().run(solve_flow).get();
  progress_finish.finish();
  internal::write_benchmark_json(run_options, ctx);
  internal::debug_report();
  SpeciesTree species = internal::to_species(solved, internal::context_taxa(ctx));
  return species;
}

}  // namespace wqfm
