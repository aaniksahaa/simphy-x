#pragma once

#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cstdint>
#include <exception>
#include <iosfwd>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "wqfm/api.hpp"
#include "wqfm/taxon.hpp"

#include <taskflow/algorithm/for_each.hpp>
#include <taskflow/algorithm/reduce.hpp>
#include <taskflow/taskflow.hpp>

#ifndef WQFM_TASKFLOW_DETERMINISTIC_WORKER_LOOKUP
#error "the deterministic Taskflow worker-lookup overlay is required"
#endif

namespace wqfm::internal {

/// Side of a bipartition or FM assignment.
enum class Side : uint8_t { Left = 0, Right = 1 };
/// Compact fixed-size bitset used for taxa, splits, and identity supports.
class Bitset {
 public:
  Bitset() = default;
  /// Creates a zero-filled bitset with the requested number of bits.
  explicit Bitset(int size);

  /// Returns the fixed number of addressable bits.
  int size() const { return size_; }
  /// Returns true when no bit is set.
  bool empty() const;
  /// Sets one bit, throwing if the index is out of range.
  void set(int pos);
  /// Clears one bit, throwing if the index is out of range.
  void reset(int pos);
  /// Tests one bit, returning false for out-of-range indices.
  bool test(int pos) const;
  /// Counts set bits.
  int count() const;
  /// Clears all bits.
  void clear();
  /// Sets every valid bit.
  void fill();
  /// Returns the bitwise complement within the fixed size.
  Bitset complement() const;
  /// Returns indices of all set bits.
  std::vector<int> positions() const;

  /// Returns whether two same-sized bitsets share any set bit.
  bool intersects(const Bitset& other) const;
  /// Returns whether the intersection is empty.
  bool intersection_empty(const Bitset& other) const;
  /// Counts set bits in the intersection.
  int intersection_count(const Bitset& other) const;
  /// Returns the set intersection.
  Bitset intersection(const Bitset& other) const;
  /// Returns the set union.
  Bitset unite(const Bitset& other) const;
  /// Returns bits present in this set but absent from another.
  Bitset difference(const Bitset& other) const;

  /// Compares bitsets by size and words.
  bool operator==(const Bitset& other) const;
  /// Returns true when bitsets differ.
  bool operator!=(const Bitset& other) const { return !(*this == other); }
  /// Provides deterministic lexicographic ordering.
  bool lex_less(const Bitset& other) const;
  /// Computes a stable hash over size and words.
  std::size_t hash() const;
  /// Exposes raw words for hot unchecked bit tests.
  const std::vector<uint64_t>& words() const { return words_; }

 private:
  int size_ = 0;
  std::vector<uint64_t> words_;

  /// Clears padding bits beyond the fixed size.
  void mask_unused_bits();
};

/// Hash functor for unordered containers keyed by Bitset.
struct BitsetHash {
  /// Returns the stable Bitset hash value.
  std::size_t operator()(const Bitset& bits) const { return bits.hash(); }
};

/// Active taxon unit in a recursive subproblem.
struct Identity {
  int64_t id = -1;
  int universe_size = 0;
  Bitset support;
  std::vector<int> support_taxa;
  std::string stable_key;
};

/// Taxon-indexed tree used by consensus construction.
struct IndexedTree {
  /// Indexed tree node with adjacency and optional taxon id.
  struct Node {
    std::vector<int> neighbors;
    int taxon = -1;
  };

  std::vector<Node> nodes;
  int root = -1;
};

/// Parent/order/interval representation of a rooted tree traversal.
struct RootedTraversal {
  std::vector<int> parent;
  std::vector<int> preorder;
  std::vector<int> postorder;
  std::vector<int> tin;
  std::vector<int> tout;
};

/// Immutable traversal and linear-space LCA data for a consensus tree.
struct ConsensusIndex {
  std::vector<int> parent;
  std::vector<int> depth;
  std::vector<int> dfs_order;
  std::vector<int> tin;
  std::vector<int> tout;
  std::vector<int> taxon_to_leaf;
  std::vector<int> first_occurrence;
  std::vector<int> euler;
  std::vector<int> rmq_nodes;
  int rmq_base = 0;
};

/// Shared one-time publication state for an immutable consensus index.
struct ConsensusIndexCache {
  /// 0 = empty, 1 = building, 2 = published.
  std::atomic<int> state{0};
  std::shared_ptr<const ConsensusIndex> value;
};

/// Rooted consensus tree used to generate candidate splits.
struct ConsensusTree {
  /// Rooted consensus node with optional taxon label and child list.
  struct Node {
    int parent = -1;
    int taxon = -1;
    std::vector<int> children;
  };

  std::vector<Node> nodes;
  int root = -1;
  std::vector<int> taxon_to_leaf;
  std::vector<int> preorder;
  std::vector<int> tin;
  std::vector<int> tout;
  // Consensus trees are immutable once published. Copies share this cache, so
  // recursive child contexts that retain consensus_ref reuse the same index.
  std::shared_ptr<ConsensusIndexCache> index_cache =
      std::make_shared<ConsensusIndexCache>();
};

/// Mutable solution tree assembled by recursive bipartitioning.
struct SolverTree {
  /// Mutable solver node storing adjacency and represented identity id.
  struct Node {
    std::vector<int> neighbors;
    int64_t identity_id = -1;
    bool active = true;
  };

  std::vector<Node> nodes;
};

/// Taxon-to-identity ownership and mass for a compressed subproblem.
struct SubproblemWeights {
  std::vector<int> taxon_to_identity;
  std::vector<double> taxon_weight;
};

/// Rooted traversal order and interval indices for a compressed gene tree.
struct GeneTraversal {
  std::vector<int> parent;
  std::vector<int> preorder;
  std::vector<int> postorder;
  std::vector<int> tin;
  std::vector<int> tout;
};

/// Gene tree representation indexed by taxon ids and optional current identity
/// leaves used by compressed subproblems.
struct GeneIndex {
  /// Indexed gene node storing adjacency plus leaf ownership metadata.
  struct Node {
    std::vector<int> neighbors;
    int taxon = -1;
    int identity = -1;
    double mass = 1.0;
    double real_count = 1.0;
  };

  std::vector<Node> nodes;
  int root = 0;
  int reference_root = 0;
  std::vector<int> taxon_to_node;
  GeneTraversal root_traversal;
  GeneTraversal reference_traversal;
};

/// Per-physical-leaf projection used when an indexed gene contributes taxon
/// splits. Exactly one occurrence represents each distinct taxon in the split
/// key; duplicate_balance is zero only when a subtree contains either all or
/// none of every taxon's physical occurrences, subject to the documented
/// additive-fingerprint noncollision premise.
struct GeneSplitLeafProjection {
  std::vector<uint32_t> unique_size;
  std::vector<std::array<uint64_t, 4>> unique_fingerprint;
  std::vector<std::array<uint64_t, 4>> duplicate_balance;
  std::vector<uint8_t> reference_marker;
  uint32_t distinct_taxa = 0;
  int reference_taxon = -1;
};

/// Candidate bipartition induced by the marker-compressed consensus tree.
struct Candidate {
  Bitset split;
  int representative_order = 0;
};

struct HatSplitCertificateIndex;

/// Shared state for one top-level inference run.
struct SolverContext {
  TaxonRegistry taxa;
  const TaxonRegistry* taxa_ref = nullptr;
  std::vector<GeneIndex> genes;
  const std::vector<GeneIndex>* genes_ref = nullptr;
  const std::vector<GeneIndex>* active_genes = nullptr;
  std::shared_ptr<const std::vector<GeneIndex>> owned_active_genes;
  std::shared_ptr<const std::vector<int64_t>> active_gene_identity_ids;
  std::shared_ptr<const std::vector<int>> active_gene_identity_taxa;
  std::shared_ptr<const HatSplitCertificateIndex> split_certificate_index;
  int active_gene_source_identity_count = 0;
  bool active_genes_match_identities = false;
  const SubproblemWeights* active_weights = nullptr;
  ConsensusTree consensus;
  const ConsensusTree* consensus_ref = nullptr;
  Options options;
  int workers = 1;
};

/// Fixed-size diagnostic summary of an ordered identity sequence.
struct BenchmarkIdentitySummary {
  int count = 0;
  uint64_t key_bytes = 0;
  uint64_t ordered_hash64 = 0;
};

/// Runtime counters emitted in benchmark JSON and debug timing output.
struct BenchmarkStats {
  std::atomic<long long> subproblems{0};
  std::atomic<long long> completed_subproblems{0};
  std::atomic<long long> max_depth{0};
  std::atomic<long long> max_active_identities{0};
  std::atomic<long long> max_active_taxa{0};
  std::atomic<long long> worst_balance_ratio_ppm{0};
  std::atomic<long long> candidate_sets{0};
  std::atomic<long long> candidates{0};
  std::atomic<long long> max_candidates{0};
  std::atomic<long long> active_consensus_candidate_sets{0};
  std::atomic<long long> active_consensus_candidates{0};
  std::atomic<long long> active_gene_split_candidates_seen{0};
  std::atomic<long long> active_gene_split_candidates_added{0};
  std::atomic<long long> active_gene_split_candidate_duplicates{0};
  std::atomic<long long> active_total_candidates_after_expansion{0};
  std::atomic<long long> attachment_calls{0};
  std::atomic<long long> attachment_internal_nodes{0};
  std::atomic<long long> attachment_four_color_passes{0};
  std::atomic<long long> attachment_edge_ties{0};
  std::atomic<long long> attachment_normalized_one_x_calls{0};
  std::atomic<long long> hat_subproblems{0};
  std::atomic<long long> hat_max_active_taxa{0};
  std::atomic<long long> hat_consensus_index_builds{0};
  std::atomic<long long> hat_consensus_index_nodes{0};
  std::atomic<long long> hat_active_local_taxon_entries{0};
  std::atomic<long long> hat_active_local_index_entries{0};
  std::atomic<long long> hat_virtual_query_selected_leaves{0};
  std::atomic<long long> hat_virtual_query_lca_queries{0};
  std::atomic<long long> hat_virtual_query_retained_nodes{0};
  std::atomic<long long> hat_localized_scoring_gene_nodes{0};
  std::atomic<long long> hat_localized_scoring_calls{0};
  // solve_hat_recursive owns these exact outer-recursion work counters.
  std::atomic<long long> hat_outer_recursive_calls{0};
  std::atomic<long long> hat_outer_nonbase_calls{0};
  std::atomic<long long> hat_outer_max_depth{0};
  std::atomic<long long> hat_outer_active_identity_mass{0};
  std::atomic<long long> hat_outer_active_support_word_mass{0};
  std::atomic<long long> hat_outer_active_taxon_mass{0};
  std::atomic<long long> hat_outer_split_identity_tests{0};
  std::atomic<long long> hat_outer_partition_identity_moves{0};
  std::atomic<long long> hat_outer_partition_support_word_moves{0};
  std::atomic<long long> hat_outer_partition_taxon_moves{0};
  std::atomic<long long> hat_outer_balanced_splits{0};
  std::atomic<long long> hat_outer_maximally_skewed_splits{0};
  std::atomic<long long> hat_outer_catalog_identity_entries{0};
  std::atomic<long long> hat_outer_catalog_taxon_entries{0};
  std::atomic<long long> hat_outer_catalog_lookup_entries{0};
  std::atomic<long long> hat_outer_tree_node_visits{0};
  std::atomic<long long> hat_outer_tree_identity_lookups{0};
  std::atomic<long long> hat_outer_tree_taxon_entries{0};
  std::atomic<long long> hat_outer_taxon_radix_entries{0};
  std::atomic<long long> hat_outer_taxon_radix_bucket_ops{0};
  std::atomic<long long> hat_outer_taxon_unique_entries{0};
  std::atomic<long long> hat_outer_root_bank_nodes{0};
  std::atomic<long long> hat_outer_root_bank_other_entries{0};

  // Exact logical payload gauges. Categories are disjoint except the
  // continuation gauges, which are role-specific subsets and are therefore
  // intentionally excluded from hat_outer_live_payload_units.
  std::atomic<long long> hat_outer_live_payload_units{0};
  std::atomic<long long> hat_outer_peak_payload_units{0};
  std::atomic<long long> hat_outer_live_identity_entries{0};
  std::atomic<long long> hat_outer_peak_identity_entries{0};
  std::atomic<long long> hat_outer_live_support_words{0};
  std::atomic<long long> hat_outer_peak_support_words{0};
  std::atomic<long long> hat_outer_live_identity_taxon_entries{0};
  std::atomic<long long> hat_outer_peak_identity_taxon_entries{0};
  std::atomic<long long> hat_outer_live_identity_key_bytes{0};
  std::atomic<long long> hat_outer_peak_identity_key_bytes{0};
  std::atomic<long long> hat_outer_live_catalog_entries{0};
  std::atomic<long long> hat_outer_peak_catalog_entries{0};
  std::atomic<long long> hat_outer_live_bank_nodes{0};
  std::atomic<long long> hat_outer_peak_bank_nodes{0};
  std::atomic<long long> hat_outer_live_bank_other_entries{0};
  std::atomic<long long> hat_outer_peak_bank_other_entries{0};
  std::atomic<long long> hat_outer_live_bank_identity_entries{0};
  std::atomic<long long> hat_outer_peak_bank_identity_entries{0};
  std::atomic<long long> hat_outer_live_continuation_identity_entries{0};
  std::atomic<long long> hat_outer_peak_continuation_identity_entries{0};
  std::atomic<long long> hat_outer_live_continuation_support_words{0};
  std::atomic<long long> hat_outer_peak_continuation_support_words{0};
  std::atomic<long long> hat_outer_live_continuation_taxon_entries{0};
  std::atomic<long long> hat_outer_peak_continuation_taxon_entries{0};
  std::atomic<long long> hat_outer_live_continuation_key_bytes{0};
  std::atomic<long long> hat_outer_peak_continuation_key_bytes{0};
  std::atomic<long long> hat_outer_live_completed_tree_nodes{0};
  std::atomic<long long> hat_outer_peak_completed_tree_nodes{0};
  std::atomic<long long> hat_outer_live_completed_tree_adjacencies{0};
  std::atomic<long long> hat_outer_peak_completed_tree_adjacencies{0};
  std::atomic<long long> hat_outer_live_taxon_scratch_entries{0};
  std::atomic<long long> hat_outer_peak_taxon_scratch_entries{0};
  std::atomic<long long> hat_outer_live_accounting_violations{0};
  std::atomic<long long> hat_consensus_edges_scored{0};
  std::atomic<long long> hat_hdt_banks_built{0};
  std::atomic<long long> hat_hdt_bank_build_ns{0};
  std::atomic<long long> hat_hdt_recolor_ops{0};
  std::atomic<long long> hat_hdt_nodes_recomputed{0};
  std::atomic<long long> hat_hdt_queries{0};
  std::atomic<long long> hat_reduced_hdt_nodes{0};
  std::atomic<long long> hat_contract_calls{0};
  std::atomic<long long> hat_extract_calls{0};
  std::atomic<long long> hat_contract_input_nodes{0};
  std::atomic<long long> hat_contract_output_nodes{0};
  std::atomic<long long> hat_contract_live_leaves{0};
  std::atomic<long long> hat_contract_max_output_per_live_ppm{0};
  std::atomic<long long> hat_extract_input_nodes{0};
  std::atomic<long long> hat_extract_output_nodes{0};
  std::atomic<long long> hat_extract_live_leaves{0};
  std::atomic<long long> hat_extract_max_output_per_live_ppm{0};
  std::atomic<long long> hat_extract_source_hdt_nodes{0};
  std::atomic<long long> hat_extract_marked_hdt_nodes{0};
  std::atomic<long long> hat_extract_hanging_components{0};
  std::atomic<long long> hat_extract_intermediate_units{0};
  std::atomic<long long> hat_extract_live_leaves_exact{0};
  std::atomic<long long> hat_contract_light_calls{0};
  std::atomic<long long> hat_contract_heavy_calls{0};
  std::atomic<long long> hat_contract_input_units{0};
  std::atomic<long long> hat_contract_output_units{0};
  std::atomic<long long> hat_contract_live_leaves_exact{0};
  std::atomic<long long> hat_contract_corridors_emitted{0};
  std::atomic<long long> hat_contract_background_units_absorbed{0};
  std::atomic<long long> hat_contract_polytomy_arms_folded{0};
  std::atomic<long long> hat_contract_cc_ops{0};
  std::atomic<long long> hat_contract_gg_ops{0};
  std::atomic<long long> hat_contract_ig_ops{0};
  std::atomic<long long> hat_contract_c_to_g_ops{0};
  std::atomic<long long> hat_contract_bound_violations{0};
  std::atomic<long long> hat_heavy_input_units{0};
  std::atomic<long long> hat_heavy_output_units{0};
  std::atomic<long long> hat_heavy_calls{0};
  std::atomic<long long> hat_heavy_shrink_violations{0};
  std::atomic<long long> hat_hdt_initial_atoms{0};
  std::atomic<long long> hat_hdt_rebuild_nodes{0};
  std::atomic<long long> hat_hdt_local_balance_violations{0};
  std::atomic<long long> hat_hdt_max_update_height{0};
  std::atomic<long long> hat_recolor_changed_leaves{0};
  std::atomic<long long> hat_recolor_marked_nodes{0};
  std::atomic<long long> hat_dp_transition_evals{0};
  std::atomic<long long> hat_dp_scalar_fmas{0};
  std::atomic<long long> hat_extract_scratch_peak_bytes{0};
  std::atomic<long long> hat_contract_scratch_peak_bytes{0};
  std::atomic<long long> hat_state_arena_peak_bytes{0};
  std::atomic<long long> hat_worker_total_peak_bytes{0};
  std::atomic<long long> hat_score_partial_peak_bytes{0};
  std::atomic<long long> hat_contract_mark_ns{0};
  std::atomic<long long> hat_contract_background_fold_ns{0};
  std::atomic<long long> hat_contract_corridor_fold_ns{0};
  std::atomic<long long> hat_contract_emit_ns{0};
  std::atomic<long long> hat_hdt_rebuild_ns{0};
  std::atomic<long long> complexity_contract_satisfied{1};
  std::atomic<long long> hat_frozen_segments{0};
  std::atomic<long long> hat_frozen_segment_source_degree_sum{0};
  std::atomic<long long> hat_recompute_child_folds{0};
  std::atomic<long long> hat_max_recompute_degree{0};
  std::atomic<long long> hat_split_score_ns{0};
  std::atomic<long long> hat_certificate_index_ns{0};
  std::atomic<long long> hat_certificate_ns{0};
  std::atomic<long long> hat_certificate_tests{0};
  std::atomic<long long> hat_certificate_accepts{0};
  std::atomic<long long> hat_certificate_rejects{0};
  std::atomic<long long> hat_certificate_negative_accepts{0};
  std::atomic<long long> hat_certificate_max_negative_gap_twice_k{0};
  std::atomic<long long> hat_certificate_singletons{0};
  std::atomic<long long> hat_certificate_distance_queries{0};
  std::atomic<long long> hat_one_x_calls{0};
  std::atomic<long long> hat_one_x_medians{0};
  std::atomic<long long> hat_one_x_ns{0};
  std::atomic<long long> hat_root_source_scoring_calls{0};
  std::atomic<long long> hat_child_gene_bank_reuses{0};
  std::atomic<long long> hat_child_bank_state_attempts{0};
  std::atomic<long long> hat_child_bank_state_accepts{0};
  std::atomic<long long> hat_child_bank_state_heavy_reuses{0};
  std::atomic<long long> hat_child_bank_state_pair_violations{0};
  std::atomic<long long> hat_child_bank_state_input_nodes{0};
  std::atomic<long long> hat_child_bank_state_output_nodes{0};
  std::atomic<long long> hat_child_bank_state_max_accepted_ratio_ppm{0};
  std::atomic<long long> hat_max_scoring_source_identities{0};
  std::atomic<long long> hat_max_scoring_source_inflation_ppm{0};
  std::atomic<long long> hat_scoring_source_inflation_violations{0};
  std::atomic<long long> hat_structural_serial_child_subproblems{0};
  std::atomic<long long> hat_structural_parallel_child_subproblems{0};
  std::atomic<long long> hat_structural_parallel_rejections{0};
  std::atomic<long long> hat_structural_parallel_physical_rejections{0};
  std::atomic<long long> hat_active_structural_parallel_slots{0};
  std::atomic<long long> hat_max_structural_parallel_slots{0};
  std::atomic<long long> hat_nni_passes_executed{0};
  std::atomic<long long> hat_nni_edges_scored{0};
  std::atomic<long long> hat_nni_full_score_passes{0};
  std::atomic<long long> hat_nni_local_score_passes{0};
  std::atomic<long long> hat_nni_local_edges_scored{0};
  std::atomic<long long> hat_nni_moves_accepted{0};
  std::atomic<long long> hat_nni_moves_rejected_by_threshold{0};
  std::atomic<long long> hat_nni_peak_proposals{0};
  std::atomic<long long> hat_nni_max_qsplit_leaves{0};
  std::atomic<long long> hat_nni_qsplit_owner_bytes_avoided{0};
  std::atomic<double> hat_nni_gain_threshold{0.0};
  std::atomic<double> hat_nni_last_accepted_gain{0.0};
  std::atomic<double> hat_nni_max_accepted_gain{0.0};
  std::atomic<double> hat_nni_verified_gain{0.0};
  std::atomic<long long> hat_nni_prefix_gate_stops{0};
  std::atomic<long long> hat_nni_prefix_gate_quartet_stops{0};
  std::atomic<long long> hat_nni_prefix_gate_split_stops{0};
  std::atomic<long long> hat_nni_split_support_delta{0};
  std::atomic<long long> hat_nni_split_support_unique_splits{0};
  std::atomic<long long> hat_nni_split_support_ns{0};
  std::atomic<long long> hat_nni_ns{0};
  std::atomic<long long> hat_peak_bank_bytes{0};
  std::atomic<long long> active_child_compressed_bank_identity_slots{0};
  std::atomic<long long> max_active_child_compressed_bank_identity_slots{0};
  std::atomic<long long> consensus_distinct_splits{0};
  std::atomic<long long> consensus_accepted_splits{0};
  std::atomic<long long> consensus_hash_low{0};
  std::atomic<long long> mre_candidates_checked{0};
  std::atomic<long long> mre_candidates_rejected{0};
  std::atomic<long long> mre_case_a_accepts{0};
  std::atomic<long long> mre_case_b_accepts{0};
  std::atomic<long long> mre_moved_children{0};
  std::atomic<long long> mre_max_child_degree{0};
  std::atomic<long long> mre_radix_passes{0};
  std::atomic<long long> mre_radix_items{0};
  std::atomic<long long> mre_child_index_operations{0};
  std::atomic<long long> mre_child_index_rotations{0};
  std::atomic<long long> mre_child_index_max_height{0};
  std::atomic<long long> mre_child_index_height_violations{0};
  std::atomic<long long> mre_child_index_enumeration_visits{0};
  std::atomic<long long> mre_topology_slots_allocated{0};
  std::atomic<long long> mre_topology_peak_live_slots{0};
  std::atomic<long long> compression_materializations{0};
  std::atomic<long long> active_compressed_banks{0};
  std::atomic<long long> max_active_compressed_banks{0};
  std::atomic<long long> active_child_compressed_banks{0};
  std::atomic<long long> max_active_child_compressed_banks{0};
  std::atomic<long long> max_child_compressed_bank_bytes{0};
  std::atomic<long long> choose_initial_ns{0};
  std::atomic<long long> attachment_range_add_ns{0};
  std::atomic<long long> attachment_score_ns{0};
  std::atomic<long long> attachment_total_ns{0};
  std::atomic<long long> parse_input_ns{0};
  std::atomic<long long> taxon_registry_ns{0};
  std::atomic<long long> taxon_registry_physical_labels{0};
  std::atomic<long long> taxon_registry_unique_labels{0};
  std::atomic<long long> taxon_registry_sort_items{0};
  std::atomic<long long> gene_index_ns{0};
  std::atomic<long long> consensus_count_ns{0};
  std::atomic<long long> consensus_assembly_ns{0};
  std::atomic<long long> startup_ns{0};
  std::atomic<long long> compression_materialization_ns{0};
  std::atomic<long long> cpu_work_units{0};
};

/// Directed-incidence endpoint moments for HAT no-dummy quartet scoring.
struct EndpointValues {
  double pi_l = 0.0;
  double pi_r = 0.0;
  double pi_lr = 0.0;
};

/// Selected edge in a solver tree for one-X attachment. A negative `b`
/// denotes the singleton-tree attachment point at node `a`.
struct AttachmentEdge {
  int a = -1;
  int b = -1;
  double rank_score = 0.0;
  double second_rank_score = 0.0;
  double rank_margin = 0.0;
  double raw_one_x_score = 0.0;
  double normalized_one_x_score = 0.0;
  int ties = 0;
  int rank = 0;
};

/// Builds a taxon registry from all input gene tree labels.
TaxonRegistry build_registry(const std::vector<GeneTree>& genes);
/// Builds traversal arrays for an indexed tree.
RootedTraversal build_rooted_traversal(const IndexedTree& tree);

/// Creates one real singleton identity per registered taxon.
std::vector<Identity> make_initial_identities(const TaxonRegistry& taxa);
/// Returns the taxon universe size implied by an identity.
int identity_universe_size(const Identity& identity);
/// Returns the maximum universe size across identities.
int identity_universe_size(const std::vector<Identity>& identities);
/// Tests whether an identity supports a taxon.
bool identity_contains_taxon(const Identity& identity, int taxon);
/// Builds taxon ownership weights for the current identity partition.
SubproblemWeights build_subproblem_weights(const std::vector<Identity>& identities,
                                           int taxon_count);

/// Returns active compressed genes when available, otherwise base genes.
const std::vector<GeneIndex>& scoring_genes(const SolverContext& ctx);
/// Returns the immutable base gene bank for a solver context.
const std::vector<GeneIndex>& base_genes(const SolverContext& ctx);
/// Returns the taxon registry for a solver context.
const TaxonRegistry& context_taxa(const SolverContext& ctx);
/// Returns the consensus tree for a solver context.
const ConsensusTree& context_consensus(const SolverContext& ctx);
/// Returns process-global benchmark counters.
BenchmarkStats& benchmark_stats();
/// Returns whether process-global telemetry collection is currently enabled.
bool benchmark_telemetry_enabled() noexcept;
/// Enables or disables telemetry helpers at a quiescent inference boundary.
void set_benchmark_telemetry_enabled(bool enabled) noexcept;
/// Resets process-global benchmark counters.
void reset_benchmark_stats();
/// Atomically adds to a benchmark counter.
void benchmark_counter_add(std::atomic<long long>& counter, long long value = 1);
/// Atomically records a benchmark maximum.
void benchmark_counter_max(std::atomic<long long>& counter, long long value);
/// Acquires one disjoint logical HAT payload category and updates total peak.
void benchmark_live_payload_acquire(std::atomic<long long>& current,
                                    std::atomic<long long>& peak,
                                    long long value);
/// Releases one disjoint logical HAT payload category.
void benchmark_live_payload_release(std::atomic<long long>& current,
                                    long long value);
/// Acquires a role-specific live gauge that is not part of payload total.
void benchmark_live_gauge_acquire(std::atomic<long long>& current,
                                  std::atomic<long long>& peak,
                                  long long value);
/// Releases a role-specific live gauge that is not part of payload total.
void benchmark_live_gauge_release(std::atomic<long long>& current,
                                  long long value);
/// Adds elapsed nanoseconds to a benchmark timing counter.
void benchmark_timer_add(std::atomic<long long>& counter,
                         std::chrono::nanoseconds elapsed);
/// Starts a structured benchmark record with a fixed-size identity summary.
void benchmark_subproblem_start(uint64_t path,
                                int depth,
                                int active_identities,
                                int active_taxa,
                                const std::vector<Identity>& identities);
/// Adds an elapsed phase total to a structured subproblem benchmark record.
void benchmark_subproblem_add_phase(uint64_t path,
                                    std::string phase,
                                    std::chrono::nanoseconds elapsed);
/// Stores fixed-size summaries of the selected initial bipartition.
void benchmark_subproblem_set_initial(uint64_t path,
                                      std::string source,
                                      int candidate_count,
                                      double score,
                                      int representative_order,
                                      const std::vector<Side>& assignment,
                                      const std::vector<Identity>& identities);
/// Stores fixed-size summaries of the final refined bipartition.
void benchmark_subproblem_set_final(uint64_t path,
                                    const std::vector<Identity>& left,
                                    const std::vector<Identity>& right);
/// Writes all structured subproblem benchmark records as a JSON array.
void write_benchmark_subproblem_records_json(std::ostream& out);

/// Returns a stable display name for a progress mode.
const char* progress_mode_name(ProgressMode mode);
/// Initializes progress reporting for one inference run.
void progress_start(const Options& options, int total_splits);
/// Emits a startup/status progress line before recursive solving begins.
void progress_startup_status(const Options& options,
                             std::string stage,
                             std::string detail = {});
/// Starts the MRE/WuHash consensus progress bar.
void progress_consensus_start(const Options& options,
                              int total_genes,
                              int taxon_count);
/// Records one gene processed by consensus split counting.
void progress_consensus_count_gene();
/// Updates consensus assembly progress.
void progress_consensus_assembly(long long checked_entries,
                                 long long total_entries,
                                 int accepted,
                                 int max_accepted);
/// Finishes the consensus progress bar.
void progress_consensus_finish(int accepted, int max_accepted);
/// Records one non-terminal recursive subproblem for progress and diagnostics.
void progress_note_subproblem(uint64_t path,
                              int depth,
                              int active_identities,
                              int active_taxa);
/// Updates the live phase label for a recursive subproblem.
void progress_note_phase(uint64_t path, std::string phase);
/// Updates the global phase label shown when no recursive subproblem is active.
void progress_note_stage(std::string stage);
/// Records the start of one HAT one-X glue scoring task.
void progress_note_one_x_start(int parent_active_size,
                               int left_taxa,
                               int right_taxa);
/// Records the end of one HAT one-X glue scoring task.
void progress_note_one_x_done(int parent_active_size);
/// Updates HAT NNI repair progress.
void progress_note_hat_nni(std::string mode,
                           int pass,
                           int pass_cap,
                           bool unbounded,
                           long long edges_scored,
                           long long moves,
                           double gain,
                           double min_gain = 0.0);
/// Records one completed recursive split and refreshes progress output.
void progress_note_split(uint64_t path,
                         int depth,
                         int active_identities,
                         int active_taxa,
                         int left_size,
                         int right_size);
/// Finishes progress output for one inference run.
void progress_finish();

/// Returns a stable display name for an MRE mode.
/// Tests whether debug trace output is enabled at a recursion depth.
bool trace_enabled_at(const Options& options, int depth);
/// Escapes a string for trace JSON output.
std::string json_escape(const std::string& value);
/// Appends one JSON object line to the configured trace file.
void trace_json_line(const Options& options, int depth, const std::string& body);
/// Returns the stable trace label for an identity.
std::string identity_trace_key(const Identity& identity);
/// Serializes identities as a JSON string array.
std::string identity_array_json(const std::vector<Identity>& identities);
/// Serializes one side of a split as identity trace keys.
std::string split_array_json(const Bitset& split,
                             const std::vector<Identity>& identities,
                             bool left_side);
/// Serializes one side of an assignment as identity trace keys.
std::string assignment_array_json(const std::vector<Side>& assignment,
                                  const std::vector<Identity>& identities,
                                  Side side);
/// Serializes all taxa represented by identities.
std::string identity_taxa_array_json(const std::vector<Identity>& identities);
/// Serializes taxa represented by one side of an assignment.
std::string assignment_taxa_array_json(const std::vector<Side>& assignment,
                                       const std::vector<Identity>& identities,
                                       Side side);

/// Builds traversal arrays for a rooted gene index.
GeneTraversal traversal(const GeneIndex& gene, int root);
/// Converts a public gene tree into internal taxon-indexed form.
GeneIndex index_gene(const GeneTree& gene, const TaxonRegistry& taxa);
/// Projects physical occurrences into canonical distinct-taxon split leaves.
GeneSplitLeafProjection make_gene_split_leaf_projection(
    const GeneIndex& gene,
    int taxon_count,
    const std::vector<std::array<uint64_t, 4>>& taxon_fingerprints);
/// Returns the current identity owner of a gene leaf node.
int leaf_owner(const GeneIndex& gene, const SubproblemWeights& weights, int node);
/// Returns the current identity mass represented by a gene leaf node.
double leaf_mass(const GeneIndex& gene, const SubproblemWeights& weights, int node);
/// Returns how many real identities a gene leaf node contributes.
double leaf_real_count(const GeneIndex& gene,
                       const std::vector<Identity>& identities,
                       const SubproblemWeights& weights,
                       int node);
/// Prunes inactive leaves and suppresses uninformative internal nodes.
GeneIndex compress_gene_for_subproblem(const GeneIndex& original,
                                       const std::vector<Identity>& identities,
                                       const SubproblemWeights& weights);
/// Prunes a previously reduced gene under an identity remapping.
GeneIndex compress_gene_for_identity_remap(
    const GeneIndex& original,
    const std::vector<Identity>& identities,
    const std::vector<int>& old_identity_to_new,
    const std::vector<double>& old_identity_scale);

/// Builds the MRE consensus tree over indexed genes.
ConsensusTree build_consensus(const std::vector<GeneIndex>& genes,
                              int taxon_count,
                              const Options& options);

struct MreComplexityProbeResult {
  bool radix_matches_oracle = false;
  bool dense_split_ids = false;
  bool avl_invariants_hold = false;
  bool range_aggregates_match = false;
  int radix_passes = 0;
  long long radix_items = 0;
  int occurrence_radix_passes = 0;
  int priority_radix_passes = 0;
  long long occurrence_radix_items = 0;
  long long priority_radix_items = 0;
  long long emitted_occurrences = 0;
  long long deduplicated_gene_splits = 0;
  long long distinct_splits = 0;
  int max_avl_height = 0;
  int avl_height_limit = 0;
  long long empty_complement_visits = 0;
  long long empty_complement_visit_limit = 0;
  bool flat_topology_matches = false;
  int topology_peak_slots = 0;
  int topology_slot_limit = 0;
};

MreComplexityProbeResult mre_complexity_probe_for_tests(int item_count);
/// Converts a parsed Newick tree into the rooted consensus representation.
ConsensusTree consensus_from_tree(const Tree& tree, const TaxonRegistry& taxa);
/// Writes the consensus tree dump when configured.
void write_consensus_dump(const Options& options, const SolverContext& ctx);
/// Builds candidate bipartitions from the consensus tree.
std::vector<Candidate> build_candidates(const ConsensusTree& consensus,
                                        const std::vector<Identity>& identities,
                                        int taxon_count);
/// Converts a split bitset into a left/right assignment vector.
std::vector<Side> assignment_from_split(const Bitset& split, int m);
/// Tests whether a recursive split makes legal progress.
bool legal_final_split_progress(const std::vector<Identity>& left,
                                const std::vector<Identity>& right);
/// Solves one HAT-v1 active-only recursive run.
SolverTree solve_hat_v1(SolverContext& ctx,
                        std::vector<Identity> identities,
                        uint64_t path,
                        int depth,
                        tf::Subflow* subflow);
/// Split shape override used only by the instrumentation regression driver.
enum class HatRecursionTestShape : uint8_t { Balanced, MaximallySkewed };
/// Runs the production HAT recursion while forcing a deterministic test shape.
SolverTree solve_hat_v1_for_test(SolverContext& ctx,
                                 std::vector<Identity> identities,
                                 HatRecursionTestShape shape);
/// Converts MRE consensus to a solver tree and optionally repairs it with HAT NNI.
SolverTree solve_hat_consensus(SolverContext& ctx,
                               bool run_nni);
/// Configures the shared Taskflow executor and returns worker count.
int configure_threads(const Options& options);
/// Returns the lazily initialized shared Taskflow executor.
tf::Executor& get_executor();
/// Returns the number of workers in the shared executor.
int get_num_workers();
/// Runs a nested Taskflow without deadlocking worker threads.
void run_inner(tf::Taskflow& taskflow);
/// Mixes three integers into a deterministic 64-bit hash.
uint64_t stable_hash_u64(uint64_t a, uint64_t b, uint64_t c);

inline thread_local bool tl_allow_inner_parallelism = true;
inline thread_local int tl_inner_parallel_worker_budget = 0;

class ScopedInnerParallelBudgetSuspend {
 public:
  ScopedInnerParallelBudgetSuspend()
      : saved_budget_(tl_inner_parallel_worker_budget) {
    tl_inner_parallel_worker_budget = 1;
  }

  ScopedInnerParallelBudgetSuspend(
      const ScopedInnerParallelBudgetSuspend&) = delete;
  ScopedInnerParallelBudgetSuspend& operator=(
      const ScopedInnerParallelBudgetSuspend&) = delete;

  ~ScopedInnerParallelBudgetSuspend() {
    tl_inner_parallel_worker_budget = saved_budget_;
  }

 private:
  int saved_budget_;
};

/// Returns the active worker budget for inner parallel loops.
inline int current_inner_parallel_worker_limit() {
  if (!tl_allow_inner_parallelism) {
    return 1;
  }
  const int workers = std::max(1, get_num_workers());
  if (tl_inner_parallel_worker_budget <= 0) {
    return workers;
  }
  return std::max(1, std::min(workers, tl_inner_parallel_worker_budget));
}

struct ParallelLoopFailure {
  int index = std::numeric_limits<int>::max();
  std::exception_ptr exception;

  void capture(int logical_index) noexcept {
    // Each slot belongs to one ascending participant, so its first failure is
    // also that participant's lowest-index failure.
    if (!exception) {
      index = logical_index;
      exception = std::current_exception();
    }
  }
};

class ParallelLoopJoinState {
 public:
  explicit ParallelLoopJoinState(int helpers) : remaining_(helpers) {}

  void complete(int count = 1) noexcept {
    if (remaining_.fetch_sub(count, std::memory_order_acq_rel) == count) {
      remaining_.notify_one();
    }
  }

  int remaining() const noexcept {
    return remaining_.load(std::memory_order_acquire);
  }

  void wait(int observed) const noexcept {
    remaining_.wait(observed, std::memory_order_acquire);
  }

 private:
  std::atomic<int> remaining_;
};

class ParallelLoopCompletion {
 public:
  explicit ParallelLoopCompletion(ParallelLoopJoinState& join) noexcept
      : join_(join) {}
  ParallelLoopCompletion(const ParallelLoopCompletion&) = delete;
  ParallelLoopCompletion& operator=(const ParallelLoopCompletion&) = delete;
  ~ParallelLoopCompletion() { join_.complete(); }

 private:
  ParallelLoopJoinState& join_;
};

template <typename F>
void run_parallel_loop_range(int begin,
                             int end,
                             F& body,
                             ParallelLoopFailure& failure) noexcept {
  // Parallel loops finish all iterations after a callback failure. This keeps
  // completion bounded and makes the lowest failing logical index observable.
  for (int i = begin; i < end; ++i) {
    try {
      body(i);
    } catch (...) {
      failure.capture(i);
    }
  }
}

inline void wait_for_parallel_loop(
    tf::Executor& exec,
    const std::shared_ptr<ParallelLoopJoinState>& join) {
  if (exec.this_worker_id() >= 0) {
    ScopedInnerParallelBudgetSuspend suspend_budget;
    exec.corun_until([&] { return join->remaining() <= 0; });
  } else {
    int observed = join->remaining();
    while (observed > 0) {
      join->wait(observed);
      observed = join->remaining();
    }
  }
}

inline std::exception_ptr select_parallel_loop_failure(
    const std::vector<ParallelLoopFailure>& failures) noexcept {
  const ParallelLoopFailure* selected = nullptr;
  for (const ParallelLoopFailure& failure : failures) {
    if (failure.exception &&
        (selected == nullptr || failure.index < selected->index)) {
      selected = &failure;
    }
  }
  return selected == nullptr ? std::exception_ptr{} : selected->exception;
}

/// Runs a small work-sharing loop on the shared executor.
template <typename F>
void parallel_for_inner(int begin, int end, F&& body, int min_tasks = 4) {
  const int n = end - begin;
  if (n <= 0) return;
  if (!tl_allow_inner_parallelism) {
    for (int i = begin; i < end; ++i) body(i);
    return;
  }
  auto& exec = get_executor();
  const int nw = std::min(current_inner_parallel_worker_limit(), n);
  if (nw <= 1 || n < min_tasks) {
    for (int i = begin; i < end; ++i) body(i);
    return;
  }
  std::vector<ParallelLoopFailure> failures(static_cast<std::size_t>(nw));
  auto join = std::make_shared<ParallelLoopJoinState>(nw - 1);
  std::exception_ptr submission_failure;
  int first_unsubmitted = nw;
  for (int t = 1; t < nw; ++t) {
    int lo = begin + static_cast<int>(static_cast<int64_t>(n) * t / nw);
    int hi = begin + static_cast<int>(static_cast<int64_t>(n) * (t + 1) / nw);
    try {
      exec.silent_async(
          [lo, hi, &body, &failure = failures[t], join]() noexcept {
            ParallelLoopCompletion completion(*join);
            run_parallel_loop_range(lo, hi, body, failure);
          });
    } catch (...) {
      submission_failure = std::current_exception();
      first_unsubmitted = t;
      join->complete(nw - t);
      break;
    }
  }
  int hi0 = begin + static_cast<int>(static_cast<int64_t>(n) / nw);
  run_parallel_loop_range(begin, hi0, body, failures[0]);
  for (int t = first_unsubmitted; t < nw; ++t) {
    int lo = begin + static_cast<int>(static_cast<int64_t>(n) * t / nw);
    int hi = begin + static_cast<int>(static_cast<int64_t>(n) * (t + 1) / nw);
    run_parallel_loop_range(lo, hi, body, failures[t]);
  }
  wait_for_parallel_loop(exec, join);
  if (std::exception_ptr failure = select_parallel_loop_failure(failures)) {
    std::rethrow_exception(failure);
  }
  if (submission_failure) {
    std::rethrow_exception(submission_failure);
  }
}

/// Like parallel_for_inner but assigns indices [begin, end) to workers via a
/// shared atomic cursor (work-stealing) instead of a static contiguous split.
/// Use when per-index cost is uneven: callers should pass a modest number of
/// indices (e.g. a few times the worker count) so each is a coarse chunk that
/// amortises the atomic, and must keep any cross-index reduction ordered by
/// index (not completion order) to stay deterministic.
template <typename F>
void parallel_for_dynamic(int begin, int end, F&& body, int min_tasks = 2) {
  const int n = end - begin;
  if (n <= 0) return;
  if (!tl_allow_inner_parallelism) {
    for (int i = begin; i < end; ++i) body(i);
    return;
  }
  auto& exec = get_executor();
  const int nw = std::min(current_inner_parallel_worker_limit(), n);
  if (nw <= 1 || n < min_tasks) {
    for (int i = begin; i < end; ++i) body(i);
    return;
  }
  std::atomic<int> cursor{begin};
  std::vector<ParallelLoopFailure> failures(static_cast<std::size_t>(nw));
  auto join = std::make_shared<ParallelLoopJoinState>(nw - 1);
  std::exception_ptr submission_failure;
  auto drain = [&](ParallelLoopFailure& failure) noexcept {
    int i;
    while ((i = cursor.fetch_add(1, std::memory_order_relaxed)) < end) {
      try {
        body(i);
      } catch (...) {
        failure.capture(i);
      }
    }
  };
  for (int t = 1; t < nw; ++t) {
    try {
      exec.silent_async([&drain, &failure = failures[t], join]() noexcept {
        ParallelLoopCompletion completion(*join);
        drain(failure);
      });
    } catch (...) {
      submission_failure = std::current_exception();
      join->complete(nw - t);
      break;
    }
  }
  drain(failures[0]);
  wait_for_parallel_loop(exec, join);
  if (std::exception_ptr failure = select_parallel_loop_failure(failures)) {
    std::rethrow_exception(failure);
  }
  if (submission_failure) {
    std::rethrow_exception(submission_failure);
  }
}

extern std::atomic<bool> g_debug_timing_enabled;

/// Returns whether debug timers and counters should record.
inline bool debug_timing_enabled() {
  return g_debug_timing_enabled.load(std::memory_order_relaxed);
}

/// Writes a debug log line when timing is enabled.
void debug_log_line_impl(std::string_view message);
/// Adds to a named debug counter.
void debug_counter_add_impl(std::string_view name, long long delta);
/// Records a named debug maximum.
void debug_counter_max_impl(std::string_view name, long long value);
/// Adds elapsed time to a named debug timer.
void debug_timer_add_impl(std::string_view name, std::chrono::nanoseconds elapsed);
/// Emits the debug timing report.
void debug_report_impl();

/// Conditionally writes a debug log line.
inline void debug_log_line(std::string_view message) {
  if (debug_timing_enabled()) {
    debug_log_line_impl(message);
  }
}

/// Conditionally adds to a debug counter.
inline void debug_counter_add(std::string_view name, long long delta = 1) {
  if (debug_timing_enabled()) {
    debug_counter_add_impl(name, delta);
  }
}

/// Conditionally records a debug maximum.
inline void debug_counter_max(std::string_view name, long long value) {
  if (debug_timing_enabled()) {
    debug_counter_max_impl(name, value);
  }
}

/// Conditionally records elapsed debug time.
inline void debug_timer_add(std::string_view name, std::chrono::nanoseconds elapsed) {
  if (debug_timing_enabled()) {
    debug_timer_add_impl(name, elapsed);
  }
}

/// Conditionally emits the debug timing report.
inline void debug_report() {
  if (debug_timing_enabled()) {
    debug_report_impl();
  }
}

/// RAII timer that records elapsed time under a debug counter name.
class DebugScope {
 public:
  explicit DebugScope(std::string_view name)
      : name_(name),
        enabled_(debug_timing_enabled()),
        start_(enabled_ ? std::chrono::steady_clock::now()
                        : std::chrono::steady_clock::time_point{}) {}
  DebugScope(const DebugScope&) = delete;
  DebugScope& operator=(const DebugScope&) = delete;
  ~DebugScope() {
    if (enabled_) {
      debug_timer_add(name_, std::chrono::steady_clock::now() - start_);
    }
  }

 private:
  std::string_view name_;
  bool enabled_ = false;
  std::chrono::steady_clock::time_point start_;
};

}  // namespace wqfm::internal
