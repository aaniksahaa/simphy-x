#pragma once

#include "hat/hat_internal.hpp"
#include "hat/hdt/hdt_contract.hpp"

#include <array>
#include <cstdint>
#include <functional>
#include <vector>

namespace wqfm::internal {

struct HdtQueryOwnerIndex {
  const std::vector<int>* query_owner_order = nullptr;
  std::vector<int> owner_query_position;
  std::vector<int> represented_owner_order;
  std::vector<int> represented_prefix;
  std::vector<int> live_atom_prefix;
};

struct HdtTraversalStats {
  long long states_visited = 0;
  long long score_reads = 0;

  long long extract_calls = 0;
  long long extract_source_hdt_nodes = 0;
  long long extract_marked_hdt_nodes = 0;
  long long extract_hanging_components = 0;
  long long extract_intermediate_units = 0;
  long long extract_live_leaves = 0;

  long long contract_calls = 0;
  long long contract_light_calls = 0;
  long long contract_heavy_calls = 0;
  long long contract_input_units = 0;
  long long contract_output_units = 0;
  long long contract_live_leaves = 0;
  long long contract_corridors_emitted = 0;
  long long contract_background_units_absorbed = 0;
  long long contract_polytomy_arms_folded = 0;
  long long contract_cc_ops = 0;
  long long contract_gg_ops = 0;
  long long contract_ig_ops = 0;
  long long contract_c_to_g_ops = 0;
  long long contract_bound_violations = 0;

  long long heavy_input_units = 0;
  long long heavy_output_units = 0;
  long long heavy_calls = 0;
  long long heavy_shrink_violations = 0;

  long long hdt_initial_atoms = 0;
  long long hdt_rebuild_nodes = 0;
  int hdt_max_recompute_degree = 0;
  long long hdt_local_balance_violations = 0;
  int hdt_max_update_height = 0;

  long long recolor_changed_leaves = 0;
  long long recolor_marked_nodes = 0;
  long long dp_transition_evals = 0;
  long long dp_scalar_fmas = 0;

  int max_state_units = 0;
  int max_state_hdt_nodes = 0;
  long long pending_state_units_peak = 0;
  long long pending_state_hdt_nodes_peak = 0;
  std::size_t extract_scratch_peak_bytes = 0;
  std::size_t contract_scratch_peak_bytes = 0;
  std::size_t state_arena_peak_bytes = 0;
  std::size_t worker_total_peak_bytes = 0;
  long long contract_mark_ns = 0;
  long long contract_background_fold_ns = 0;
  long long contract_corridor_fold_ns = 0;
  long long contract_emit_ns = 0;
  long long hdt_rebuild_ns = 0;
  bool complexity_contract_satisfied = true;
};

void benchmark_merge_hdt_traversal_stats(const HdtTraversalStats& stats);

void hdt_account_initial_state(
    HdtInitialStateResult& initial,
    HdtTraversalStats& stats);

HdtQueryOwnerIndex build_hdt_query_owner_index(
    const HdtState& base,
    const HatQueryTree& query,
    const std::vector<int>& query_owner_order);

HdtLiveSelector hdt_query_node_selector(
    const HdtQueryOwnerIndex& index,
    const HatQueryTree& query,
    int query_node);

using HdtActiveSplitVisitor =
    std::function<void(int query_node, HdtState& state)>;

void hdt_traverse_active_splits(
    HdtState state,
    const HatQueryTree& query,
    const std::vector<int>& query_owner_order,
    const std::vector<unsigned char>& score_node,
    HdtExtractWorkspace& workspace,
    HdtTraversalStats& stats,
    const HdtActiveSplitVisitor& visitor,
    int state_reuse_bound_factor = 4);

using HdtOneXVisitor = std::function<void(
    int median_node,
    const std::array<int, 3>& branch_nodes,
    const std::array<std::uint8_t, 3>& branch_colors,
    HdtState& state)>;

void hdt_traverse_one_x_medians(
    HdtState state,
    const HatQueryTree& query,
    const std::vector<int>& query_owner_order,
    HdtExtractWorkspace& workspace,
    HdtTraversalStats& stats,
    const HdtOneXVisitor& visitor,
    int state_reuse_bound_factor = 4);

}  // namespace wqfm::internal
