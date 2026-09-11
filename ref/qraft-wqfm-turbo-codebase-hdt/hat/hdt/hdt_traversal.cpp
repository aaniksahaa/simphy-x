#include "hat/hdt/hdt_traversal.hpp"

#include <algorithm>
#include <memory>
#include <stdexcept>
#include <utility>

namespace wqfm::internal {
namespace {

struct StateCounters {
  long long changed = 0;
  long long marked = 0;
  long long transitions = 0;
  long long thread_scalar_fmas = 0;
};

StateCounters state_counters(const HdtState& state) {
  return StateCounters{state.recolor_changed_leaves,
                       state.recolor_marked_nodes,
                       state.transition_evals,
                       full_bdp_thread_scalar_fmas()};
}

void add_counter_delta(HdtTraversalStats& stats,
                       const HdtState& state,
                       const StateCounters& before) {
  stats.recolor_changed_leaves += state.recolor_changed_leaves - before.changed;
  stats.recolor_marked_nodes += state.recolor_marked_nodes - before.marked;
  stats.dp_transition_evals += state.transition_evals - before.transitions;
  stats.dp_scalar_fmas +=
      full_bdp_thread_scalar_fmas() - before.thread_scalar_fmas;
}

void note_hdt(const HdtState& state,
              HdtTraversalStats& stats,
              bool count_scalar_fmas = true) {
  const HdtBuildStats& hdt = state.hdt.stats;
  stats.hdt_initial_atoms += hdt.initial_atoms;
  stats.hdt_rebuild_nodes += hdt.reachable_nodes;
  stats.hdt_max_recompute_degree =
      std::max(stats.hdt_max_recompute_degree, hdt.max_recompute_degree);
  stats.hdt_local_balance_violations += hdt.local_balance_violations;
  stats.hdt_max_update_height =
      std::max(stats.hdt_max_update_height, hdt.max_update_height);
  stats.dp_transition_evals +=
      hdt.gg_ops + hdt.ig_ops + hdt.cc_ops + hdt.c_to_g_ops;
  if (count_scalar_fmas) {
    stats.dp_scalar_fmas += hdt.scalar_fmas +
                            state.components.construction_scalar_fmas;
  }
  stats.max_state_units = std::max(
      stats.max_state_units, state.components.component_units());
  stats.max_state_hdt_nodes =
      std::max(stats.max_state_hdt_nodes, hdt.reachable_nodes);
  if (!hdt.node_count_ok || !hdt.atom_bound_ok ||
      hdt.max_recompute_degree > 2 || hdt.local_balance_violations != 0) {
    stats.complexity_contract_satisfied = false;
    throw std::logic_error(
        "HDT state violated a structural invariant");
  }
}

void note_contract(const HdtContractStats& contract,
                   const HdtState& state,
                   HdtTraversalStats& stats) {
  ++stats.contract_calls;
  if (contract.mode == HdtContractMode::HeavyReplace) {
    ++stats.contract_heavy_calls;
  } else {
    ++stats.contract_light_calls;
  }
  stats.contract_input_units += contract.input_units;
  stats.contract_output_units += contract.output_units;
  stats.contract_live_leaves += contract.selected_live_atoms;
  stats.contract_corridors_emitted += contract.corridors_emitted;
  stats.contract_background_units_absorbed +=
      contract.background_units_absorbed;
  stats.contract_polytomy_arms_folded += contract.polytomy_arms_folded;
  stats.contract_cc_ops += contract.cc_ops;
  stats.contract_gg_ops += contract.gg_ops;
  stats.contract_ig_ops += contract.ig_ops;
  stats.contract_c_to_g_ops += contract.c_to_g_ops;
  stats.dp_transition_evals +=
      contract.cc_ops + contract.gg_ops + contract.ig_ops +
      contract.c_to_g_ops;
  stats.dp_scalar_fmas += contract.scalar_fmas;
  stats.contract_scratch_peak_bytes = std::max(
      stats.contract_scratch_peak_bytes, contract.scratch_peak_bytes);
  stats.contract_mark_ns += contract.mark_ns;
  stats.contract_background_fold_ns += contract.background_fold_ns;
  stats.contract_corridor_fold_ns += contract.corridor_fold_ns;
  stats.contract_emit_ns += contract.emit_ns;
  stats.hdt_rebuild_ns += contract.hdt_rebuild_ns;
  if (!contract.size_bound_ok) {
    ++stats.contract_bound_violations;
    stats.complexity_contract_satisfied = false;
  }
  note_hdt(state, stats, false);
}

void note_contract(const HdtContractResult& result,
                   HdtTraversalStats& stats) {
  note_contract(result.stats, result.state, stats);
}

void note_extract(const HdtExtractContractResult& result,
                  HdtTraversalStats& stats,
                  bool heavy = false) {
  ++stats.extract_calls;
  stats.extract_source_hdt_nodes += result.extract.source_hdt_nodes;
  stats.extract_marked_hdt_nodes += result.extract.marked_hdt_nodes;
  stats.extract_hanging_components += result.extract.hanging_components;
  stats.extract_intermediate_units += result.extract.intermediate_units;
  stats.extract_live_leaves += result.extract.selected_live_atoms;

  // Account the immediate Contract and final rebuild without copying the
  // returned state.
  ++stats.contract_calls;
  if (heavy) {
    ++stats.contract_heavy_calls;
  } else {
    ++stats.contract_light_calls;
  }
  stats.contract_input_units += result.contract.input_units;
  stats.contract_output_units += result.contract.output_units;
  stats.contract_live_leaves += result.contract.selected_live_atoms;
  stats.contract_corridors_emitted += result.contract.corridors_emitted;
  stats.contract_background_units_absorbed +=
      result.contract.background_units_absorbed;
  stats.contract_polytomy_arms_folded +=
      result.contract.polytomy_arms_folded;
  stats.contract_cc_ops += result.contract.cc_ops;
  stats.contract_gg_ops += result.contract.gg_ops;
  stats.contract_ig_ops += result.contract.ig_ops;
  stats.contract_c_to_g_ops += result.contract.c_to_g_ops;
  stats.dp_transition_evals +=
      result.contract.cc_ops + result.contract.gg_ops +
      result.contract.ig_ops + result.contract.c_to_g_ops;
  stats.dp_scalar_fmas += result.contract.scalar_fmas;
  stats.contract_scratch_peak_bytes = std::max(
      stats.contract_scratch_peak_bytes,
      result.contract.scratch_peak_bytes);
  stats.contract_mark_ns += result.contract.mark_ns;
  stats.contract_background_fold_ns += result.contract.background_fold_ns;
  stats.contract_corridor_fold_ns += result.contract.corridor_fold_ns;
  stats.contract_emit_ns += result.contract.emit_ns;
  stats.hdt_rebuild_ns += result.contract.hdt_rebuild_ns;
  if (!result.contract.size_bound_ok) {
    ++stats.contract_bound_violations;
    stats.complexity_contract_satisfied = false;
  }
  note_hdt(result.state, stats, false);
}

std::pair<int, int> query_owner_range(const HatQueryTree& query,
                                      int node,
                                      int owner_order_size) {
  if (node < 0 ||
      static_cast<std::size_t>(node) >= query.subtree_leaf_first.size() ||
      static_cast<std::size_t>(node) >= query.subtree_leaf_count.size()) {
    return {0, 0};
  }
  const int begin = std::clamp(
      query.subtree_leaf_first[static_cast<std::size_t>(node)], 0,
      owner_order_size);
  const int end = std::clamp(
      checked_hdt_int_add(
          begin, query.subtree_leaf_count[static_cast<std::size_t>(node)],
          "query owner interval overflow"),
      begin,
      owner_order_size);
  return {begin, end};
}

void recolor_query_node(HdtState& state,
                        const HatQueryTree& query,
                        const std::vector<int>& owner_order,
                        int node,
                        std::uint8_t color,
                        HdtTraversalStats& stats) {
  const int owner_count = checked_hdt_int_size(
      owner_order.size(), "query exceeds integer owner domain");
  const auto [begin, end] =
      query_owner_range(query, node, owner_count);
  const StateCounters before = state_counters(state);
  hdt_recolor_owner_range(state, owner_order, begin, end, color);
  add_counter_delta(stats, state, before);
}

void recolor_light_children(HdtState& state,
                            const HatQueryTree& query,
                            const std::vector<int>& owner_order,
                            int node,
                            int heavy,
                            HdtTraversalStats& stats,
                            std::vector<std::pair<int, int>>& owner_ranges) {
  owner_ranges.clear();
  const int owner_count = checked_hdt_int_size(
      owner_order.size(), "query exceeds integer owner domain");
  const int child_begin = query.child_first[static_cast<std::size_t>(node)];
  const int child_end = child_begin +
                        query.child_count[static_cast<std::size_t>(node)];
  for (int child_pos = child_begin; child_pos < child_end; ++child_pos) {
    const int child = query.child_nodes[static_cast<std::size_t>(child_pos)];
    if (child == heavy) continue;
    const auto [begin, end] =
        query_owner_range(query, child, owner_count);
    owner_ranges.emplace_back(begin, end);
  }
  const StateCounters before = state_counters(state);
  hdt_recolor_owner_ranges(state, owner_order, owner_ranges, 0);
  add_counter_delta(stats, state, before);
}

bool active_should_score(const HatQueryTree& query,
                         const std::vector<unsigned char>& score_node,
                         int node) {
  return node != query.root &&
         (score_node.empty() ||
          (node >= 0 && static_cast<std::size_t>(node) < score_node.size() &&
           score_node[static_cast<std::size_t>(node)] != 0));
}

bool is_query_leaf(const HatQueryTree& query, int node) {
  return node < 0 || static_cast<std::size_t>(node) >= query.nodes.size() ||
         query.nodes[static_cast<std::size_t>(node)].taxon_active >= 0;
}

bool can_reuse_state_for_child(const HdtState& state,
                               const HdtLiveSelector& selector,
                               int initial_units,
                               int reuse_bound_factor) {
  const int live = hdt_selected_atom_count(state, selector);
  const int bound = hdt_contract_size_bound(live);
  return static_cast<long long>(initial_units) <=
         static_cast<long long>(reuse_bound_factor) * bound;
}

void maybe_heavy_contract(HdtState& state,
                          const HdtLiveSelector& selector,
                          HdtExtractWorkspace& workspace,
                          HdtTraversalStats& stats,
                          int reuse_bound_factor) {
  const int live = hdt_selected_atom_count(state, selector);
  const int bound = hdt_contract_size_bound(live);
  const int input_units = state.components.component_units();
  if (static_cast<long long>(input_units) <=
      static_cast<long long>(reuse_bound_factor) * bound) {
    return;
  }
  HdtExtractContractResult contracted =
      hdt_extract_contract(state, selector, workspace);
  ++stats.heavy_calls;
  stats.heavy_input_units += input_units;
  stats.heavy_output_units += contracted.contract.output_units;
  if (2LL * static_cast<long long>(contracted.contract.output_units) >=
      input_units) {
    ++stats.heavy_shrink_violations;
    stats.complexity_contract_satisfied = false;
    throw std::logic_error("heavy Contract did not halve R");
  }
  note_extract(contracted, stats, true);
  state = std::move(contracted.state);
}

struct TraversalTask {
  int node = -1;
  HdtState state;
};

struct SavedLight {
  int node = -1;
  HdtState state;
};

void finalize_memory_stats(const HdtExtractWorkspace& workspace,
                           HdtTraversalStats& stats) {
  if (stats.pending_state_units_peak < 0 ||
      stats.pending_state_hdt_nodes_peak < 0) {
    throw std::logic_error("pending-state accounting underflowed");
  }
  const std::size_t state_bytes = checked_hdt_size_mul(
      static_cast<std::size_t>(stats.pending_state_units_peak),
      sizeof(HdtComponentUnit) +
          2u * kFullBdpMaxDimension * sizeof(double),
      "state-arena byte count overflow");
  const std::size_t hdt_bytes = checked_hdt_size_mul(
      static_cast<std::size_t>(stats.pending_state_hdt_nodes_peak),
      sizeof(HdtNode) + sizeof(HdtExpandedAtom) +
          2u * kFullBdpMaxDimension * sizeof(double),
      "HDT-arena byte count overflow");
  stats.extract_scratch_peak_bytes =
      std::max(stats.extract_scratch_peak_bytes, workspace.peak_bytes);
  stats.state_arena_peak_bytes = std::max(
      stats.state_arena_peak_bytes,
      checked_hdt_size_add(state_bytes, hdt_bytes,
                           "state-arena byte count overflow"));
  stats.worker_total_peak_bytes = std::max(
      stats.worker_total_peak_bytes,
      checked_hdt_size_add(
          checked_hdt_size_add(stats.state_arena_peak_bytes,
                               stats.extract_scratch_peak_bytes,
                               "worker byte count overflow"),
          stats.contract_scratch_peak_bytes,
          "worker byte count overflow"));
}

}  // namespace

void hdt_account_initial_state(
    HdtInitialStateResult& initial,
    HdtTraversalStats& stats) {
  if (!initial.contracted) return;
  note_contract(initial.contract, initial.state, stats);
  initial.state.initial_hdt_accounted = true;
}

void benchmark_merge_hdt_traversal_stats(
    const HdtTraversalStats& stats) {
  BenchmarkStats& out = benchmark_stats();
  benchmark_counter_add(out.hat_contract_calls, stats.contract_calls);
  benchmark_counter_add(out.hat_extract_calls, stats.extract_calls);
  benchmark_counter_add(out.hat_contract_input_nodes,
                        stats.contract_input_units);
  benchmark_counter_add(out.hat_contract_output_nodes,
                        stats.contract_output_units);
  benchmark_counter_add(out.hat_contract_live_leaves,
                        stats.contract_live_leaves);
  benchmark_counter_add(out.hat_extract_input_nodes,
                        stats.extract_source_hdt_nodes);
  benchmark_counter_add(out.hat_extract_output_nodes,
                        stats.extract_intermediate_units);
  benchmark_counter_add(out.hat_extract_live_leaves,
                        stats.extract_live_leaves);
  benchmark_counter_add(out.hat_reduced_hdt_nodes, stats.hdt_rebuild_nodes);
  benchmark_counter_max(out.hat_max_recompute_degree,
                        stats.hdt_max_recompute_degree);

  benchmark_counter_add(out.hat_extract_source_hdt_nodes,
                        stats.extract_source_hdt_nodes);
  benchmark_counter_add(out.hat_extract_marked_hdt_nodes,
                        stats.extract_marked_hdt_nodes);
  benchmark_counter_add(out.hat_extract_hanging_components,
                        stats.extract_hanging_components);
  benchmark_counter_add(out.hat_extract_intermediate_units,
                        stats.extract_intermediate_units);
  benchmark_counter_add(out.hat_extract_live_leaves_exact,
                        stats.extract_live_leaves);
  benchmark_counter_add(out.hat_contract_light_calls,
                        stats.contract_light_calls);
  benchmark_counter_add(out.hat_contract_heavy_calls,
                        stats.contract_heavy_calls);
  benchmark_counter_add(out.hat_contract_input_units,
                        stats.contract_input_units);
  benchmark_counter_add(out.hat_contract_output_units,
                        stats.contract_output_units);
  benchmark_counter_add(out.hat_contract_live_leaves_exact,
                        stats.contract_live_leaves);
  benchmark_counter_add(out.hat_contract_corridors_emitted,
                        stats.contract_corridors_emitted);
  benchmark_counter_add(out.hat_contract_background_units_absorbed,
                        stats.contract_background_units_absorbed);
  benchmark_counter_add(out.hat_contract_polytomy_arms_folded,
                        stats.contract_polytomy_arms_folded);
  benchmark_counter_add(out.hat_contract_cc_ops, stats.contract_cc_ops);
  benchmark_counter_add(out.hat_contract_gg_ops, stats.contract_gg_ops);
  benchmark_counter_add(out.hat_contract_ig_ops, stats.contract_ig_ops);
  benchmark_counter_add(out.hat_contract_c_to_g_ops,
                        stats.contract_c_to_g_ops);
  benchmark_counter_add(out.hat_contract_bound_violations,
                        stats.contract_bound_violations);
  benchmark_counter_add(out.hat_heavy_input_units, stats.heavy_input_units);
  benchmark_counter_add(out.hat_heavy_output_units, stats.heavy_output_units);
  benchmark_counter_add(out.hat_heavy_calls, stats.heavy_calls);
  benchmark_counter_add(out.hat_heavy_shrink_violations,
                        stats.heavy_shrink_violations);
  benchmark_counter_add(out.hat_hdt_initial_atoms, stats.hdt_initial_atoms);
  benchmark_counter_add(out.hat_hdt_rebuild_nodes, stats.hdt_rebuild_nodes);
  benchmark_counter_add(out.hat_hdt_local_balance_violations,
                        stats.hdt_local_balance_violations);
  benchmark_counter_max(out.hat_hdt_max_update_height,
                        stats.hdt_max_update_height);
  benchmark_counter_add(out.hat_recolor_changed_leaves,
                        stats.recolor_changed_leaves);
  benchmark_counter_add(out.hat_recolor_marked_nodes,
                        stats.recolor_marked_nodes);
  benchmark_counter_add(out.hat_dp_transition_evals,
                        stats.dp_transition_evals);
  benchmark_counter_add(out.hat_dp_scalar_fmas, stats.dp_scalar_fmas);
  benchmark_counter_max(
      out.hat_extract_scratch_peak_bytes,
      checked_hdt_long_long_size(
          stats.extract_scratch_peak_bytes,
          "Extract scratch telemetry exceeds its counter domain"));
  benchmark_counter_max(
      out.hat_contract_scratch_peak_bytes,
      checked_hdt_long_long_size(
          stats.contract_scratch_peak_bytes,
          "Contract scratch telemetry exceeds its counter domain"));
  benchmark_counter_max(
      out.hat_state_arena_peak_bytes,
      checked_hdt_long_long_size(
          stats.state_arena_peak_bytes,
          "state telemetry exceeds its counter domain"));
  benchmark_counter_max(
      out.hat_worker_total_peak_bytes,
      checked_hdt_long_long_size(
          stats.worker_total_peak_bytes,
          "worker telemetry exceeds its counter domain"));
  benchmark_counter_add(out.hat_contract_mark_ns, stats.contract_mark_ns);
  benchmark_counter_add(out.hat_contract_background_fold_ns,
                        stats.contract_background_fold_ns);
  benchmark_counter_add(out.hat_contract_corridor_fold_ns,
                        stats.contract_corridor_fold_ns);
  benchmark_counter_add(out.hat_contract_emit_ns, stats.contract_emit_ns);
  benchmark_counter_add(out.hat_hdt_rebuild_ns, stats.hdt_rebuild_ns);
  if (!stats.complexity_contract_satisfied) {
    out.complexity_contract_satisfied.store(0, std::memory_order_relaxed);
  }
}

HdtQueryOwnerIndex build_hdt_query_owner_index(
    const HdtState& base,
    const HatQueryTree& query,
    const std::vector<int>& query_owner_order) {
  if (base.components.owner_count < 0) {
    throw std::invalid_argument("negative query owner domain");
  }
  const int query_owner_count = checked_hdt_int_size(
      query_owner_order.size(),
      "query exceeds its integer owner domain");
  const int owner_offset_count = checked_hdt_int_size(
      base.hdt.owner_offsets.size(),
      "owner index exceeds its integer domain");
  HdtQueryOwnerIndex index;
  index.query_owner_order = &query_owner_order;
  index.owner_query_position.assign(
      static_cast<std::size_t>(base.components.owner_count), -1);
  const std::size_t prefix_size = checked_hdt_size_add(
      query_owner_order.size(), 1u, "query prefix capacity overflow");
  index.represented_prefix.assign(prefix_size, 0);
  index.live_atom_prefix.assign(prefix_size, 0);
  index.represented_owner_order.reserve(query_owner_order.size());

  for (int position = 0; position < query_owner_count; ++position) {
    const int owner = query_owner_order[static_cast<std::size_t>(position)];
    int live_atoms = 0;
    if (owner >= 0 && owner < base.components.owner_count) {
      int& prior = index.owner_query_position[static_cast<std::size_t>(owner)];
      if (prior >= 0) {
        throw std::invalid_argument(
            "query owner appears at more than one leaf position");
      }
      prior = position;
      const int owner_slot = base.hdt.owner_index_by_position
                                 ? position - base.hdt.owner_position_begin
                                 : owner;
      if (owner_slot >= 0 && owner_slot < owner_offset_count - 1) {
        live_atoms =
            base.hdt.owner_offsets[static_cast<std::size_t>(owner_slot) + 1u] -
            base.hdt.owner_offsets[static_cast<std::size_t>(owner_slot)];
      }
      if (live_atoms > 0) index.represented_owner_order.push_back(owner);
    }
    index.represented_prefix[static_cast<std::size_t>(position) + 1u] =
        checked_hdt_int_size(
            index.represented_owner_order.size(),
            "represented-owner index exceeds integer domain");
    index.live_atom_prefix[static_cast<std::size_t>(position) + 1u] =
        checked_hdt_int_add(
            index.live_atom_prefix[static_cast<std::size_t>(position)],
            live_atoms, "live-atom prefix overflow");
  }

  if (query.root >= 0) {
    const auto [begin, end] = query_owner_range(
        query, query.root, query_owner_count);
    if (begin != 0 || end != query_owner_count) {
      throw std::invalid_argument(
          "query root leaf interval does not cover owner order");
    }
  }
  return index;
}

HdtLiveSelector hdt_query_node_selector(
    const HdtQueryOwnerIndex& index,
    const HatQueryTree& query,
    int query_node) {
  if (index.query_owner_order == nullptr) {
    throw std::invalid_argument("query owner index is unbound");
  }
  const auto [begin, end] = query_owner_range(
      query, query_node, static_cast<int>(index.query_owner_order->size()));
  return hdt_interval_owner_selector(
      index.represented_owner_order, index.owner_query_position,
      index.represented_prefix[static_cast<std::size_t>(begin)],
      index.represented_prefix[static_cast<std::size_t>(end)], begin, end,
      index.live_atom_prefix[static_cast<std::size_t>(end)] -
          index.live_atom_prefix[static_cast<std::size_t>(begin)]);
}

void hdt_traverse_active_splits(
    HdtState state,
    const HatQueryTree& query,
    const std::vector<int>& query_owner_order,
    const std::vector<unsigned char>& score_node,
    HdtExtractWorkspace& workspace,
    HdtTraversalStats& stats,
    const HdtActiveSplitVisitor& visitor,
    int state_reuse_bound_factor) {
  if (state_reuse_bound_factor < 1) {
    throw std::invalid_argument("state reuse bound factor must be positive");
  }
  if (query.root < 0 ||
      static_cast<std::size_t>(query.root) >= query.nodes.size()) {
    return;
  }
  const HdtQueryOwnerIndex owner_index =
      build_hdt_query_owner_index(state, query, query_owner_order);
  if (!state.initial_hdt_accounted) note_hdt(state, stats);
  state.initial_hdt_accounted = false;

  enum class Phase : std::uint8_t {
    Enter,
    ExtractLights,
    AfterHeavy,
    InPlaceLights,
  };
  struct Frame {
    int node = -1;
    int heavy = -1;
    int child_begin = 0;
    int child_end = 0;
    int child_cursor = 0;
    int initial_units = 0;
    Phase phase = Phase::Enter;
    bool shared_state = false;
    bool has_in_place_light = false;
    bool retained_heavy_child = false;
    std::unique_ptr<HdtState> owned_state;
    HdtState* state = nullptr;
  };

  std::vector<Frame> frames;
  frames.reserve(query.nodes.size());
  std::vector<std::pair<int, int>> recolor_ranges;
  long long pending_units = 0;
  long long pending_hdt_nodes = 0;
  auto note_pending = [&] {
    stats.pending_state_units_peak =
        std::max(stats.pending_state_units_peak, pending_units);
    stats.pending_state_hdt_nodes_peak =
        std::max(stats.pending_state_hdt_nodes_peak, pending_hdt_nodes);
  };
  auto push_owned = [&](int node, HdtState next) {
    pending_units += next.components.component_units();
    pending_hdt_nodes += next.hdt.stats.reachable_nodes;
    frames.emplace_back();
    Frame& frame = frames.back();
    frame.node = node;
    frame.owned_state =
        std::make_unique<HdtState>(std::move(next));
    frame.state = frame.owned_state.get();
    note_pending();
  };
  auto push_shared = [&](int node, HdtState& current) {
    frames.emplace_back();
    Frame& frame = frames.back();
    frame.node = node;
    frame.shared_state = true;
    frame.state = &current;
  };
  auto pop_frame = [&] {
    Frame& frame = frames.back();
    if (frame.owned_state != nullptr) {
      pending_units -= frame.state->components.component_units();
      pending_hdt_nodes -= frame.state->hdt.stats.reachable_nodes;
    }
    frames.pop_back();
  };
  push_owned(query.root, std::move(state));
  while (!frames.empty()) {
    Frame& frame = frames.back();
    HdtState& current = *frame.state;
    if (frame.phase == Phase::Enter) {
      ++stats.states_visited;
      stats.max_state_units = std::max(
          stats.max_state_units, current.components.component_units());
      if (active_should_score(query, score_node, frame.node)) {
        const StateCounters before = state_counters(current);
        visitor(frame.node, current);
        add_counter_delta(stats, current, before);
        ++stats.score_reads;
      }
      if (is_query_leaf(query, frame.node)) {
        if (frame.shared_state) {
          recolor_query_node(current, query, query_owner_order, frame.node, 0,
                             stats);
        }
        pop_frame();
        continue;
      }
      frame.heavy =
          query.heavy_child[static_cast<std::size_t>(frame.node)];
      frame.child_begin =
          query.child_first[static_cast<std::size_t>(frame.node)];
      frame.child_end = frame.child_begin +
                        query.child_count[static_cast<std::size_t>(frame.node)];
      frame.child_cursor = frame.child_begin;
      frame.initial_units = current.components.component_units();
      frame.phase = Phase::ExtractLights;
      continue;
    }

    if (frame.phase == Phase::ExtractLights) {
      bool descended = false;
      while (frame.child_cursor < frame.child_end) {
        const int child = query.child_nodes[
            static_cast<std::size_t>(frame.child_cursor++)];
        if (child == frame.heavy ||
            (is_query_leaf(query, child) &&
             !active_should_score(query, score_node, child))) {
          continue;
        }
        const HdtLiveSelector selector =
            hdt_query_node_selector(owner_index, query, child);
        if (can_reuse_state_for_child(
                current, selector, frame.initial_units,
                state_reuse_bound_factor)) {
          frame.has_in_place_light = true;
          continue;
        }
        HdtExtractContractResult extracted =
            hdt_extract_contract(current, selector, workspace);
        note_extract(extracted, stats);
        push_owned(child, std::move(extracted.state));
        descended = true;
        break;
      }
      if (descended) continue;

      recolor_light_children(current, query, query_owner_order, frame.node,
                             frame.heavy, stats, recolor_ranges);
      frame.phase = Phase::AfterHeavy;
      const bool visit_heavy =
          frame.heavy >= 0 &&
          (!is_query_leaf(query, frame.heavy) ||
           active_should_score(query, score_node, frame.heavy));
      if (!visit_heavy) {
        if (frame.heavy >= 0 &&
            (frame.shared_state || frame.has_in_place_light)) {
          recolor_query_node(current, query, query_owner_order, frame.heavy, 0,
                             stats);
        }
        continue;
      }

      const HdtLiveSelector selector =
          hdt_query_node_selector(owner_index, query, frame.heavy);
      const int live = hdt_selected_atom_count(current, selector);
      const int bound = hdt_contract_size_bound(live);
      const int input_units = current.components.component_units();
      if ((frame.shared_state || frame.has_in_place_light) &&
          static_cast<long long>(input_units) >
              static_cast<long long>(state_reuse_bound_factor) * bound) {
        HdtExtractContractResult contracted =
            hdt_extract_contract(current, selector, workspace);
        ++stats.heavy_calls;
        stats.heavy_input_units += input_units;
        stats.heavy_output_units += contracted.contract.output_units;
        if (2LL * static_cast<long long>(contracted.contract.output_units) >=
            input_units) {
          ++stats.heavy_shrink_violations;
          stats.complexity_contract_satisfied = false;
          throw std::logic_error(
              "retained heavy Contract did not halve R");
        }
        note_extract(contracted, stats, true);
        frame.retained_heavy_child = true;
        push_owned(frame.heavy, std::move(contracted.state));
      } else {
        if (!frame.has_in_place_light) {
          maybe_heavy_contract(current, selector, workspace, stats,
                               state_reuse_bound_factor);
        }
        push_shared(frame.heavy, current);
      }
      continue;
    }

    if (frame.phase == Phase::AfterHeavy) {
      if (frame.retained_heavy_child) {
        recolor_query_node(current, query, query_owner_order, frame.heavy, 0,
                           stats);
        frame.retained_heavy_child = false;
      }
      if (!frame.has_in_place_light) {
        pop_frame();
        continue;
      }
      frame.child_cursor = frame.child_begin;
      frame.phase = Phase::InPlaceLights;
      continue;
    }

    bool descended = false;
    while (frame.child_cursor < frame.child_end) {
      const int child =
          query.child_nodes[static_cast<std::size_t>(frame.child_cursor++)];
      if (child == frame.heavy ||
          (is_query_leaf(query, child) &&
           !active_should_score(query, score_node, child))) {
        continue;
      }
      const HdtLiveSelector selector =
          hdt_query_node_selector(owner_index, query, child);
      if (!can_reuse_state_for_child(
              current, selector, frame.initial_units,
              state_reuse_bound_factor)) {
        continue;
      }
      recolor_query_node(current, query, query_owner_order, child, 1, stats);
      push_shared(child, current);
      descended = true;
      break;
    }
    if (!descended) pop_frame();
  }
  finalize_memory_stats(workspace, stats);
}

void hdt_traverse_one_x_medians(
    HdtState state,
    const HatQueryTree& query,
    const std::vector<int>& query_owner_order,
    HdtExtractWorkspace& workspace,
    HdtTraversalStats& stats,
    const HdtOneXVisitor& visitor,
    int state_reuse_bound_factor) {
  if (state_reuse_bound_factor < 1) {
    throw std::invalid_argument("state reuse bound factor must be positive");
  }
  if (query.root < 0 ||
      static_cast<std::size_t>(query.root) >= query.nodes.size()) {
    return;
  }
  const HdtQueryOwnerIndex owner_index =
      build_hdt_query_owner_index(state, query, query_owner_order);
  if (!state.initial_hdt_accounted) note_hdt(state, stats);
  state.initial_hdt_accounted = false;

  const HdtLiveSelector root_selector =
      hdt_query_node_selector(owner_index, query, query.root);
  const int root_live = hdt_selected_atom_count(state, root_selector);
  const int root_bound = hdt_contract_size_bound(root_live);
  if (root_live != hdt_live_atom_count(state) ||
      state.components.component_units() > root_bound) {
    HdtContractResult root_contract = hdt_contract(
        state, root_selector, HdtContractMode::LightCopy);
    note_contract(root_contract, stats);
    state = std::move(root_contract.state);
  }

  enum class Phase : std::uint8_t {
    Enter,
    ExtractLights,
    AfterHeavy,
    InPlaceLights,
  };
  struct Frame {
    int node = -1;
    int heavy = -1;
    int child_begin = 0;
    int child_end = 0;
    int child_cursor = 0;
    int initial_units = 0;
    Phase phase = Phase::Enter;
    bool shared_state = false;
    bool has_in_place_light = false;
    bool retained_heavy_child = false;
    std::unique_ptr<HdtState> owned_state;
    HdtState* state = nullptr;
  };

  std::vector<Frame> frames;
  frames.reserve(query.nodes.size());
  std::vector<std::pair<int, int>> recolor_ranges;
  long long pending_units = 0;
  long long pending_hdt_nodes = 0;
  auto note_pending = [&] {
    stats.pending_state_units_peak =
        std::max(stats.pending_state_units_peak, pending_units);
    stats.pending_state_hdt_nodes_peak =
        std::max(stats.pending_state_hdt_nodes_peak, pending_hdt_nodes);
  };
  auto push_owned = [&](int node, HdtState next) {
    pending_units += next.components.component_units();
    pending_hdt_nodes += next.hdt.stats.reachable_nodes;
    frames.emplace_back();
    Frame& frame = frames.back();
    frame.node = node;
    frame.owned_state =
        std::make_unique<HdtState>(std::move(next));
    frame.state = frame.owned_state.get();
    note_pending();
  };
  auto push_shared = [&](int node, HdtState& current) {
    frames.emplace_back();
    Frame& frame = frames.back();
    frame.node = node;
    frame.shared_state = true;
    frame.state = &current;
  };
  auto pop_frame = [&] {
    Frame& frame = frames.back();
    if (frame.owned_state != nullptr) {
      pending_units -= frame.state->components.component_units();
      pending_hdt_nodes -= frame.state->hdt.stats.reachable_nodes;
    }
    frames.pop_back();
  };
  push_owned(query.root, std::move(state));
  while (!frames.empty()) {
    Frame& frame = frames.back();
    HdtState& current = *frame.state;
    if (frame.phase == Phase::Enter) {
      ++stats.states_visited;
      stats.max_state_units = std::max(
          stats.max_state_units, current.components.component_units());
      if (is_query_leaf(query, frame.node)) {
        if (frame.shared_state) {
          recolor_query_node(current, query, query_owner_order, frame.node, 0,
                             stats);
        }
        pop_frame();
        continue;
      }

      frame.child_begin =
          query.child_first[static_cast<std::size_t>(frame.node)];
      frame.child_end = frame.child_begin +
                        query.child_count[static_cast<std::size_t>(frame.node)];
      const int child_count = frame.child_end - frame.child_begin;
      if (frame.node == query.root && child_count == 3) {
        const int child0 = query.child_nodes[
            static_cast<std::size_t>(frame.child_begin)];
        const int child1 = query.child_nodes[
            static_cast<std::size_t>(frame.child_begin + 1)];
        const int child2 = query.child_nodes[
            static_cast<std::size_t>(frame.child_begin + 2)];
        recolor_query_node(current, query, query_owner_order, child0, 0, stats);
        recolor_query_node(current, query, query_owner_order, child2, 2, stats);
        const StateCounters before = state_counters(current);
        visitor(frame.node, std::array<int, 3>{child0, child1, child2},
                std::array<std::uint8_t, 3>{0, 1, 2}, current);
        add_counter_delta(stats, current, before);
        ++stats.score_reads;
        recolor_query_node(current, query, query_owner_order, child0, 1, stats);
        recolor_query_node(current, query, query_owner_order, child2, 1, stats);
      } else if (frame.node != query.root && child_count == 2) {
        const HatMedianBranchColoring coloring =
            hat_binary_median_branch_coloring(query, frame.node);
        recolor_query_node(current, query, query_owner_order,
                           coloring.temporary_color_two_node, 2, stats);
        const StateCounters before = state_counters(current);
        visitor(frame.node, coloring.branch_nodes, coloring.branch_colors,
                current);
        add_counter_delta(stats, current, before);
        ++stats.score_reads;
        recolor_query_node(current, query, query_owner_order,
                           coloring.temporary_color_two_node, 1, stats);
      }

      frame.heavy =
          query.heavy_child[static_cast<std::size_t>(frame.node)];
      frame.child_cursor = frame.child_begin;
      frame.initial_units = current.components.component_units();
      frame.phase = Phase::ExtractLights;
      continue;
    }

    if (frame.phase == Phase::ExtractLights) {
      bool descended = false;
      while (frame.child_cursor < frame.child_end) {
        const int child = query.child_nodes[
            static_cast<std::size_t>(frame.child_cursor++)];
        if (child == frame.heavy || is_query_leaf(query, child)) continue;
        const HdtLiveSelector selector =
            hdt_query_node_selector(owner_index, query, child);
        if (can_reuse_state_for_child(
                current, selector, frame.initial_units,
                state_reuse_bound_factor)) {
          frame.has_in_place_light = true;
          continue;
        }
        HdtExtractContractResult extracted =
            hdt_extract_contract(current, selector, workspace);
        note_extract(extracted, stats);
        push_owned(child, std::move(extracted.state));
        descended = true;
        break;
      }
      if (descended) continue;

      recolor_light_children(current, query, query_owner_order, frame.node,
                             frame.heavy, stats, recolor_ranges);
      frame.phase = Phase::AfterHeavy;
      if (frame.heavy < 0 || is_query_leaf(query, frame.heavy)) {
        if (frame.heavy >= 0 &&
            (frame.shared_state || frame.has_in_place_light)) {
          recolor_query_node(current, query, query_owner_order, frame.heavy, 0,
                             stats);
        }
        continue;
      }

      const HdtLiveSelector selector =
          hdt_query_node_selector(owner_index, query, frame.heavy);
      const int live = hdt_selected_atom_count(current, selector);
      const int bound = hdt_contract_size_bound(live);
      const int input_units = current.components.component_units();
      if ((frame.shared_state || frame.has_in_place_light) &&
          static_cast<long long>(input_units) >
              static_cast<long long>(state_reuse_bound_factor) * bound) {
        HdtExtractContractResult contracted =
            hdt_extract_contract(current, selector, workspace);
        ++stats.heavy_calls;
        stats.heavy_input_units += input_units;
        stats.heavy_output_units += contracted.contract.output_units;
        if (2LL * static_cast<long long>(contracted.contract.output_units) >=
            input_units) {
          ++stats.heavy_shrink_violations;
          stats.complexity_contract_satisfied = false;
          throw std::logic_error(
              "retained heavy Contract did not halve R");
        }
        note_extract(contracted, stats, true);
        frame.retained_heavy_child = true;
        push_owned(frame.heavy, std::move(contracted.state));
      } else {
        if (!frame.has_in_place_light) {
          maybe_heavy_contract(current, selector, workspace, stats,
                               state_reuse_bound_factor);
        }
        push_shared(frame.heavy, current);
      }
      continue;
    }

    if (frame.phase == Phase::AfterHeavy) {
      if (frame.retained_heavy_child) {
        recolor_query_node(current, query, query_owner_order, frame.heavy, 0,
                           stats);
        frame.retained_heavy_child = false;
      }
      if (!frame.has_in_place_light) {
        pop_frame();
        continue;
      }
      frame.child_cursor = frame.child_begin;
      frame.phase = Phase::InPlaceLights;
      continue;
    }

    bool descended = false;
    while (frame.child_cursor < frame.child_end) {
      const int child =
          query.child_nodes[static_cast<std::size_t>(frame.child_cursor++)];
      if (child == frame.heavy || is_query_leaf(query, child)) continue;
      const HdtLiveSelector selector =
          hdt_query_node_selector(owner_index, query, child);
      if (!can_reuse_state_for_child(
              current, selector, frame.initial_units,
              state_reuse_bound_factor)) {
        continue;
      }
      recolor_query_node(current, query, query_owner_order, child, 1, stats);
      push_shared(child, current);
      descended = true;
      break;
    }
    if (!descended) pop_frame();
  }
  finalize_memory_stats(workspace, stats);
}

}  // namespace wqfm::internal
