#include "hat/hdt/hdt_contract.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <span>
#include <stdexcept>
#include <utility>

namespace wqfm::internal {
namespace {

struct GraphTraversal {
  std::vector<int> parent;
  std::vector<int> parent_edge;
  std::vector<int> preorder;
  std::size_t construction_scratch_peak_bytes = 0;
};

template <typename T>
std::size_t vector_capacity_bytes(const std::vector<T>& values,
                                  const char* message) {
  return checked_hdt_size_mul(values.capacity(), sizeof(T), message);
}

void note_ephemeral_scratch(HdtContractStats& stats,
                            std::size_t bytes) {
  stats.ephemeral_scratch_peak_bytes =
      std::max(stats.ephemeral_scratch_peak_bytes, bytes);
}

void finalize_contract_scratch(HdtContractStats& stats,
                               std::size_t persistent_bytes) {
  stats.persistent_scratch_bytes = persistent_bytes;
  stats.scratch_peak_bytes = checked_hdt_size_add(
      persistent_bytes, stats.ephemeral_scratch_peak_bytes,
      "Contract scratch byte count exceeds size_t");
}

int edge_other(const HdtIncidenceEdge& edge, int unit) {
  return edge.a == unit ? edge.b : edge.a;
}

HdtPortKey edge_port(const HdtIncidenceEdge& edge, int unit) {
  return edge.a == unit ? edge.a_port : edge.b_port;
}

bool is_selected_contract_anchor(const HdtComponentUnit& unit,
                                 const HdtLiveSelector& selector) {
  return unit.kind == HdtComponentKind::IAnchor ||
         (unit.kind == HdtComponentKind::LiveLeaf &&
          selector.contains_owner(unit.owner));
}

// Base component graphs inherit the source rooting.  In a sparse bank that
// root can be an absent leaf, represented by a zero-mass C1Opaque unit.  Such a
// unit cannot be retained as a skeleton anchor because its payload belongs to
// the background quotient.  Follow only the unique root tail to the nearest
// structural anchor.  Choosing an arbitrary anchor elsewhere would reverse a
// source-oriented ordinary edge and invalidate the CC transition contract.
int effective_contract_root(const HdtComponentGraph& graph,
                            const HdtLiveSelector& selector) {
  int previous = -1;
  int current = graph.root_unit;
  for (int steps = 0; steps < graph.component_units(); ++steps) {
    if (current < 0 || current >= graph.component_units()) {
      throw std::logic_error("Contract input has an invalid root tail");
    }
    const HdtComponentUnit& unit =
        graph.units[static_cast<std::size_t>(current)];
    if (is_selected_contract_anchor(unit, selector)) return current;

    int next = -1;
    for (int i = 0; i < unit.incidence_count; ++i) {
      const int edge_index = graph.incidence_edges[static_cast<std::size_t>(
          unit.incidence_begin + i)];
      const int candidate = edge_other(
          graph.edges[static_cast<std::size_t>(edge_index)], current);
      if (candidate == previous) continue;
      if (next >= 0) {
        throw std::logic_error(
            "Contract opaque root is not a unique structural tail");
      }
      next = candidate;
    }
    if (next < 0) {
      throw std::logic_error(
          "Contract opaque root tail contains no selected anchor");
    }
    previous = current;
    current = next;
  }
  throw std::logic_error("Contract opaque root tail contains a cycle");
}

GraphTraversal root_component_graph(const HdtComponentGraph& graph,
                                    int root_unit) {
  GraphTraversal traversal;
  const int n = graph.component_units();
  traversal.parent.assign(static_cast<std::size_t>(n), -2);
  traversal.parent_edge.assign(static_cast<std::size_t>(n), -1);
  if (root_unit < 0) return traversal;
  traversal.parent[static_cast<std::size_t>(root_unit)] = -1;
  std::vector<int> stack{root_unit};
  while (!stack.empty()) {
    const int unit = stack.back();
    stack.pop_back();
    traversal.preorder.push_back(unit);
    const HdtComponentUnit& component =
        graph.units[static_cast<std::size_t>(unit)];
    for (int i = component.incidence_count - 1; i >= 0; --i) {
      const int edge_index = graph.incidence_edges[static_cast<std::size_t>(
          component.incidence_begin + i)];
      const int next = edge_other(
          graph.edges[static_cast<std::size_t>(edge_index)], unit);
      if (traversal.parent[static_cast<std::size_t>(next)] != -2) continue;
      traversal.parent[static_cast<std::size_t>(next)] = unit;
      traversal.parent_edge[static_cast<std::size_t>(next)] = edge_index;
      stack.push_back(next);
    }
  }
  if (static_cast<int>(traversal.preorder.size()) != n) {
    throw std::logic_error("Contract input graph is disconnected");
  }
  traversal.construction_scratch_peak_bytes = vector_capacity_bytes(
      stack, "Contract traversal scratch byte count exceeds size_t");
  return traversal;
}

HdtComponentUnit background_unit(const HdtComponentUnit& source) {
  HdtComponentUnit unit = source;
  unit.owner = -1;
  if (unit.kind == HdtComponentKind::LiveLeaf) {
    unit.kind = HdtComponentKind::C1Opaque;
  }
  unit.current_color = unit.background_color;
  unit.current_color_mask = unit.background_color_mask;
  unit.dp_offset = -1;
  unit.background_dp_offset = -1;
  unit.dp_dimension = 0;
  unit.dp_orientation = FullBdpOrientation::None;
  unit.current_equals_background = true;
  if (unit.local_background_g != nullptr &&
      !unit.local_background_g->current_equals_background) {
    FullBdpDual frozen = *unit.local_background_g;
    frozen.current = frozen.background;
    frozen.current_equals_background = true;
    unit.local_background_g =
        std::make_shared<const FullBdpDual>(std::move(frozen));
  }
  unit.incidence_begin = 0;
  unit.incidence_count = 0;
  return unit;
}

struct FoldedFragment {
  FullBdpDual dp;
  HdtOutputKind output_kind = HdtOutputKind::Closed;
  HdtContractStats stats;
};

FoldedFragment fold_background_fragment(
    const HdtComponentGraph& source,
    const std::vector<int>& selected_units,
    std::span<const int> boundary_edges,
    std::vector<int>& old_to_new) {
  if (selected_units.empty()) {
    throw std::logic_error("cannot fold an empty component fragment");
  }
  HdtComponentGraph fragment;
  fragment.palette = source.palette;
  fragment.source = source.source;
  fragment.owner_count = 0;
  fragment.background_restricted = source.background_restricted;
  fragment.units.reserve(selected_units.size());
  std::size_t payload_dimensions = 0;
  for (int old_unit : selected_units) {
    payload_dimensions = checked_hdt_size_add(
        payload_dimensions,
        source.units[static_cast<std::size_t>(old_unit)].dp_dimension,
        "fragment payload capacity overflow");
  }
  (void)checked_hdt_int_size(
      payload_dimensions,
      "fragment payload exceeds its integer offset domain");
  fragment.current_payload.reserve(payload_dimensions);
  fragment.background_payload.reserve(payload_dimensions);
  for (int old_unit : selected_units) {
    old_to_new[static_cast<std::size_t>(old_unit)] =
        checked_hdt_int_size(fragment.units.size(),
                             "fragment exceeds integer node domain");
    const HdtComponentUnit& source_unit =
        source.units[static_cast<std::size_t>(old_unit)];
    FullBdpDual payload;
    const bool has_payload =
        source_unit.kind != HdtComponentKind::IAnchor;
    if (has_payload) {
      payload = hdt_component_dp(source, old_unit);
      payload.current = payload.background;
      payload.current_equals_background = true;
    }
    fragment.units.push_back(background_unit(source_unit));
    if (has_payload) {
      hdt_component_store_dp(fragment, fragment.units.back(), payload);
    }
  }

  for (int old_unit : selected_units) {
    const HdtComponentUnit& component =
        source.units[static_cast<std::size_t>(old_unit)];
    for (int i = 0; i < component.incidence_count; ++i) {
      const int edge_index = source.incidence_edges[static_cast<std::size_t>(
          component.incidence_begin + i)];
      const HdtIncidenceEdge& edge =
          source.edges[static_cast<std::size_t>(edge_index)];
      if (edge.a != old_unit) continue;
      const int mapped_a = old_to_new[static_cast<std::size_t>(edge.a)];
      const int mapped_b = old_to_new[static_cast<std::size_t>(edge.b)];
      if (mapped_a >= 0 && mapped_b >= 0) {
        fragment.edges.push_back(HdtIncidenceEdge{
            mapped_a, mapped_b, edge.a_port, edge.b_port});
      }
    }
  }

  int preferred_root = -1;
  for (int edge_index : boundary_edges) {
    const HdtIncidenceEdge& edge =
        source.edges[static_cast<std::size_t>(edge_index)];
    const bool a_inside = old_to_new[static_cast<std::size_t>(edge.a)] >= 0;
    const bool b_inside = old_to_new[static_cast<std::size_t>(edge.b)] >= 0;
    if (a_inside == b_inside) {
      throw std::logic_error("fragment boundary does not cross its unit set");
    }
    const int old_inside = a_inside ? edge.a : edge.b;
    const HdtPortKey port = a_inside ? edge.a_port : edge.b_port;
    const int mapped = old_to_new[static_cast<std::size_t>(old_inside)];
    fragment.external_ports.push_back(HdtExternalPort{mapped, port});
    if (preferred_root < 0 ||
        port.orientation == HdtPortOrientation::Up) {
      preferred_root = mapped;
    }
  }
  fragment.root_unit = preferred_root >= 0 ? preferred_root : 0;
  rebuild_component_graph_indices(fragment);
  Hdt hdt;
  build_frozen_fold_hdt_into(
      hdt, fragment, boundary_edges.empty());
  const HdtNode& root =
      hdt.nodes[static_cast<std::size_t>(hdt.root)];
  FoldedFragment folded;
  folded.dp = hdt_node_dp(hdt, hdt.root);
  folded.output_kind = root.output_kind;
  folded.stats.fragment_hdt_nodes = hdt.stats.reachable_nodes;
  folded.stats.gg_ops = hdt.stats.gg_ops;
  folded.stats.ig_ops = hdt.stats.ig_ops;
  folded.stats.cc_ops = hdt.stats.cc_ops;
  folded.stats.c_to_g_ops = hdt.stats.c_to_g_ops;
  std::size_t scratch_bytes = 0;
  checked_hdt_accumulate_bytes(
      scratch_bytes, fragment.units.capacity(), sizeof(HdtComponentUnit),
      "Contract fragment byte count exceeds size_t");
  checked_hdt_accumulate_bytes(
      scratch_bytes, fragment.current_payload.capacity(), sizeof(double),
      "Contract fragment byte count exceeds size_t");
  checked_hdt_accumulate_bytes(
      scratch_bytes, fragment.background_payload.capacity(), sizeof(double),
      "Contract fragment byte count exceeds size_t");
  checked_hdt_accumulate_bytes(
      scratch_bytes, fragment.edges.capacity(), sizeof(HdtIncidenceEdge),
      "Contract fragment byte count exceeds size_t");
  checked_hdt_accumulate_bytes(
      scratch_bytes, fragment.external_ports.capacity(),
      sizeof(HdtExternalPort),
      "Contract fragment byte count exceeds size_t");
  checked_hdt_accumulate_bytes(
      scratch_bytes, fragment.incidence_edges.capacity(), sizeof(int),
      "Contract fragment byte count exceeds size_t");
  checked_hdt_accumulate_bytes(
      scratch_bytes, hdt.nodes.capacity(), sizeof(HdtNode),
      "Contract fragment byte count exceeds size_t");
  scratch_bytes = checked_hdt_size_add(
      scratch_bytes, hdt.current_payload.capacity_bytes(),
      "Contract fragment byte count exceeds size_t");
  scratch_bytes = checked_hdt_size_add(
      scratch_bytes, hdt.background_payload.capacity_bytes(),
      "Contract fragment byte count exceeds size_t");
  checked_hdt_accumulate_bytes(
      scratch_bytes, hdt.atoms.capacity(), sizeof(HdtExpandedAtom),
      "Contract fragment byte count exceeds size_t");
  checked_hdt_accumulate_bytes(
      scratch_bytes, hdt.atom_ports.capacity(), sizeof(HdtPortKey),
      "Contract fragment byte count exceeds size_t");
  checked_hdt_accumulate_bytes(
      scratch_bytes, hdt.owner_offsets.capacity(), sizeof(int),
      "Contract fragment byte count exceeds size_t");
  checked_hdt_accumulate_bytes(
      scratch_bytes, hdt.owner_atom_nodes.capacity(), sizeof(int),
      "Contract fragment byte count exceeds size_t");
  checked_hdt_accumulate_bytes(
      scratch_bytes, hdt.recompute_stack.capacity(), sizeof(int),
      "Contract fragment byte count exceeds size_t");
  checked_hdt_accumulate_bytes(
      scratch_bytes, hdt.dirty.capacity(), sizeof(std::uint8_t),
      "Contract fragment byte count exceeds size_t");
  scratch_bytes = checked_hdt_size_add(
      scratch_bytes, hdt.stats.build_scratch_peak_bytes,
      "Contract fragment byte count exceeds size_t");
  folded.stats.scratch_peak_bytes = scratch_bytes;
  for (int old_unit : selected_units) {
    old_to_new[static_cast<std::size_t>(old_unit)] = -1;
  }
  return folded;
}

void add_fold_stats(HdtContractStats& total,
                    const FoldedFragment& folded,
                    std::size_t caller_scratch_bytes) {
  total.fragment_hdt_nodes += folded.stats.fragment_hdt_nodes;
  total.gg_ops += folded.stats.gg_ops;
  total.ig_ops += folded.stats.ig_ops;
  total.cc_ops += folded.stats.cc_ops;
  total.c_to_g_ops += folded.stats.c_to_g_ops;
  note_ephemeral_scratch(
      total,
      checked_hdt_size_add(caller_scratch_bytes,
                           folded.stats.scratch_peak_bytes,
                           "Contract fold scratch byte count exceeds size_t"));
}

struct CollectedBranch {
  std::vector<int> units;
  std::size_t traversal_scratch_peak_bytes = 0;
};

CollectedBranch collect_branch_units(const HdtComponentGraph& graph,
                                     const GraphTraversal& traversal,
                                     int branch_root) {
  CollectedBranch result;
  std::vector<int> stack{branch_root};
  while (!stack.empty()) {
    const int unit = stack.back();
    stack.pop_back();
    result.units.push_back(unit);
    const HdtComponentUnit& component =
        graph.units[static_cast<std::size_t>(unit)];
    for (int i = component.incidence_count - 1; i >= 0; --i) {
      const int edge_index = graph.incidence_edges[static_cast<std::size_t>(
          component.incidence_begin + i)];
      const int next = edge_other(
          graph.edges[static_cast<std::size_t>(edge_index)], unit);
      if (traversal.parent[static_cast<std::size_t>(next)] == unit) {
        stack.push_back(next);
      }
    }
  }
  result.traversal_scratch_peak_bytes = vector_capacity_bytes(
      stack, "Contract branch traversal byte count exceeds size_t");
  return result;
}

FullBdpDual branch_as_g(int palette,
                        const FoldedFragment& branch,
                        int expected_center) {
  (void)expected_center;
  if (branch.output_kind == HdtOutputKind::G) return branch.dp;
  FullBdpDual c1 = branch.dp;
  if (branch.output_kind == HdtOutputKind::C1D) {
    c1 = full_bdp_reverse(palette, c1);
  } else if (branch.output_kind != HdtOutputKind::C1U) {
    throw std::logic_error("background branch did not fold to C1/G");
  }
  return full_bdp_c_to_g(palette, c1);
}

void connect_output_edge(HdtComponentGraph& output,
                         int a,
                         int b,
                         HdtPortKey a_port,
                         HdtPortKey b_port) {
  output.edges.push_back(HdtIncidenceEdge{a, b, a_port, b_port});
}

HdtComponentKind opaque_kind(HdtOutputKind kind) {
  switch (kind) {
    case HdtOutputKind::Closed:
      return HdtComponentKind::ClosedOpaque;
    case HdtOutputKind::C1U:
    case HdtOutputKind::C1D:
      return HdtComponentKind::C1Opaque;
    case HdtOutputKind::C2:
      return HdtComponentKind::C2Opaque;
    case HdtOutputKind::G:
      return HdtComponentKind::GOpaque;
    case HdtOutputKind::I: break;
  }
  throw std::logic_error("cannot make an opaque I component");
}

struct EmittedPartitionUnit {
  int port_begin = 0;
  int port_count = 0;
};

FullBdpDual frozen_dual(const FullBdpVec& background) {
  return FullBdpDual{background, background, true};
}

}  // namespace

int HdtLiveSelector::selected_owner_count() const {
  if (owner_order != nullptr) {
    if (owner_begin < 0 || owner_end < owner_begin) {
      throw std::logic_error("invalid selector owner interval");
    }
    return owner_end - owner_begin;
  }
  return checked_hdt_int_size(
      owned_owners.size(), "selector exceeds integer owner domain");
}

int HdtLiveSelector::selected_owner_at(int index) const {
  if (index < 0 || index >= selected_owner_count()) {
    throw std::out_of_range("selector owner index out of range");
  }
  return owner_order != nullptr
             ? (*owner_order)[static_cast<std::size_t>(owner_begin + index)]
             : owned_owners[static_cast<std::size_t>(index)];
}

bool HdtLiveSelector::contains_owner(int owner) const {
  if (owner < 0) return false;
  if (owner_position != nullptr) {
    if (static_cast<std::size_t>(owner) >= owner_position->size()) return false;
    const int position = (*owner_position)[static_cast<std::size_t>(owner)];
    return position >= membership_begin && position < membership_end;
  }
  return static_cast<std::size_t>(owner) < dense_membership.size() &&
         dense_membership[static_cast<std::size_t>(owner)] != 0;
}

HdtOwnerIndexView owner_index_view(
    const HdtLiveSelector& selector) {
  if (selector.owner_position == nullptr) return {};
  return HdtOwnerIndexView{selector.owner_position,
                                 selector.membership_begin,
                                 selector.membership_end};
}

HdtLiveSelector hdt_owner_selector(
    int owner_count,
    const std::vector<int>& owners) {
  if (owner_count < 0) {
    throw std::invalid_argument("negative owner domain");
  }
  HdtLiveSelector selector;
  selector.dense_membership.assign(
      static_cast<std::size_t>(std::max(0, owner_count)), 0);
  for (int owner : owners) {
    if (owner < 0 || owner >= owner_count ||
        selector.dense_membership[static_cast<std::size_t>(owner)] != 0) {
      continue;
    }
    selector.dense_membership[static_cast<std::size_t>(owner)] = 1;
    selector.owned_owners.push_back(owner);
  }
  return selector;
}

HdtLiveSelector hdt_interval_owner_selector(
    const std::vector<int>& represented_owner_order,
    const std::vector<int>& owner_position,
    int owner_begin,
    int owner_end,
    int membership_begin,
    int membership_end,
    int exact_selected_atom_count) {
  const int represented_owner_count = checked_hdt_int_size(
      represented_owner_order.size(),
      "interval selector exceeds integer owner domain");
  if (owner_begin < 0 || owner_end < owner_begin ||
      owner_end > represented_owner_count || membership_begin < 0 ||
      membership_end < membership_begin || exact_selected_atom_count < -1) {
    throw std::invalid_argument("invalid interval selector");
  }
  HdtLiveSelector selector;
  selector.owner_order = &represented_owner_order;
  selector.owner_position = &owner_position;
  selector.owner_begin = owner_begin;
  selector.owner_end = owner_end;
  selector.membership_begin = membership_begin;
  selector.membership_end = membership_end;
  selector.exact_selected_atom_count = exact_selected_atom_count;
  return selector;
}

int hdt_selected_atom_count(
    const HdtState& state,
    const HdtLiveSelector& selector) {
  if (selector.exact_selected_atom_count >= 0) {
    return selector.exact_selected_atom_count;
  }
  int count = 0;
  const HdtComponentGraph& graph = state.components;
  for (const HdtComponentUnit& unit : graph.units) {
    if (unit.kind == HdtComponentKind::LiveLeaf &&
        selector.contains_owner(unit.owner)) {
      ++count;
    }
  }
  return count;
}

int hdt_contract_size_bound(int selected_live_atoms) {
  if (selected_live_atoms < 0) {
    throw std::invalid_argument("negative selected-atom count");
  }
  return checked_hdt_int_mul_add(
      std::max(1, selected_live_atoms), 4, 3,
      "Contract bound exceeds integer domain");
}

HdtContractResult hdt_contract(
    const HdtState& input,
    const HdtLiveSelector& selector,
    HdtContractMode mode) {
  using Clock = std::chrono::steady_clock;
  const auto contract_started = Clock::now();
  const long long scalar_fmas_before = full_bdp_thread_scalar_fmas();
  const HdtComponentGraph& graph = input.components;
  HdtContractResult result;
  HdtContractStats& stats = result.stats;
  stats.mode = mode;
  stats.input_units = graph.component_units();
  stats.selected_live_atoms =
      hdt_selected_atom_count(input, selector);
  stats.size_bound = hdt_contract_size_bound(stats.selected_live_atoms);
  if (graph.units.empty()) {
    stats.size_bound_ok = true;
    stats.mark_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(
                        Clock::now() - contract_started)
                        .count();
    stats.scalar_fmas =
        full_bdp_thread_scalar_fmas() - scalar_fmas_before;
    return result;
  }

  bool all_live_selected = true;
  int live_atoms = 0;
  for (const HdtComponentUnit& unit : graph.units) {
    if (unit.kind != HdtComponentKind::LiveLeaf) continue;
    ++live_atoms;
    all_live_selected =
        all_live_selected && selector.contains_owner(unit.owner);
  }
  const bool input_root_is_anchor = is_selected_contract_anchor(
      graph.units[static_cast<std::size_t>(graph.root_unit)], selector);
  if (mode == HdtContractMode::LightCopy && all_live_selected &&
      live_atoms == stats.selected_live_atoms &&
      graph.component_units() <= stats.size_bound && input_root_is_anchor) {
    stats.mark_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(
                        Clock::now() - contract_started)
                        .count();
    const auto rebuild_started = Clock::now();
    result.state = build_hdt_state(graph, owner_index_view(selector));
    stats.hdt_rebuild_ns =
        std::chrono::duration_cast<std::chrono::nanoseconds>(
            Clock::now() - rebuild_started)
            .count();
    stats.output_units = graph.component_units();
    stats.size_bound_ok = true;
    note_ephemeral_scratch(
        stats, result.state.hdt.stats.build_scratch_peak_bytes);
    finalize_contract_scratch(stats, 0);
    stats.scalar_fmas =
        full_bdp_thread_scalar_fmas() - scalar_fmas_before;
    return result;
  }

  const int contract_root = stats.selected_live_atoms == 0
                                ? graph.root_unit
                                : effective_contract_root(graph, selector);
  const GraphTraversal traversal = root_component_graph(graph, contract_root);
  note_ephemeral_scratch(stats,
                         traversal.construction_scratch_peak_bytes);
  const int n = graph.component_units();
  std::vector<int> live_count(static_cast<std::size_t>(n), 0);
  for (auto it = traversal.preorder.rbegin(); it != traversal.preorder.rend();
       ++it) {
    const int unit = *it;
    const HdtComponentUnit& component =
        graph.units[static_cast<std::size_t>(unit)];
    int live = component.kind == HdtComponentKind::LiveLeaf &&
                       selector.contains_owner(component.owner)
                   ? 1
                   : 0;
    const HdtComponentUnit& node =
        graph.units[static_cast<std::size_t>(unit)];
    for (int i = 0; i < node.incidence_count; ++i) {
      const int edge_index = graph.incidence_edges[static_cast<std::size_t>(
          node.incidence_begin + i)];
      const int child = edge_other(
          graph.edges[static_cast<std::size_t>(edge_index)], unit);
      if (traversal.parent[static_cast<std::size_t>(child)] == unit) {
        live += live_count[static_cast<std::size_t>(child)];
      }
    }
    live_count[static_cast<std::size_t>(unit)] = live;
  }
  if (live_count[static_cast<std::size_t>(contract_root)] !=
      stats.selected_live_atoms) {
    throw std::logic_error("Contract selected-live count mismatch");
  }

  std::vector<int> old_to_new(static_cast<std::size_t>(n), -1);
  std::vector<int> fragment_map(static_cast<std::size_t>(n), -1);
  if (stats.selected_live_atoms == 0) {
    const auto background_started = Clock::now();
    stats.mark_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(
                        background_started - contract_started)
                        .count();
    std::vector<int> all_units = traversal.preorder;
    const FoldedFragment folded = fold_background_fragment(
        graph, all_units, std::span<const int>{}, fragment_map);
    add_fold_stats(
        stats, folded,
        vector_capacity_bytes(
            all_units, "Contract all-unit byte count exceeds size_t"));
    const auto background_finished = Clock::now();
    stats.background_fold_ns =
        std::chrono::duration_cast<std::chrono::nanoseconds>(
            background_finished - background_started)
            .count();
    HdtComponentGraph output;
    output.palette = graph.palette;
    output.source = graph.source;
    output.owner_count = graph.owner_count;
    output.background_restricted = graph.background_restricted;
    output.root_unit = 0;
    HdtComponentUnit closed;
    closed.kind = HdtComponentKind::ClosedOpaque;
    closed.stable_id = 0;
    closed.source_vertex_id = graph.source ? graph.source->root_vertex : -1;
    output.units.push_back(std::move(closed));
    hdt_component_store_dp(output, output.units.back(), folded.dp);
    const auto emit_started = Clock::now();
    rebuild_component_graph_indices(output);
    stats.emit_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(
                        Clock::now() - emit_started)
                        .count();
    const auto rebuild_started = Clock::now();
    result.state = build_hdt_state(std::move(output),
                                         owner_index_view(selector));
    stats.hdt_rebuild_ns =
        std::chrono::duration_cast<std::chrono::nanoseconds>(
            Clock::now() - rebuild_started)
            .count();
    stats.output_units = 1;
    stats.anchors = 1;
    stats.background_units_absorbed = n;
    stats.size_bound_ok = true;
    note_ephemeral_scratch(
        stats, result.state.hdt.stats.build_scratch_peak_bytes);
    std::size_t main_scratch = 0;
    for (std::size_t capacity : {
             traversal.parent.capacity(), traversal.parent_edge.capacity(),
             traversal.preorder.capacity(), live_count.capacity(),
             old_to_new.capacity(), fragment_map.capacity()}) {
      checked_hdt_accumulate_bytes(
          main_scratch, capacity, sizeof(int),
          "Contract scratch byte count exceeds size_t");
    }
    finalize_contract_scratch(stats, main_scratch);
    stats.scalar_fmas =
        full_bdp_thread_scalar_fmas() - scalar_fmas_before;
    return result;
  }

  std::vector<std::uint8_t> marked_edge(graph.edges.size(), 0);
  std::vector<int> marked_degree(static_cast<std::size_t>(n), 0);
  for (int unit : traversal.preorder) {
    if (unit == contract_root ||
        live_count[static_cast<std::size_t>(unit)] == 0) {
      continue;
    }
    const int edge_index = traversal.parent_edge[static_cast<std::size_t>(unit)];
    marked_edge[static_cast<std::size_t>(edge_index)] = 1;
    const HdtIncidenceEdge& edge =
        graph.edges[static_cast<std::size_t>(edge_index)];
    ++marked_degree[static_cast<std::size_t>(edge.a)];
    ++marked_degree[static_cast<std::size_t>(edge.b)];
  }

  std::vector<std::uint8_t> anchor(static_cast<std::size_t>(n), 0);
  for (int unit : traversal.preorder) {
    if (unit != contract_root &&
        live_count[static_cast<std::size_t>(unit)] == 0) {
      continue;
    }
    const HdtComponentUnit& component =
        graph.units[static_cast<std::size_t>(unit)];
    const bool selected_live =
        component.kind == HdtComponentKind::LiveLeaf &&
        selector.contains_owner(component.owner);
    if (unit == contract_root || selected_live ||
        marked_degree[static_cast<std::size_t>(unit)] != 2) {
      anchor[static_cast<std::size_t>(unit)] = 1;
    }
  }
  stats.mark_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(
                      Clock::now() - contract_started)
                      .count();

  HdtComponentGraph output;
  output.palette = graph.palette;
  output.source = graph.source;
  output.owner_count = graph.owner_count;
  output.background_restricted = graph.background_restricted;
  for (int unit : traversal.preorder) {
    if (anchor[static_cast<std::size_t>(unit)] == 0) continue;
    HdtComponentUnit retained =
        graph.units[static_cast<std::size_t>(unit)];
    if (retained.kind != HdtComponentKind::IAnchor &&
        retained.kind != HdtComponentKind::LiveLeaf) {
      throw std::logic_error(
          "Contract retained a non-anchor component kind=" +
          std::to_string(static_cast<int>(retained.kind)) + " unit=" +
          std::to_string(unit) + " root=" +
          std::to_string(contract_root) + " marked_degree=" +
          std::to_string(marked_degree[static_cast<std::size_t>(unit)]));
    }
    retained.incidence_begin = 0;
    retained.incidence_count = 0;
    retained.dp_offset = -1;
    retained.background_dp_offset = -1;
    if (retained.kind == HdtComponentKind::IAnchor) {
      retained.dp_dimension = 0;
      retained.dp_orientation = FullBdpOrientation::None;
    }
    old_to_new[static_cast<std::size_t>(unit)] =
        checked_hdt_int_size(output.units.size(),
                             "Contract exceeds integer node domain");
    output.units.push_back(std::move(retained));
    ++stats.anchors;
  }
  output.root_unit = old_to_new[static_cast<std::size_t>(contract_root)];

  // Fold every unmarked branch incident to a retained I into one local G value.
  const auto background_started = Clock::now();
  for (int anchor_unit : traversal.preorder) {
    if (anchor[static_cast<std::size_t>(anchor_unit)] == 0) continue;
    HdtComponentUnit& retained = output.units[static_cast<std::size_t>(
        old_to_new[static_cast<std::size_t>(anchor_unit)])];
    if (retained.kind != HdtComponentKind::IAnchor) continue;
    bool have_g = retained.local_background_g != nullptr;
    FullBdpDual local = have_g ? *retained.local_background_g : FullBdpDual{};
    const HdtComponentUnit& original =
        graph.units[static_cast<std::size_t>(anchor_unit)];
    for (int i = 0; i < original.incidence_count; ++i) {
      const int edge_index = graph.incidence_edges[static_cast<std::size_t>(
          original.incidence_begin + i)];
      if (marked_edge[static_cast<std::size_t>(edge_index)] != 0) continue;
      const int branch_root = edge_other(
          graph.edges[static_cast<std::size_t>(edge_index)], anchor_unit);
      if (traversal.parent[static_cast<std::size_t>(branch_root)] != anchor_unit) {
        continue;
      }
      const CollectedBranch branch =
          collect_branch_units(graph, traversal, branch_root);
      const std::size_t branch_units_bytes = vector_capacity_bytes(
          branch.units, "Contract branch byte count exceeds size_t");
      note_ephemeral_scratch(
          stats,
          checked_hdt_size_add(
              branch_units_bytes, branch.traversal_scratch_peak_bytes,
              "Contract branch traversal byte count exceeds size_t"));
      const std::array<int, 1> boundary{edge_index};
      const FoldedFragment folded = fold_background_fragment(
          graph, branch.units, boundary, fragment_map);
      add_fold_stats(stats, folded, branch_units_bytes);
      stats.background_units_absorbed = checked_hdt_int_add(
          stats.background_units_absorbed,
          checked_hdt_int_size(
              branch.units.size(),
              "Contract branch exceeds integer node domain"),
          "Contract absorbed-unit count overflow");
      const FullBdpDual arm = branch_as_g(
          graph.palette, folded, retained.source_vertex_id);
      ++stats.polytomy_arms_folded;
      local = have_g ? full_bdp_gg(graph.palette, local, arm) : arm;
      have_g = true;
    }
    if (have_g) {
      local.current = local.background;
      local.current_equals_background = true;
      retained.local_background_g =
          std::make_shared<const FullBdpDual>(std::move(local));
    } else {
      retained.local_background_g.reset();
    }
  }
  stats.background_fold_ns =
      std::chrono::duration_cast<std::chrono::nanoseconds>(
          Clock::now() - background_started)
          .count();

  // Emit one direct edge or one folded C2 unit for every anchor-to-anchor path.
  const auto corridor_started = Clock::now();
  for (int top_anchor : traversal.preorder) {
    if (anchor[static_cast<std::size_t>(top_anchor)] == 0) continue;
    const HdtComponentUnit& top_component =
        graph.units[static_cast<std::size_t>(top_anchor)];
    for (int i = 0; i < top_component.incidence_count; ++i) {
      const int first_edge = graph.incidence_edges[static_cast<std::size_t>(
          top_component.incidence_begin + i)];
      if (marked_edge[static_cast<std::size_t>(first_edge)] == 0) continue;
      const int first_child = edge_other(
          graph.edges[static_cast<std::size_t>(first_edge)], top_anchor);
      if (traversal.parent[static_cast<std::size_t>(first_child)] != top_anchor) {
        continue;
      }
      int previous = top_anchor;
      int current = first_child;
      int last_edge = first_edge;
      std::vector<int> corridor_units;
      while (anchor[static_cast<std::size_t>(current)] == 0) {
        corridor_units.push_back(current);
        int next = -1;
        const HdtComponentUnit& component =
            graph.units[static_cast<std::size_t>(current)];
        for (int q = 0; q < component.incidence_count; ++q) {
          const int edge_index = graph.incidence_edges[static_cast<std::size_t>(
              component.incidence_begin + q)];
          if (marked_edge[static_cast<std::size_t>(edge_index)] == 0) continue;
          const int candidate = edge_other(
              graph.edges[static_cast<std::size_t>(edge_index)], current);
          if (candidate == previous) continue;
          next = candidate;
          last_edge = edge_index;
          break;
        }
        if (next < 0) {
          throw std::logic_error("Contract corridor ended before an anchor");
        }
        previous = current;
        current = next;
      }
      const int bottom_anchor = current;
      if (corridor_units.empty()) {
        const HdtIncidenceEdge& edge =
            graph.edges[static_cast<std::size_t>(first_edge)];
        connect_output_edge(
            output, old_to_new[static_cast<std::size_t>(top_anchor)],
            old_to_new[static_cast<std::size_t>(bottom_anchor)],
            edge_port(edge, top_anchor), edge_port(edge, bottom_anchor));
        continue;
      }

      std::vector<int> fragment_units = corridor_units;
      for (int path_unit : corridor_units) {
        const HdtComponentUnit& component =
            graph.units[static_cast<std::size_t>(path_unit)];
        for (int q = 0; q < component.incidence_count; ++q) {
          const int edge_index = graph.incidence_edges[static_cast<std::size_t>(
              component.incidence_begin + q)];
          if (marked_edge[static_cast<std::size_t>(edge_index)] != 0) continue;
          const int branch_root = edge_other(
              graph.edges[static_cast<std::size_t>(edge_index)], path_unit);
          if (traversal.parent[static_cast<std::size_t>(branch_root)] != path_unit) {
            continue;
          }
          const CollectedBranch branch =
              collect_branch_units(graph, traversal, branch_root);
          const std::size_t corridor_bytes = vector_capacity_bytes(
              corridor_units, "Contract corridor byte count exceeds size_t");
          const std::size_t fragment_before_bytes = vector_capacity_bytes(
              fragment_units,
              "Contract corridor fragment byte count exceeds size_t");
          const std::size_t branch_units_bytes = vector_capacity_bytes(
              branch.units, "Contract branch byte count exceeds size_t");
          note_ephemeral_scratch(
              stats,
              checked_hdt_size_add(
                  checked_hdt_size_add(
                      corridor_bytes, fragment_before_bytes,
                      "Contract corridor assembly byte count exceeds size_t"),
                  checked_hdt_size_add(
                      branch_units_bytes,
                      branch.traversal_scratch_peak_bytes,
                      "Contract branch traversal byte count exceeds size_t"),
                  "Contract corridor assembly byte count exceeds size_t"));
          fragment_units.insert(fragment_units.end(), branch.units.begin(),
                                branch.units.end());
          note_ephemeral_scratch(
              stats,
              checked_hdt_size_add(
                  checked_hdt_size_add(
                      corridor_bytes,
                      vector_capacity_bytes(
                          fragment_units,
                          "Contract corridor fragment byte count exceeds size_t"),
                      "Contract corridor assembly byte count exceeds size_t"),
                  branch_units_bytes,
                  "Contract corridor assembly byte count exceeds size_t"));
          ++stats.polytomy_arms_folded;
        }
      }
      const std::array<int, 2> boundary{first_edge, last_edge};
      const FoldedFragment folded = fold_background_fragment(
          graph, fragment_units, boundary, fragment_map);
      const std::size_t corridor_fold_caller_bytes = checked_hdt_size_add(
          vector_capacity_bytes(
              corridor_units, "Contract corridor byte count exceeds size_t"),
          vector_capacity_bytes(
              fragment_units,
              "Contract corridor fragment byte count exceeds size_t"),
          "Contract corridor fold byte count exceeds size_t");
      add_fold_stats(stats, folded, corridor_fold_caller_bytes);
      if (folded.output_kind != HdtOutputKind::C2) {
        throw std::logic_error("Contract corridor did not fold to C2");
      }
      HdtComponentUnit corridor;
      corridor.kind = HdtComponentKind::C2Opaque;
      corridor.stable_id =
          (static_cast<std::uint64_t>(top_anchor) << 32u) |
          static_cast<std::uint32_t>(bottom_anchor);
      corridor.owner = -1;
      const int corridor_index = checked_hdt_int_size(
          output.units.size(),
          "Contract exceeds integer node domain");
      output.units.push_back(std::move(corridor));
      hdt_component_store_dp(output, output.units.back(), folded.dp);
      const HdtIncidenceEdge& top_edge =
          graph.edges[static_cast<std::size_t>(first_edge)];
      const HdtIncidenceEdge& bottom_edge =
          graph.edges[static_cast<std::size_t>(last_edge)];
      const int top_inside = edge_other(top_edge, top_anchor);
      const int bottom_inside = edge_other(bottom_edge, bottom_anchor);
      connect_output_edge(
          output, old_to_new[static_cast<std::size_t>(top_anchor)],
          corridor_index, edge_port(top_edge, top_anchor),
          edge_port(top_edge, top_inside));
      connect_output_edge(
          output, corridor_index,
          old_to_new[static_cast<std::size_t>(bottom_anchor)],
          edge_port(bottom_edge, bottom_inside),
          edge_port(bottom_edge, bottom_anchor));
      ++stats.corridors_emitted;
      stats.background_units_absorbed = checked_hdt_int_add(
          stats.background_units_absorbed,
          checked_hdt_int_size(
              fragment_units.size(),
              "Contract fragment exceeds integer node domain"),
          "Contract absorbed-unit count overflow");
    }
  }
  stats.corridor_fold_ns =
      std::chrono::duration_cast<std::chrono::nanoseconds>(
          Clock::now() - corridor_started)
          .count();

  const auto emit_started = Clock::now();
  rebuild_component_graph_indices(output);
#ifndef NDEBUG
  std::string error;
  if (!validate_component_graph(output, &error)) {
    throw std::logic_error("Contract output invalid: " + error);
  }
#endif
  stats.output_units = output.component_units();
  stats.size_bound_ok = stats.output_units <= stats.size_bound;
  if (!stats.size_bound_ok) {
    throw std::logic_error("Contract output exceeded 4*live+3 bound");
  }
  stats.emit_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(
                      Clock::now() - emit_started)
                      .count();
  const auto rebuild_started = Clock::now();
  result.state = build_hdt_state(std::move(output),
                                       owner_index_view(selector));
  stats.hdt_rebuild_ns =
      std::chrono::duration_cast<std::chrono::nanoseconds>(
          Clock::now() - rebuild_started)
          .count();
  note_ephemeral_scratch(
      stats, result.state.hdt.stats.build_scratch_peak_bytes);
  std::size_t main_scratch = 0;
  for (std::size_t capacity : {
           traversal.parent.capacity(), traversal.parent_edge.capacity(),
           traversal.preorder.capacity(), live_count.capacity(),
           old_to_new.capacity(), fragment_map.capacity(),
           marked_degree.capacity()}) {
    checked_hdt_accumulate_bytes(
        main_scratch, capacity, sizeof(int),
        "Contract scratch byte count exceeds size_t");
  }
  checked_hdt_accumulate_bytes(
      main_scratch, marked_edge.capacity(), sizeof(std::uint8_t),
      "Contract scratch byte count exceeds size_t");
  checked_hdt_accumulate_bytes(
      main_scratch, anchor.capacity(), sizeof(std::uint8_t),
      "Contract scratch byte count exceeds size_t");
  finalize_contract_scratch(stats, main_scratch);
  stats.scalar_fmas =
      full_bdp_thread_scalar_fmas() - scalar_fmas_before;
  return result;
}

HdtInitialStateResult hdt_build_initial_state(
    HdtComponentGraph graph) {
  HdtInitialStateResult initial;
  const bool has_absent_leaf =
      std::any_of(graph.units.begin(), graph.units.end(),
                  [](const HdtComponentUnit& unit) {
                    return unit.kind == HdtComponentKind::C1Opaque &&
                           unit.weight == 0.0;
                  });
  if (!has_absent_leaf) {
    initial.state = build_hdt_state(std::move(graph));
    return initial;
  }

  const long long construction_scalar_fmas = graph.construction_scalar_fmas;
  std::vector<int> owners;
  owners.reserve(graph.units.size());
  for (const HdtComponentUnit& unit : graph.units) {
    if (unit.kind == HdtComponentKind::LiveLeaf) {
      owners.push_back(unit.owner);
    }
  }
  const HdtLiveSelector selector =
      hdt_owner_selector(graph.owner_count, owners);
  HdtState unbuilt;
  unbuilt.components = std::move(graph);
  HdtContractResult compacted = hdt_contract(
      unbuilt, selector, HdtContractMode::InitialCompact);
  compacted.stats.scalar_fmas += construction_scalar_fmas;
  initial.state = std::move(compacted.state);
  initial.contract = compacted.stats;
  initial.contracted = true;
  return initial;
}

HdtExtractContractResult hdt_extract_contract(
    const HdtState& parent,
    const HdtLiveSelector& selector,
    HdtExtractWorkspace& workspace) {
  HdtExtractContractResult result;
  const Hdt& hdt = parent.hdt;
  result.extract.source_hdt_nodes = hdt.stats.reachable_nodes;
  result.extract.selected_live_atoms =
      hdt_selected_atom_count(parent, selector);
  if (hdt.root < 0) return result;

  ++workspace.epoch;
  if (workspace.epoch == 0) {
    std::fill(workspace.hdt_mark_epoch.begin(),
              workspace.hdt_mark_epoch.end(), 0);
    std::fill(workspace.ordinary_epoch.begin(),
              workspace.ordinary_epoch.end(), 0);
    std::fill(workspace.center_epoch.begin(), workspace.center_epoch.end(), 0);
    ++workspace.epoch;
  }
  const std::uint64_t epoch = workspace.epoch;
  workspace.hdt_mark_epoch.resize(hdt.nodes.size(), 0);
  auto mark_to_root = [&](int node) {
    while (node >= 0 &&
           workspace.hdt_mark_epoch[static_cast<std::size_t>(node)] != epoch) {
      workspace.hdt_mark_epoch[static_cast<std::size_t>(node)] = epoch;
      ++result.extract.marked_hdt_nodes;
      node = hdt.nodes[static_cast<std::size_t>(node)].parent;
    }
  };
  if (hdt.owner_index_by_position && selector.owner_position != nullptr) {
    if (hdt.owner_offsets.empty()) {
      throw std::logic_error("owner-position index has no offsets");
    }
    const int indexed_slots = checked_hdt_int_size(
        hdt.owner_offsets.size() - 1u,
        "owner index exceeds integer domain");
    const int indexed_end = checked_hdt_int_add(
        hdt.owner_position_begin, indexed_slots,
        "owner-position interval overflow");
    const int begin_position =
        std::max(hdt.owner_position_begin, selector.membership_begin);
    const int end_position =
        std::min(indexed_end, selector.membership_end);
    for (int position = begin_position; position < end_position; ++position) {
      const int slot = position - hdt.owner_position_begin;
      const int begin = hdt.owner_offsets[static_cast<std::size_t>(slot)];
      const int end = hdt.owner_offsets[static_cast<std::size_t>(slot) + 1u];
      for (int p = begin; p < end; ++p) {
        mark_to_root(hdt.owner_atom_nodes[static_cast<std::size_t>(p)]);
      }
    }
  } else if (hdt.owner_index_by_position) {
    for (int node : hdt.owner_atom_nodes) {
      const HdtNode& atom_node =
          hdt.nodes[static_cast<std::size_t>(node)];
      const HdtExpandedAtom& atom = hdt.atoms[static_cast<std::size_t>(
          atom_node.atom_descriptor)];
      if (selector.contains_owner(atom.owner)) mark_to_root(node);
    }
  } else {
    for (int i = 0; i < selector.selected_owner_count(); ++i) {
      const int owner = selector.selected_owner_at(i);
      const int owner_offset_count = checked_hdt_int_size(
          hdt.owner_offsets.size(),
          "owner index exceeds integer domain");
      if (owner < 0 || owner >= owner_offset_count - 1) {
        continue;
      }
      const int begin = hdt.owner_offsets[static_cast<std::size_t>(owner)];
      const int end =
          hdt.owner_offsets[static_cast<std::size_t>(owner) + 1u];
      for (int p = begin; p < end; ++p) {
        mark_to_root(hdt.owner_atom_nodes[static_cast<std::size_t>(p)]);
      }
    }
  }
  if (hdt.root_anchor_atom_node >= 0) mark_to_root(hdt.root_anchor_atom_node);

  const std::size_t emitted_bound = std::max<std::size_t>(
      1u, checked_hdt_size_mul(
              static_cast<std::size_t>(result.extract.marked_hdt_nodes), 2u,
              "Extract emission bound overflow"));
  const std::size_t retained_bound = static_cast<std::size_t>(
      hdt_contract_size_bound(result.extract.selected_live_atoms));
  const std::size_t initial_capacity =
      std::min(emitted_bound, retained_bound);
  HdtComponentGraph intermediate;
  intermediate.palette = parent.components.palette;
  intermediate.source = parent.components.source;
  intermediate.owner_count = parent.components.owner_count;
  intermediate.background_restricted =
      parent.components.background_restricted;
  intermediate.units.reserve(initial_capacity);
  const std::size_t max_dimension = static_cast<std::size_t>(
      full_bdp_max_dimension(parent.components.palette));
  const std::size_t payload_capacity = checked_hdt_size_mul(
      initial_capacity, max_dimension,
      "Extract payload capacity overflow");
  intermediate.current_payload.reserve(payload_capacity);
  intermediate.background_payload.reserve(payload_capacity);
  std::vector<EmittedPartitionUnit> emitted;
  std::vector<HdtPortKey> emitted_ports;
  emitted.reserve(emitted_bound);
  const std::size_t emitted_port_bound = checked_hdt_size_mul(
      emitted_bound, 2u, "Extract port capacity overflow");
  emitted_ports.reserve(emitted_port_bound);
  int extracted_root_unit = -1;
  int extracted_source_root_unit = -1;
  auto emit_atom = [&](const HdtNode& node) {
    const HdtExpandedAtom& atom =
        hdt.atoms[static_cast<std::size_t>(node.atom_descriptor)];
    const bool target = selector.contains_owner(atom.owner);
    const bool is_parent_root_atom =
        atom.component_unit == parent.components.root_unit &&
        atom.component_kind == parent.components.units[static_cast<std::size_t>(
                                   parent.components.root_unit)]
                                   .kind;
    const bool is_parent_root_anchor =
        is_parent_root_atom &&
        (atom.component_kind == HdtComponentKind::IAnchor ||
         (atom.component_kind == HdtComponentKind::LiveLeaf && target));
    const int emitted_unit = checked_hdt_int_size(
        intermediate.units.size(),
        "Extract exceeds integer node domain");
    FullBdpDual materialized;
    const FullBdpDual* base_dp = nullptr;
    if (!(atom.component_kind == HdtComponentKind::LiveLeaf && target) &&
        atom.component_kind != HdtComponentKind::IAnchor &&
        atom.component_unit >= 0 &&
        atom.component_unit < parent.components.component_units()) {
      const HdtComponentUnit& source = parent.components.units[
          static_cast<std::size_t>(atom.component_unit)];
      if (source.kind == atom.component_kind) {
        materialized =
            hdt_component_dp(parent.components, atom.component_unit);
        base_dp = &materialized;
      }
    }
    if (!(atom.component_kind == HdtComponentKind::LiveLeaf && target) &&
        atom.component_kind != HdtComponentKind::IAnchor &&
        base_dp == nullptr) {
      materialized = hdt_node_dp(
          hdt, static_cast<int>(&node - hdt.nodes.data()));
      base_dp = &materialized;
    }
    HdtComponentUnit unit;
    unit.kind = atom.component_kind;
    unit.stable_id = static_cast<std::uint64_t>(node.atom_descriptor);
    unit.source_vertex_id = atom.source_vertex_id;
    unit.weight = atom.weight;
    unit.background_color = atom.background_color;
    unit.current_color = target ? atom.current_color : atom.background_color;
    unit.owner = target ? atom.owner : -1;
    if (atom.component_kind == HdtComponentKind::LiveLeaf && !target) {
      unit.kind = HdtComponentKind::C1Opaque;
    } else if (unit.kind == HdtComponentKind::LiveLeaf) {
      unit.dp_dimension = static_cast<std::uint8_t>(full_bdp_dimension(
          intermediate.palette, FullBdpInterface::C1));
      unit.dp_orientation = FullBdpOrientation::U;
      unit.current_color_mask = unit.weight == 0.0
                                    ? 0
                                    : static_cast<std::uint8_t>(
                                          1u << unit.current_color);
      unit.background_color_mask = unit.weight == 0.0
                                       ? 0
                                       : static_cast<std::uint8_t>(
                                             1u << unit.background_color);
      unit.current_equals_background =
          unit.current_color == unit.background_color || unit.weight == 0.0;
    }
    EmittedPartitionUnit out;
    out.port_begin = checked_hdt_int_size(
        emitted_ports.size(), "Extract exceeds its port index domain");
    for (int i = 0; i < atom.port_count; ++i) {
      emitted_ports.push_back(hdt.atom_ports[static_cast<std::size_t>(
          atom.port_begin + i)]);
    }
    out.port_count = checked_hdt_int_size(
                         emitted_ports.size(),
                         "Extract exceeds its port index domain") -
                     out.port_begin;
    intermediate.units.push_back(std::move(unit));
    if (base_dp != nullptr) {
      if (target) {
        hdt_component_store_dp(
            intermediate, intermediate.units.back(), *base_dp);
      } else {
        const FullBdpDual frozen = frozen_dual(base_dp->background);
        hdt_component_store_dp(
            intermediate, intermediate.units.back(), frozen);
      }
    }
    if (is_parent_root_atom) extracted_source_root_unit = emitted_unit;
    if (is_parent_root_anchor) extracted_root_unit = emitted_unit;
    emitted.push_back(out);
  };
  auto emit_opaque = [&](int node_index) {
    const HdtNode& node =
        hdt.nodes[static_cast<std::size_t>(node_index)];
    HdtComponentUnit unit;
    unit.kind = opaque_kind(node.output_kind);
    unit.stable_id =
        (std::uint64_t{1} << 63u) | static_cast<std::uint64_t>(node_index);
    unit.source_vertex_id = node.source_center_key;
    unit.owner = -1;
    EmittedPartitionUnit out;
    out.port_begin = checked_hdt_int_size(
        emitted_ports.size(), "Extract exceeds its port index domain");
    for (int i = 0; i < node.output_port_count; ++i) {
      emitted_ports.push_back(
          node.output_ports[static_cast<std::size_t>(i)]);
    }
    out.port_count = checked_hdt_int_size(
                         emitted_ports.size(),
                         "Extract exceeds its port index domain") -
                     out.port_begin;
    intermediate.units.push_back(std::move(unit));
    const FullBdpDual frozen =
        frozen_dual(hdt_node_background(hdt, node_index));
    hdt_component_store_dp(
        intermediate, intermediate.units.back(), frozen);
    emitted.push_back(out);
    ++result.extract.hanging_components;
  };

  std::vector<int> stack{hdt.root};
  while (!stack.empty()) {
    const int node_index = stack.back();
    stack.pop_back();
    const HdtNode& node =
        hdt.nodes[static_cast<std::size_t>(node_index)];
    if (node.operation == HdtOperation::Atom) {
      emit_atom(node);
      continue;
    }
    if (workspace.hdt_mark_epoch[static_cast<std::size_t>(node_index)] !=
        epoch) {
      emit_opaque(node_index);
      continue;
    }
    // Push right first so the stable left operand is emitted first.
    stack.push_back(node.right);
    stack.push_back(node.left);
  }

  const int edge_slots = intermediate.source
                             ? intermediate.source->edge_slots
                             : checked_hdt_int_mul_add(
                                   intermediate.component_units(), 2, 1,
                                   "Extract edge workspace overflow");
  const int vertex_slots = intermediate.source
                               ? checked_hdt_int_size(
                                     intermediate.source->parent.size(),
                                     "Extract source exceeds vertex domain")
                               : checked_hdt_int_add(
                                     intermediate.component_units(), 1,
                                     "Extract vertex workspace overflow");
  if (edge_slots < 0 || vertex_slots < 0) {
    throw std::logic_error("Extract has negative workspace domain");
  }
  workspace.ordinary_epoch.resize(static_cast<std::size_t>(edge_slots), 0);
  workspace.ordinary_unit.resize(static_cast<std::size_t>(edge_slots), -1);
  workspace.ordinary_port.resize(static_cast<std::size_t>(edge_slots));
  workspace.center_epoch.resize(static_cast<std::size_t>(vertex_slots), 0);
  workspace.center_anchor.resize(static_cast<std::size_t>(vertex_slots), -1);
  std::vector<int> touched_ordinary;
  touched_ordinary.reserve(checked_hdt_size_mul(
      emitted.size(), 2u, "Extract touched-edge capacity overflow"));
  const int emitted_count = checked_hdt_int_size(
      emitted.size(), "Extract exceeds integer node domain");
  for (int unit_index = 0; unit_index < emitted_count; ++unit_index) {
    const HdtComponentUnit& unit =
        intermediate.units[static_cast<std::size_t>(unit_index)];
    if (unit.kind == HdtComponentKind::IAnchor &&
        unit.source_vertex_id >= 0 &&
        unit.source_vertex_id < vertex_slots) {
      const std::size_t source =
          static_cast<std::size_t>(unit.source_vertex_id);
      if (workspace.center_epoch[source] == epoch) {
        throw std::logic_error("duplicate extracted I anchor");
      }
      workspace.center_epoch[source] = epoch;
      workspace.center_anchor[source] = unit_index;
    }
  }

  for (int unit_index = 0; unit_index < emitted_count; ++unit_index) {
    const EmittedPartitionUnit& item =
        emitted[static_cast<std::size_t>(unit_index)];
    const HdtComponentUnit& unit =
        intermediate.units[static_cast<std::size_t>(unit_index)];
    for (int port_index = 0; port_index < item.port_count; ++port_index) {
      const HdtPortKey& port = emitted_ports[static_cast<std::size_t>(
          item.port_begin + port_index)];
      if (port.kind == HdtPortKind::Center) {
        if (unit.kind == HdtComponentKind::IAnchor) continue;
        if (port.source_id < 0 || port.source_id >= vertex_slots ||
            workspace.center_epoch[static_cast<std::size_t>(port.source_id)] !=
                epoch ||
            workspace.center_anchor[static_cast<std::size_t>(port.source_id)] <
                0) {
          throw std::logic_error("extracted G lacks its I center anchor");
        }
        const int anchor = workspace.center_anchor[
            static_cast<std::size_t>(port.source_id)];
        const HdtPortKey center{HdtPortKind::Center, port.source_id,
                                      HdtPortOrientation::Center};
        connect_output_edge(intermediate, anchor, unit_index, center, center);
        continue;
      }
      if (port.source_id < 0 || port.source_id >= edge_slots) {
        throw std::logic_error("extracted ordinary source key out of range");
      }
      const std::size_t source = static_cast<std::size_t>(port.source_id);
      if (workspace.ordinary_epoch[source] != epoch) {
        workspace.ordinary_epoch[source] = epoch;
        workspace.ordinary_unit[source] = unit_index;
        workspace.ordinary_port[source] = port;
        touched_ordinary.push_back(port.source_id);
      } else if (workspace.ordinary_unit[source] >= 0) {
        const int first = workspace.ordinary_unit[source];
        connect_output_edge(intermediate, first, unit_index,
                            workspace.ordinary_port[source], port);
        workspace.ordinary_unit[source] = -2;
      } else {
        throw std::logic_error("ordinary key appears more than twice in Extract");
      }
    }
  }
  for (int source_id : touched_ordinary) {
    const int unmatched = workspace.ordinary_unit[
        static_cast<std::size_t>(source_id)];
    if (unmatched >= 0 &&
        intermediate.units[static_cast<std::size_t>(unmatched)].kind ==
            HdtComponentKind::IAnchor) {
      // C_TO_G consumes the original I--C edge. When the transformed G side is
      // emitted as one opaque component it exposes only its center key, while
      // the structural I atom still remembers the pre-transform ordinary key.
      // A matching emitted C/I endpoint reconnects above; an unmatched I-only
      // occurrence is therefore exactly a consumed transform edge.
      continue;
    }
    if (unmatched >= 0) {
      throw std::logic_error(
          "unmatched ordinary key after Extract: source=" +
          std::to_string(source_id) + " unit=" +
          std::to_string(unmatched));
    }
  }

  intermediate.root_unit = extracted_root_unit >= 0
                               ? extracted_root_unit
                               : extracted_source_root_unit;
  if (intermediate.root_unit < 0 && intermediate.source) {
    const int source_root = intermediate.source->root_vertex;
    if (source_root >= 0 && source_root < vertex_slots) {
      if (workspace.center_epoch[static_cast<std::size_t>(source_root)] ==
          epoch) {
        intermediate.root_unit = workspace.center_anchor[
            static_cast<std::size_t>(source_root)];
      }
    }
  }
  if (intermediate.root_unit < 0 && !intermediate.units.empty()) {
    for (int unit = 0; unit < intermediate.component_units(); ++unit) {
      if (intermediate.units[static_cast<std::size_t>(unit)].kind ==
          HdtComponentKind::IAnchor) {
        intermediate.root_unit = unit;
        break;
      }
    }
  }
  if (intermediate.root_unit < 0 && !intermediate.units.empty()) {
    for (int unit = 0; unit < intermediate.component_units(); ++unit) {
      const HdtComponentUnit& candidate =
          intermediate.units[static_cast<std::size_t>(unit)];
      if (candidate.kind == HdtComponentKind::LiveLeaf &&
          selector.contains_owner(candidate.owner)) {
        intermediate.root_unit = unit;
        break;
      }
    }
  }
  if (intermediate.root_unit < 0 && !intermediate.units.empty()) {
    if (result.extract.selected_live_atoms != 0) {
      throw std::logic_error("Extract emitted no selected root anchor");
    }
    intermediate.root_unit = 0;
  }
  rebuild_component_graph_indices(intermediate);
#ifndef NDEBUG
  std::string error;
  if (!validate_component_graph(intermediate, &error)) {
    throw std::logic_error("Extract intermediate invalid: " + error);
  }
#endif
  result.extract.intermediate_units = intermediate.component_units();
  std::size_t workspace_bytes = 0;
  checked_hdt_accumulate_bytes(
      workspace_bytes, workspace.hdt_mark_epoch.capacity(),
      sizeof(std::uint64_t), "Extract workspace byte count exceeds size_t");
  checked_hdt_accumulate_bytes(
      workspace_bytes, workspace.ordinary_epoch.capacity(),
      sizeof(std::uint64_t), "Extract workspace byte count exceeds size_t");
  checked_hdt_accumulate_bytes(
      workspace_bytes, workspace.ordinary_unit.capacity(), sizeof(int),
      "Extract workspace byte count exceeds size_t");
  checked_hdt_accumulate_bytes(
      workspace_bytes, workspace.ordinary_port.capacity(),
      sizeof(HdtPortKey),
      "Extract workspace byte count exceeds size_t");
  checked_hdt_accumulate_bytes(
      workspace_bytes, workspace.center_epoch.capacity(),
      sizeof(std::uint64_t), "Extract workspace byte count exceeds size_t");
  checked_hdt_accumulate_bytes(
      workspace_bytes, workspace.center_anchor.capacity(), sizeof(int),
      "Extract workspace byte count exceeds size_t");
  workspace.peak_bytes = std::max(workspace.peak_bytes, workspace_bytes);
  const int identity_bound =
      hdt_contract_size_bound(result.extract.selected_live_atoms);
  const bool root_is_anchor =
      intermediate.root_unit >= 0 &&
      is_selected_contract_anchor(
          intermediate.units[static_cast<std::size_t>(intermediate.root_unit)],
          selector);
  if (intermediate.component_units() <= identity_bound && root_is_anchor) {
    const long long scalar_fmas_before = full_bdp_thread_scalar_fmas();
    result.contract.mode = HdtContractMode::LightCopy;
    result.contract.input_units = intermediate.component_units();
    result.contract.output_units = intermediate.component_units();
    result.contract.selected_live_atoms =
        result.extract.selected_live_atoms;
    result.contract.size_bound = identity_bound;
    result.contract.size_bound_ok = true;
    const auto rebuild_started = std::chrono::steady_clock::now();
    result.state = build_hdt_state(std::move(intermediate),
                                         owner_index_view(selector));
    result.contract.hdt_rebuild_ns =
        std::chrono::duration_cast<std::chrono::nanoseconds>(
            std::chrono::steady_clock::now() - rebuild_started)
            .count();
    result.contract.scalar_fmas =
        full_bdp_thread_scalar_fmas() - scalar_fmas_before;
    return result;
  }
  HdtState intermediate_state;
  intermediate_state.components = std::move(intermediate);
  HdtContractResult contracted = hdt_contract(
      intermediate_state, selector, HdtContractMode::LightCopy);
  result.contract = contracted.stats;
  result.state = std::move(contracted.state);
  return result;
}

HdtExtractContractResult hdt_extract_contract(
    const HdtState& parent,
    const HdtLiveSelector& selector) {
  HdtExtractWorkspace workspace;
  return hdt_extract_contract(parent, selector, workspace);
}

}  // namespace wqfm::internal
