#include "hat/hdt/hdt_state.hpp"

#include <algorithm>
#include <bit>
#include <cmath>
#include <deque>
#include <stdexcept>
#include <type_traits>
#include <utility>

namespace wqfm::internal {
namespace {

enum class ActiveKind : std::uint8_t {
  I,
  C,
  G,
};

struct ActivePortList {
  std::array<HdtPortKey, 2> values{};
  std::uint8_t count = 0;
};

struct ActiveComponent {
  ActiveKind kind = ActiveKind::I;
  int hdt_root = -1;
  int source_center = -1;
  std::uint64_t stable_order = 0;
  ActivePortList external_ports;
};

void append_external_port(ActivePortList& ports,
                          const HdtPortKey& port,
                          bool unique = false) {
  if (unique) {
    for (int i = 0; i < ports.count; ++i) {
      if (ports.values[static_cast<std::size_t>(i)] == port) return;
    }
  }
  if (ports.count >= ports.values.size()) {
    throw std::logic_error("active component exceeds two external ports");
  }
  ports.values[ports.count++] = port;
}

struct ActiveEdge {
  int a = -1;
  int b = -1;
  HdtPortKey a_port;
  HdtPortKey b_port;
};

struct ActiveGraph {
  std::vector<ActiveComponent> components;
  std::vector<ActiveEdge> edges;
  int root = -1;
};

struct ActiveAdjacency {
  std::vector<int> offsets;
  std::vector<int> edge_indices;
};

struct MutableActiveLink {
  int child = -1;
  int previous = -1;
  int next = -1;
  HdtPortKey parent_port;
  HdtPortKey child_port;
};

struct MutableActiveNode {
  int hdt_root = -1;
  int source_center = -1;
  int parent = -1;
  int parent_link = -1;
  int first_child = -1;
  int last_child = -1;
  ActiveKind kind = ActiveKind::I;
};

struct MutableActiveTree {
  std::vector<MutableActiveNode> nodes;
  std::vector<ActivePortList> external_ports;
  std::vector<MutableActiveLink> links;
  int root = -1;
  int active_nodes = 0;
};

enum class PlannedOperation : std::uint8_t {
  Survivor,
  GG,
  IG,
  CC,
};

struct PlannedGroup {
  PlannedOperation operation = PlannedOperation::Survivor;
  int a = -1;
  int b = -1;
};

struct RoundWorkspace {
  std::vector<int> adjacency_cursor;
  std::vector<int> parent;
  std::vector<int> parent_edge;
  std::vector<int> preorder;
  std::vector<int> traversal_stack;
  std::vector<std::uint8_t> used;
  std::vector<int> mutable_used_epoch;
  std::vector<PlannedGroup> compositions;
  std::vector<PlannedGroup> ig_compositions;
  std::vector<PlannedGroup> groups;
  std::vector<int> group_scratch;
  std::vector<int> group_of;
  std::vector<int> center_edge_seen;
};

template <typename T>
void add_vector_capacity_bytes(std::size_t& bytes,
                               const std::vector<T>& values) {
  checked_hdt_accumulate_bytes(
      bytes, values.capacity(), sizeof(T),
      "HDT builder scratch byte count exceeds size_t");
}

std::size_t active_graph_capacity_bytes(const ActiveGraph& active) {
  std::size_t bytes = 0;
  add_vector_capacity_bytes(bytes, active.components);
  add_vector_capacity_bytes(bytes, active.edges);
  return bytes;
}

std::size_t live_builder_scratch_bytes(const ActiveGraph& active,
                                       const ActiveAdjacency& adjacency,
                                       const MutableActiveTree& tree,
                                       const RoundWorkspace& workspace) {
  std::size_t bytes = active_graph_capacity_bytes(active);
  add_vector_capacity_bytes(bytes, adjacency.offsets);
  add_vector_capacity_bytes(bytes, adjacency.edge_indices);
  add_vector_capacity_bytes(bytes, tree.nodes);
  add_vector_capacity_bytes(bytes, tree.external_ports);
  add_vector_capacity_bytes(bytes, tree.links);
  add_vector_capacity_bytes(bytes, workspace.adjacency_cursor);
  add_vector_capacity_bytes(bytes, workspace.parent);
  add_vector_capacity_bytes(bytes, workspace.parent_edge);
  add_vector_capacity_bytes(bytes, workspace.preorder);
  add_vector_capacity_bytes(bytes, workspace.traversal_stack);
  add_vector_capacity_bytes(bytes, workspace.used);
  add_vector_capacity_bytes(bytes, workspace.mutable_used_epoch);
  add_vector_capacity_bytes(bytes, workspace.compositions);
  add_vector_capacity_bytes(bytes, workspace.ig_compositions);
  add_vector_capacity_bytes(bytes, workspace.groups);
  add_vector_capacity_bytes(bytes, workspace.group_scratch);
  add_vector_capacity_bytes(bytes, workspace.group_of);
  add_vector_capacity_bytes(bytes, workspace.center_edge_seen);
  return bytes;
}

void build_active_adjacency(const ActiveGraph& graph,
                            ActiveAdjacency& adjacency,
                            std::vector<int>& cursor) {
  (void)checked_hdt_int_size(
      graph.components.size(), "active graph exceeds node domain");
  const int edge_count = checked_hdt_int_size(
      graph.edges.size(), "active graph exceeds edge domain");
  adjacency.offsets.assign(
      checked_hdt_size_add(graph.components.size(), 1u,
                           "active offsets overflow"),
      0);
  for (const ActiveEdge& edge : graph.edges) {
    ++adjacency.offsets[static_cast<std::size_t>(edge.a) + 1u];
    ++adjacency.offsets[static_cast<std::size_t>(edge.b) + 1u];
  }
  for (std::size_t i = 1; i < adjacency.offsets.size(); ++i) {
    adjacency.offsets[i] += adjacency.offsets[i - 1u];
  }
  const std::size_t incidence_count = checked_hdt_size_mul(
      graph.edges.size(), 2u, "active incidence storage overflow");
  (void)checked_hdt_int_size(
      incidence_count, "active incidence exceeds integer domain");
  adjacency.edge_indices.assign(incidence_count, -1);
  cursor = adjacency.offsets;
  for (int edge_index = 0; edge_index < edge_count; ++edge_index) {
    const ActiveEdge& edge = graph.edges[static_cast<std::size_t>(edge_index)];
    adjacency.edge_indices[static_cast<std::size_t>(
        cursor[static_cast<std::size_t>(edge.a)]++)] = edge_index;
    adjacency.edge_indices[static_cast<std::size_t>(
        cursor[static_cast<std::size_t>(edge.b)]++)] = edge_index;
  }
}

int other_endpoint(const ActiveEdge& edge, int component) {
  return edge.a == component ? edge.b : edge.a;
}

HdtPortKey endpoint_port(const ActiveEdge& edge, int component) {
  return edge.a == component ? edge.a_port : edge.b_port;
}

void rooted_active_tree(const ActiveGraph& graph,
                        const ActiveAdjacency& adjacency,
                        std::vector<int>& parent,
                        std::vector<int>& parent_edge,
                        std::vector<int>& preorder,
                        std::vector<int>& stack) {
  const int n = checked_hdt_int_size(
      graph.components.size(), "active graph exceeds node domain");
  parent.assign(static_cast<std::size_t>(n), -2);
  parent_edge.assign(static_cast<std::size_t>(n), -1);
  preorder.clear();
  preorder.reserve(static_cast<std::size_t>(n));
  if (graph.root < 0 || graph.root >= n) return;
  parent[static_cast<std::size_t>(graph.root)] = -1;
  stack.clear();
  stack.reserve(static_cast<std::size_t>(n));
  stack.push_back(graph.root);
  while (!stack.empty()) {
    const int component = stack.back();
    stack.pop_back();
    preorder.push_back(component);
    const int begin = adjacency.offsets[static_cast<std::size_t>(component)];
    const int end = adjacency.offsets[static_cast<std::size_t>(component) + 1u];
    for (int p = end - 1; p >= begin; --p) {
      const int edge_index =
          adjacency.edge_indices[static_cast<std::size_t>(p)];
      const int next = other_endpoint(
          graph.edges[static_cast<std::size_t>(edge_index)], component);
      if (parent[static_cast<std::size_t>(next)] != -2) continue;
      parent[static_cast<std::size_t>(next)] = component;
      parent_edge[static_cast<std::size_t>(next)] = edge_index;
      stack.push_back(next);
    }
  }
  if (static_cast<int>(preorder.size()) != n) {
    throw std::logic_error("HDT active graph disconnected");
  }
}

HdtOutputKind output_kind_of(const FullBdpVec& value) {
  switch (value.interface) {
    case FullBdpInterface::Closed: return HdtOutputKind::Closed;
    case FullBdpInterface::C1:
      return value.orientation == FullBdpOrientation::U
                 ? HdtOutputKind::C1U
                 : HdtOutputKind::C1D;
    case FullBdpInterface::C2: return HdtOutputKind::C2;
    case FullBdpInterface::G: return HdtOutputKind::G;
  }
  throw std::logic_error("unknown boundary-DP output interface");
}

HdtOutputKind output_kind_of(const HdtComponentUnit& unit) {
  switch (unit.kind) {
    case HdtComponentKind::IAnchor:
      return HdtOutputKind::I;
    case HdtComponentKind::ClosedOpaque:
      return HdtOutputKind::Closed;
    case HdtComponentKind::LiveLeaf:
      return HdtOutputKind::C1U;
    case HdtComponentKind::C1Opaque:
      return unit.dp_orientation == FullBdpOrientation::D
                 ? HdtOutputKind::C1D
                 : HdtOutputKind::C1U;
    case HdtComponentKind::C2Opaque:
      return HdtOutputKind::C2;
    case HdtComponentKind::GOpaque:
      return HdtOutputKind::G;
  }
  throw std::logic_error("unknown component output kind");
}

std::pair<FullBdpInterface, FullBdpOrientation> descriptor_for_kind(
    HdtOutputKind kind) {
  switch (kind) {
    case HdtOutputKind::Closed:
      return {FullBdpInterface::Closed, FullBdpOrientation::None};
    case HdtOutputKind::C1U:
      return {FullBdpInterface::C1, FullBdpOrientation::U};
    case HdtOutputKind::C1D:
      return {FullBdpInterface::C1, FullBdpOrientation::D};
    case HdtOutputKind::C2:
      return {FullBdpInterface::C2, FullBdpOrientation::None};
    case HdtOutputKind::G:
      return {FullBdpInterface::G, FullBdpOrientation::None};
    case HdtOutputKind::I: break;
  }
  throw std::logic_error("I node has no boundary-DP payload");
}

int packed_payload_width(int palette,
                         FullBdpInterface interface,
                         std::uint8_t color_mask) {
  return std::max(
      1, full_bdp_active_coordinate_count(palette, interface, color_mask));
}

std::uint8_t packed_payload_class_index(int width) {
  switch (width) {
    case 1: return 0;
    case 2: return 1;
    case 3: return 2;
    case 4: return 3;
    case 5: return 4;
    case 6: return 5;
    case 7: return 6;
    case 8: return 7;
    case 13: return 8;
    case 15: return 9;
    case 20: return 10;
    case 22: return 11;
    case 23: return 12;
    case 24: return 13;
    case 32: return 14;
    case 52: return 15;
    case 56: return 16;
  }
  throw std::logic_error("unreachable packed-payload width");
}

constexpr std::array<int, 17> kPackedPayloadWidths{
    1, 2, 3, 4, 5, 6, 7, 8, 13, 15, 20, 22, 23, 24, 32, 52, 56};

int packed_payload_width(std::uint8_t payload_class) {
  if (payload_class >= kPackedPayloadWidths.size()) {
    throw std::logic_error("invalid packed-payload class");
  }
  return kPackedPayloadWidths[payload_class];
}

int allocate_packed_payload(HdtPackedPayloadArena& arena,
                            std::uint8_t class_index) {
  const int width = packed_payload_width(class_index);
  HdtPackedPayloadClass& storage =
      arena.classes[static_cast<std::size_t>(class_index)];
  if (storage.free_head >= 0) {
    const int offset = storage.free_head;
    storage.free_head =
        static_cast<int>(storage.values[static_cast<std::size_t>(offset)]);
    return offset;
  }
  const int offset = checked_hdt_int_size(
      storage.values.size(),
      "packed payload exceeds its integer offset domain");
  const std::size_t next_size = checked_hdt_size_add(
      storage.values.size(), static_cast<std::size_t>(width),
      "packed payload capacity overflow");
  (void)checked_hdt_int_size(
      next_size, "packed payload exceeds its integer offset domain");
  storage.values.resize(next_size);
  return offset;
}

void release_packed_payload(HdtPackedPayloadArena& arena,
                            std::uint8_t class_index,
                            int offset) {
  if (offset < 0) return;
  const int width = packed_payload_width(class_index);
  HdtPackedPayloadClass& storage =
      arena.classes[static_cast<std::size_t>(class_index)];
  if (static_cast<std::size_t>(offset) + width >
      storage.values.size()) {
    throw std::logic_error("released payload outside arena");
  }
  storage.values[static_cast<std::size_t>(offset)] =
      static_cast<double>(storage.free_head);
  storage.free_head = offset;
}

void pack_full_value(int palette,
                     const FullBdpVec& value,
                     std::uint8_t color_mask,
                     double* output) {
  std::uint64_t active = full_bdp_active_coordinate_bits(
      palette, value.interface, color_mask);
  int packed = 0;
  while (active != 0) {
    const int coordinate = std::countr_zero(active);
    output[packed++] = value.values[static_cast<std::size_t>(coordinate)];
    active &= active - 1;
  }
}

const double* node_current_data(const Hdt& hdt,
                                const HdtNode& node) noexcept {
  const HdtPackedPayloadClass& storage =
      hdt.current_payload.classes[node.current_payload_class];
  return storage.values.data() + node.dp_offset;
}

double* node_current_data(Hdt& hdt,
                          const HdtNode& node) noexcept {
  HdtPackedPayloadClass& storage =
      hdt.current_payload.classes[node.current_payload_class];
  return storage.values.data() + node.dp_offset;
}

const double* node_background_data(const Hdt& hdt,
                                   const HdtNode& node) noexcept {
  const bool aliases_current =
      hdt.background_alias_current || node.background_dp_offset < 0;
  const HdtPackedPayloadArena& arena =
      aliases_current ? hdt.current_payload : hdt.background_payload;
  const std::uint8_t payload_class =
      aliases_current ? node.current_payload_class
                      : node.background_payload_class;
  const int offset =
      aliases_current ? node.dp_offset : node.background_dp_offset;
  return arena.classes[payload_class].values.data() + offset;
}

double* node_background_data(Hdt& hdt,
                             const HdtNode& node) noexcept {
  const bool aliases_current =
      hdt.background_alias_current || node.background_dp_offset < 0;
  HdtPackedPayloadArena& arena =
      aliases_current ? hdt.current_payload : hdt.background_payload;
  const std::uint8_t payload_class =
      aliases_current ? node.current_payload_class
                      : node.background_payload_class;
  const int offset =
      aliases_current ? node.dp_offset : node.background_dp_offset;
  return arena.classes[payload_class].values.data() + offset;
}

FullBdpConstView node_current_view(const Hdt& hdt,
                                   const HdtNode& node) {
  const auto [interface, orientation] = descriptor_for_kind(node.output_kind);
  return FullBdpConstView{
      interface, orientation, node.dp_dimension,
      node_current_data(hdt, node)};
}

FullBdpConstView node_background_view(const Hdt& hdt,
                                      const HdtNode& node) {
  const auto [interface, orientation] = descriptor_for_kind(node.output_kind);
  return FullBdpConstView{
      interface, orientation, node.dp_dimension,
      node_background_data(hdt, node)};
}

void allocate_node_background(Hdt& hdt, HdtNode& node) {
  if (hdt.background_alias_current) {
    throw std::logic_error("cannot allocate a frozen HDT background");
  }
  if (node.background_dp_offset >= 0 || node.dp_dimension == 0) return;
  const FullBdpInterface interface =
      descriptor_for_kind(node.output_kind).first;
  node.background_payload_class = packed_payload_class_index(
      packed_payload_width(hdt.palette, interface,
                           node.background_color_mask));
  node.background_dp_offset =
      allocate_packed_payload(hdt.background_payload,
                              node.background_payload_class);
}

void promote_node_background(Hdt& hdt, HdtNode& node) {
  if (node.background_dp_offset >= 0 || node.dp_dimension == 0) return;
  allocate_node_background(hdt, node);
  const FullBdpInterface interface =
      descriptor_for_kind(node.output_kind).first;
  const int active = full_bdp_active_coordinate_count(
      hdt.palette, interface, node.current_color_mask);
  std::copy_n(node_current_data(hdt, node), active,
              node_background_data(hdt, node));
}

void allocate_node_payload(Hdt& hdt,
                           HdtNode& node,
                           HdtOutputKind output_kind) {
  node.output_kind = output_kind;
  if (output_kind == HdtOutputKind::I) {
    node.dp_offset = -1;
    node.background_dp_offset = -1;
    node.current_payload_class = 0xff;
    node.background_payload_class = 0xff;
    node.dp_dimension = 0;
    node.current_equals_background = true;
    return;
  }
  const auto [interface, orientation] = descriptor_for_kind(output_kind);
  (void)orientation;
  const int dimension = full_bdp_dimension(hdt.palette, interface);
  node.current_payload_class = packed_payload_class_index(
      packed_payload_width(hdt.palette, interface, node.current_color_mask));
  node.dp_offset =
      allocate_packed_payload(hdt.current_payload,
                              node.current_payload_class);
  node.background_dp_offset = -1;
  node.background_payload_class = 0xff;
  node.dp_dimension = static_cast<std::uint8_t>(dimension);
}

void store_node_payload(Hdt& hdt,
                        HdtNode& node,
                        const FullBdpDual& value,
                        HdtOutputKind output_kind) {
  node.current_color_mask =
      full_bdp_value_color_mask(hdt.palette, value.current);
  node.background_color_mask = value.current_equals_background
                                   ? node.current_color_mask
                                   : full_bdp_value_color_mask(
                                         hdt.palette, value.background);
  node.current_equals_background = value.current_equals_background;
  allocate_node_payload(hdt, node, output_kind);
  if (output_kind == HdtOutputKind::I) return;
  pack_full_value(hdt.palette, value.current, node.current_color_mask,
                  node_current_data(hdt, node));
  if (hdt.background_alias_current) {
    if (!value.current_equals_background) {
      throw std::logic_error("non-frozen atom in frozen HDT fold");
    }
  } else if (!value.current_equals_background) {
    allocate_node_background(hdt, node);
    pack_full_value(hdt.palette, value.background, node.background_color_mask,
                    node_background_data(hdt, node));
  }
}

template <class Eval>
void evaluate_unary_payload_color_mask(Hdt& hdt,
                                       HdtNode& output,
                                       HdtOutputKind output_kind,
                                       int input_index,
                                       Eval&& eval) {
  const HdtNode& input =
      hdt.nodes[static_cast<std::size_t>(input_index)];
  output.current_color_mask = input.current_color_mask;
  output.background_color_mask = input.background_color_mask;
  allocate_node_payload(hdt, output, output_kind);
  eval(output.current_color_mask,
       node_current_data(hdt, input), node_current_data(hdt, output));
  output.current_equals_background = input.current_equals_background;
  if (!input.current_equals_background) {
    if (hdt.background_alias_current) {
      throw std::logic_error("non-frozen unary node in frozen HDT fold");
    }
    allocate_node_background(hdt, output);
    eval(output.background_color_mask,
         node_background_data(hdt, input),
         node_background_data(hdt, output));
  }
}

template <class Eval>
void evaluate_binary_payload_color_mask(Hdt& hdt,
                                        HdtNode& output,
                                        HdtOutputKind output_kind,
                                        int left_index,
                                        int right_index,
                                        Eval&& eval) {
  const HdtNode& left =
      hdt.nodes[static_cast<std::size_t>(left_index)];
  const HdtNode& right =
      hdt.nodes[static_cast<std::size_t>(right_index)];
  output.current_color_mask =
      left.current_color_mask | right.current_color_mask;
  output.background_color_mask =
      left.background_color_mask | right.background_color_mask;
  allocate_node_payload(hdt, output, output_kind);
  eval(left.current_color_mask, right.current_color_mask,
       node_current_data(hdt, left), node_current_data(hdt, right),
       node_current_data(hdt, output));
  output.current_equals_background =
      left.current_equals_background && right.current_equals_background;
  if (!output.current_equals_background) {
    if (hdt.background_alias_current) {
      throw std::logic_error("non-frozen binary node in frozen HDT fold");
    }
    allocate_node_background(hdt, output);
    const double* left_background = node_background_data(hdt, left);
    const double* right_background = node_background_data(hdt, right);
    eval(left.background_color_mask, right.background_color_mask,
         left_background, right_background,
         node_background_data(hdt, output));
  }
}

template <class Eval>
void replace_unary_payload_color_mask(Hdt& hdt,
                                      HdtNode& node,
                                      HdtOutputKind output_kind,
                                      Eval&& eval) {
  const HdtNode input = node;
  const FullBdpInterface output_interface =
      descriptor_for_kind(output_kind).first;
  const int output_dimension =
      full_bdp_dimension(hdt.palette, output_interface);
  if (output_dimension > input.dp_dimension) {
    throw std::logic_error("in-place masked unary transform grew its type");
  }
  std::array<double, kFullBdpMaxDimension> current;
  eval(input.current_color_mask,
       node_current_data(hdt, input), current.data());
  std::array<double, kFullBdpMaxDimension> background;
  if (!input.current_equals_background) {
    eval(input.background_color_mask,
         node_background_data(hdt, input), background.data());
    if (node.background_dp_offset < 0) {
      throw std::logic_error("masked unary node lacks background payload");
    }
  }

  const std::uint8_t next_current_class = packed_payload_class_index(
      packed_payload_width(hdt.palette, output_interface,
                           input.current_color_mask));
  const bool relocate_current =
      input.current_payload_class != next_current_class;
  const std::uint8_t next_background_class =
      input.current_equals_background
          ? std::uint8_t{0xff}
          : packed_payload_class_index(packed_payload_width(
                hdt.palette, output_interface, input.background_color_mask));
  const bool relocate_background =
      !input.current_equals_background &&
      input.background_payload_class != next_background_class;
  if (relocate_current) {
    node.dp_offset =
        allocate_packed_payload(hdt.current_payload, next_current_class);
    node.current_payload_class = next_current_class;
  }
  if (relocate_background) {
    node.background_dp_offset = allocate_packed_payload(
        hdt.background_payload, next_background_class);
    node.background_payload_class = next_background_class;
  }

  node.output_kind = output_kind;
  node.dp_dimension = static_cast<std::uint8_t>(output_dimension);
  node.current_equals_background = input.current_equals_background;
  node.current_color_mask = input.current_color_mask;
  node.background_color_mask = input.background_color_mask;
  const int current_count = full_bdp_active_coordinate_count(
      hdt.palette, output_interface, input.current_color_mask);
  std::copy_n(current.data(), current_count,
              node_current_data(hdt, node));
  if (!input.current_equals_background) {
    const int background_count = full_bdp_active_coordinate_count(
        hdt.palette, output_interface, input.background_color_mask);
    std::copy_n(background.data(), background_count,
                node_background_data(hdt, node));
  }
  if (relocate_current) {
    release_packed_payload(hdt.current_payload,
                           input.current_payload_class, input.dp_offset);
  }
  if (relocate_background) {
    release_packed_payload(hdt.background_payload,
                           input.background_payload_class,
                           input.background_dp_offset);
  }
}

FullBdpVec materialize_view(FullBdpConstView view,
                            int palette,
                            std::uint8_t color_mask) {
  FullBdpVec value;
  value.interface = view.interface;
  value.orientation = view.orientation;
  value.dimension = view.dimension;
  std::fill_n(value.values.data(), view.dimension, 0.0);
  std::uint64_t active = full_bdp_active_coordinate_bits(
      palette, view.interface, color_mask);
  int packed = 0;
  while (active != 0) {
    const int coordinate = std::countr_zero(active);
    value.values[static_cast<std::size_t>(coordinate)] =
        view.values[packed++];
    active &= active - 1;
  }
  return value;
}

FullBdpDual ig_output(int palette,
                      const FullBdpDual& g,
                      const ActivePortList& ports) {
  if (ports.count == 0) return full_bdp_ig_ff(palette, g);
  if (ports.count == 1) {
    FullBdpDual value = full_bdp_ig_tf(palette, g);
    if (ports.values[0].orientation == HdtPortOrientation::Down) {
      value = full_bdp_reverse(palette, value);
    } else if (ports.values[0].orientation !=
               HdtPortOrientation::Up) {
      throw std::logic_error("IG one-port output is not ordinary");
    }
    return value;
  }
  if (ports.count == 2) {
    bool have_up = false;
    bool have_down = false;
    for (const HdtPortKey& port : ports.values) {
      have_up = have_up || port.orientation == HdtPortOrientation::Up;
      have_down = have_down ||
                  port.orientation == HdtPortOrientation::Down;
    }
    if (!have_up || !have_down) {
      throw std::logic_error("IG C2 output lacks U/D ports");
    }
    return full_bdp_ig_tt(palette, g);
  }
  throw std::logic_error("IG output exceeds two ordinary ports");
}

void note_node_shape(Hdt& hdt, const HdtNode& node) {
  const int degree = (node.left >= 0 ? 1 : 0) + (node.right >= 0 ? 1 : 0);
  hdt.stats.max_recompute_degree =
      std::max(hdt.stats.max_recompute_degree, degree);
  hdt.stats.max_update_height =
      std::max(hdt.stats.max_update_height, node.height);
  if (node.atom_count <= 1) return;
  // Since 2*9^7 < 10^7, log_(10/9)(n) < 7*ceil(log2(n)).  This
  // integer certificate is slightly looser than the SODA round constant but
  // still proves the required local O(log n) update height.
  const int ceil_log2 = static_cast<int>(
      std::bit_width(static_cast<unsigned int>(node.atom_count - 1)));
  const int allowed = 7 * ceil_log2 + 1;
  if (node.height > allowed) ++hdt.stats.local_balance_violations;
}

int add_atom(Hdt& hdt,
             HdtExpandedAtom atom,
             const FullBdpDual& dp,
             HdtOutputKind output_kind,
             int source_center) {
  const int atom_index = checked_hdt_int_size(
      hdt.atoms.size(), "HDT exceeds its atom index domain");
  hdt.atoms.push_back(std::move(atom));
  HdtNode node;
  node.atom_descriptor = atom_index;
  node.operation = HdtOperation::Atom;
  node.output_kind = output_kind;
  node.source_center_key = source_center;
  store_node_payload(hdt, node, dp, output_kind);
  node.live_owner = hdt.atoms[static_cast<std::size_t>(atom_index)].owner;
  const int node_index = checked_hdt_int_size(
      hdt.nodes.size(), "HDT exceeds its node index domain");
  note_node_shape(hdt, node);
  hdt.nodes.push_back(std::move(node));
  return node_index;
}

int add_component_atom(Hdt& hdt,
                       HdtExpandedAtom atom,
                       const HdtComponentGraph& graph,
                       int component_index,
                       int source_center) {
  const HdtComponentUnit& component =
      graph.units[static_cast<std::size_t>(component_index)];
  const HdtOutputKind output_kind = output_kind_of(component);
  const int atom_index = checked_hdt_int_size(
      hdt.atoms.size(), "HDT exceeds its atom index domain");
  hdt.atoms.push_back(std::move(atom));
  HdtNode node;
  node.atom_descriptor = atom_index;
  node.operation = HdtOperation::Atom;
  node.source_center_key = source_center;
  node.current_color_mask = component.current_color_mask;
  node.background_color_mask = component.background_color_mask;
  allocate_node_payload(hdt, node, output_kind);
  if (output_kind != HdtOutputKind::I) {
    if (component.dp_dimension != node.dp_dimension) {
      throw std::logic_error("component/HDT atom payload mismatch");
    }
    if (component.kind == HdtComponentKind::LiveLeaf) {
      full_bdp_leaf_color_mask_into(
          hdt.palette, component.current_color, component.weight,
          component.current_color_mask,
          node_current_data(hdt, node));
      node.current_equals_background = component.current_equals_background;
      if (hdt.background_alias_current) {
        if (!component.current_equals_background) {
          throw std::logic_error("non-frozen component in frozen HDT fold");
        }
      } else if (!component.current_equals_background) {
        allocate_node_background(hdt, node);
        full_bdp_leaf_color_mask_into(
            hdt.palette, component.background_color, component.weight,
            component.background_color_mask,
            node_background_data(hdt, node));
      }
    } else {
      if (component.dp_offset < 0) {
        throw std::logic_error("component/HDT atom payload missing");
      }
      const FullBdpInterface interface =
          descriptor_for_kind(output_kind).first;
      FullBdpVec current;
      current.interface = interface;
      current.dimension = node.dp_dimension;
      std::copy_n(graph.current_payload.data() + component.dp_offset,
                  node.dp_dimension, current.values.data());
      pack_full_value(hdt.palette, current, component.current_color_mask,
                      node_current_data(hdt, node));
      if (hdt.background_alias_current) {
        if (!component.current_equals_background) {
          throw std::logic_error("non-frozen component in frozen HDT fold");
        }
      } else if (!component.current_equals_background) {
        node.current_equals_background = false;
        allocate_node_background(hdt, node);
        const double* background =
            graph.background_payload.data() +
            component.background_dp_offset;
        FullBdpVec background_value;
        background_value.interface = interface;
        background_value.dimension = node.dp_dimension;
        std::copy_n(background, node.dp_dimension,
                    background_value.values.data());
        pack_full_value(hdt.palette, background_value,
                        component.background_color_mask,
                        node_background_data(hdt, node));
      }
      node.current_equals_background = component.current_equals_background;
    }
  }
  node.live_owner = hdt.atoms[static_cast<std::size_t>(atom_index)].owner;
  const int node_index = checked_hdt_int_size(
      hdt.nodes.size(), "HDT exceeds its node index domain");
  note_node_shape(hdt, node);
  hdt.nodes.push_back(std::move(node));
  return node_index;
}

void evaluate_transformed_base_color_mask_into(const Hdt& hdt,
                                               const HdtNode& node,
                                               double* output) {
  if (node.operation == HdtOperation::Atom) {
    const HdtExpandedAtom& atom =
        hdt.atoms[static_cast<std::size_t>(node.atom_descriptor)];
    if (atom.component_kind == HdtComponentKind::LiveLeaf) {
      full_bdp_leaf_color_mask_into(
          hdt.palette, atom.current_color, atom.weight,
          atom.current_color_mask, output);
      return;
    }
    throw std::logic_error("immutable transformed atom was marked dirty");
  }
  const HdtNode& left =
      hdt.nodes[static_cast<std::size_t>(node.left)];
  const HdtNode& right =
      hdt.nodes[static_cast<std::size_t>(node.right)];
  switch (node.operation) {
    case HdtOperation::GG:
      throw std::logic_error("G node has an invalid C post-transform");
    case HdtOperation::IG:
      full_bdp_ig_tf_color_mask_into(
          hdt.palette, right.current_color_mask,
          node_current_data(hdt, right), output);
      return;
    case HdtOperation::CC:
      full_bdp_cc_color_mask_into(
          hdt.palette, descriptor_for_kind(left.output_kind).first,
          descriptor_for_kind(right.output_kind).first,
          left.current_color_mask, right.current_color_mask,
          node_current_data(hdt, left), node_current_data(hdt, right), output);
      return;
    case HdtOperation::Atom: break;
  }
  throw std::logic_error("invalid transformed HDT operation");
}

void recompute_node(Hdt& hdt, int node_index) {
  HdtNode& node = hdt.nodes[static_cast<std::size_t>(node_index)];
  bool next_equals_background = false;
  std::uint8_t next_color_mask = 0;
  if (node.operation == HdtOperation::Atom) {
    const HdtExpandedAtom& atom =
        hdt.atoms[static_cast<std::size_t>(node.atom_descriptor)];
    next_equals_background =
        atom.current_color == atom.background_color || atom.weight == 0.0;
    next_color_mask = atom.current_color_mask;
  } else {
    next_equals_background =
        hdt.nodes[static_cast<std::size_t>(node.left)]
            .current_equals_background &&
        hdt.nodes[static_cast<std::size_t>(node.right)]
            .current_equals_background;
    next_color_mask =
        hdt.nodes[static_cast<std::size_t>(node.left)].current_color_mask |
        hdt.nodes[static_cast<std::size_t>(node.right)].current_color_mask;
  }
  if (!next_equals_background && node.background_dp_offset < 0) {
    promote_node_background(hdt, node);
  }
  if (node.output_kind == HdtOutputKind::I) {
    node.current_equals_background = next_equals_background;
    node.current_color_mask = next_color_mask;
    return;
  }

  const FullBdpInterface output_interface =
      descriptor_for_kind(node.output_kind).first;
  const int old_dp_offset = node.dp_offset;
  const std::uint8_t old_payload_class = node.current_payload_class;
  const std::uint8_t next_payload_class = packed_payload_class_index(
      packed_payload_width(hdt.palette, output_interface, next_color_mask));
  const bool relocated = old_payload_class != next_payload_class;
  if (relocated) {
    node.dp_offset =
        allocate_packed_payload(hdt.current_payload, next_payload_class);
    node.current_payload_class = next_payload_class;
  }
  node.current_color_mask = next_color_mask;
  double* const output = node_current_data(hdt, node);

  if (node.post_transform != HdtTransform::None) {
    std::array<double, kFullBdpMaxDimension> base;
    evaluate_transformed_base_color_mask_into(hdt, node, base.data());
    if (node.post_transform == HdtTransform::CToG) {
      full_bdp_c_to_g_color_mask_into(
          hdt.palette, next_color_mask, base.data(), output);
    } else if (node.post_transform == HdtTransform::CloseRoot) {
      full_bdp_close_color_mask_into(
          hdt.palette, next_color_mask, base.data(), output);
    } else {
      throw std::logic_error("unknown HDT transform");
    }
  } else if (node.operation == HdtOperation::Atom) {
    const HdtExpandedAtom& atom =
        hdt.atoms[static_cast<std::size_t>(node.atom_descriptor)];
    if (atom.component_kind == HdtComponentKind::LiveLeaf) {
      full_bdp_leaf_color_mask_into(
          hdt.palette, atom.current_color, atom.weight, next_color_mask,
          output);
    }
  } else {
    const HdtNode& left =
        hdt.nodes[static_cast<std::size_t>(node.left)];
    const HdtNode& right =
        hdt.nodes[static_cast<std::size_t>(node.right)];
    switch (node.operation) {
      case HdtOperation::GG:
        full_bdp_gg_color_mask_into(
            hdt.palette, left.current_color_mask, right.current_color_mask,
            node_current_data(hdt, left), node_current_data(hdt, right),
            output);
        break;
      case HdtOperation::IG:
        if (node.output_kind == HdtOutputKind::Closed) {
          full_bdp_ig_ff_color_mask_into(
              hdt.palette, right.current_color_mask,
              node_current_data(hdt, right), output);
        } else if (node.output_kind == HdtOutputKind::C2) {
          full_bdp_ig_tt_color_mask_into(
              hdt.palette, right.current_color_mask,
              node_current_data(hdt, right), output);
        } else if (node.output_kind == HdtOutputKind::C1U) {
          full_bdp_ig_tf_color_mask_into(
              hdt.palette, right.current_color_mask,
              node_current_data(hdt, right), output);
        } else {
          std::array<double, kFullBdpMaxDimension> c1_u;
          full_bdp_ig_tf_color_mask_into(
              hdt.palette, right.current_color_mask,
              node_current_data(hdt, right), c1_u.data());
          full_bdp_reverse_c1_color_mask_into(
              hdt.palette, right.current_color_mask, c1_u.data(),
              output);
        }
        break;
      case HdtOperation::CC:
        full_bdp_cc_color_mask_into(
            hdt.palette,
            descriptor_for_kind(left.output_kind).first,
            descriptor_for_kind(right.output_kind).first,
            left.current_color_mask, right.current_color_mask,
            node_current_data(hdt, left), node_current_data(hdt, right),
            output);
        break;
      case HdtOperation::Atom: break;
    }
  }
  if (relocated) {
    release_packed_payload(hdt.current_payload, old_payload_class,
                           old_dp_offset);
  }
  node.current_equals_background = next_equals_background;
}

ActiveGraph initial_active_graph(const HdtComponentGraph& graph,
                                 std::size_t expanded_atom_count,
                                 Hdt& hdt,
                                 std::size_t& scratch_peak_bytes) {
  ActiveGraph active;
  const int unit_count = graph.component_units();
  if (expanded_atom_count < graph.units.size()) {
    throw std::logic_error("active atom capacity underflow");
  }
  const std::size_t local_background_count =
      expanded_atom_count - graph.units.size();
  const std::size_t component_capacity = expanded_atom_count;
  const std::size_t edge_capacity = checked_hdt_size_add(
      graph.edges.size(), local_background_count,
      "active edge capacity overflow");
  (void)checked_hdt_int_size(
      component_capacity, "active graph exceeds node domain");
  (void)checked_hdt_int_size(
      edge_capacity, "active graph exceeds edge domain");
  active.components.reserve(component_capacity);
  active.edges.reserve(edge_capacity);
  std::vector<int> unit_to_active(graph.units.size(), -1);
  std::vector<ActivePortList> external_by_unit(graph.units.size());
  for (const HdtExternalPort& external : graph.external_ports) {
    append_external_port(
        external_by_unit[static_cast<std::size_t>(external.unit)],
        external.port);
  }

  for (int unit_index = 0; unit_index < unit_count; ++unit_index) {
    const HdtComponentUnit& unit =
        graph.units[static_cast<std::size_t>(unit_index)];
    HdtExpandedAtom atom;
    atom.component_kind = unit.kind;
    atom.component_unit = unit_index;
    atom.owner = unit.owner;
    atom.source_vertex_id = unit.source_vertex_id;
    atom.weight = unit.weight;
    atom.current_color = unit.current_color;
    atom.background_color = unit.background_color;
    atom.current_color_mask = unit.current_color_mask;
    atom.background_color_mask = unit.background_color_mask;
    atom.port_begin = checked_hdt_int_size(
        hdt.atom_ports.size(), "HDT exceeds its port index domain");
    for (int i = 0; i < unit.incidence_count; ++i) {
      const int edge_index = graph.incidence_edges[static_cast<std::size_t>(
          unit.incidence_begin + i)];
      const HdtIncidenceEdge& edge =
          graph.edges[static_cast<std::size_t>(edge_index)];
      hdt.atom_ports.push_back(edge.a == unit_index ? edge.a_port : edge.b_port);
    }
    const ActivePortList& unit_external =
        external_by_unit[static_cast<std::size_t>(unit_index)];
    for (int i = 0; i < unit_external.count; ++i) {
      hdt.atom_ports.push_back(
          unit_external.values[static_cast<std::size_t>(i)]);
    }
    atom.port_count = checked_hdt_int_size(
                          hdt.atom_ports.size(),
                          "HDT exceeds its port index domain") -
                      atom.port_begin;

    ActiveComponent component;
    component.stable_order = unit.stable_id * 2u;
    component.source_center = unit.source_vertex_id;
    component.external_ports =
        external_by_unit[static_cast<std::size_t>(unit_index)];
    if (unit.kind == HdtComponentKind::IAnchor) {
      component.kind = ActiveKind::I;
    } else if (unit.kind == HdtComponentKind::GOpaque) {
      component.kind = ActiveKind::G;
    } else {
      component.kind = ActiveKind::C;
    }
    component.hdt_root = add_component_atom(
        hdt, std::move(atom), graph, unit_index, component.source_center);
    if (unit_index == graph.root_unit &&
        unit.kind == HdtComponentKind::IAnchor) {
      hdt.root_anchor_atom_node = component.hdt_root;
    }
    unit_to_active[static_cast<std::size_t>(unit_index)] =
        checked_hdt_int_size(active.components.size(),
                             "active graph exceeds node domain");
    active.components.push_back(component);
  }

  for (const HdtIncidenceEdge& edge : graph.edges) {
    active.edges.push_back(ActiveEdge{
        unit_to_active[static_cast<std::size_t>(edge.a)],
        unit_to_active[static_cast<std::size_t>(edge.b)], edge.a_port,
        edge.b_port});
  }
  active.root = unit_to_active[static_cast<std::size_t>(graph.root_unit)];

  for (int unit_index = 0; unit_index < unit_count; ++unit_index) {
    const HdtComponentUnit& unit =
        graph.units[static_cast<std::size_t>(unit_index)];
    if (unit.kind != HdtComponentKind::IAnchor ||
        unit.local_background_g == nullptr) {
      continue;
    }
    HdtExpandedAtom atom;
    atom.component_kind = HdtComponentKind::GOpaque;
    atom.component_unit = unit_index;
    atom.source_vertex_id = unit.source_vertex_id;
    atom.port_begin = checked_hdt_int_size(
        hdt.atom_ports.size(), "HDT exceeds its port index domain");
    hdt.atom_ports.push_back(HdtPortKey{
        HdtPortKind::Center, unit.source_vertex_id,
        HdtPortOrientation::Center});
    atom.port_count = 1;
    ActiveComponent local;
    local.kind = ActiveKind::G;
    local.source_center = unit.source_vertex_id;
    local.stable_order = unit.stable_id * 2u + 1u;
    local.hdt_root = add_atom(hdt, std::move(atom), *unit.local_background_g,
                              HdtOutputKind::G,
                              unit.source_vertex_id);
    const int local_index = checked_hdt_int_size(
        active.components.size(), "active graph exceeds node domain");
    active.components.push_back(local);
    const int anchor = unit_to_active[static_cast<std::size_t>(unit_index)];
    const HdtPortKey center{HdtPortKind::Center,
                                  unit.source_vertex_id,
                                  HdtPortOrientation::Center};
    active.edges.push_back(ActiveEdge{anchor, local_index, center, center});
  }
  std::size_t initial_scratch = active_graph_capacity_bytes(active);
  add_vector_capacity_bytes(initial_scratch, unit_to_active);
  add_vector_capacity_bytes(initial_scratch, external_by_unit);
  scratch_peak_bytes = std::max(scratch_peak_bytes, initial_scratch);
  return active;
}

MutableActiveTree make_mutable_active_tree(
    const ActiveGraph& active,
    const ActiveAdjacency& adjacency,
    Hdt& hdt,
    RoundWorkspace& workspace) {
  rooted_active_tree(active, adjacency, workspace.parent,
                     workspace.parent_edge, workspace.preorder,
                     workspace.traversal_stack);
  MutableActiveTree tree;
  tree.nodes.resize(active.components.size());
  tree.external_ports.reserve(active.components.size());
  tree.links.reserve(active.edges.size());
  tree.root = active.root;
  tree.active_nodes = checked_hdt_int_size(
      active.components.size(), "mutable tree exceeds node domain");
  for (int node = 0; node < tree.active_nodes; ++node) {
    const ActiveComponent& component =
        active.components[static_cast<std::size_t>(node)];
    MutableActiveNode& mutable_node =
        tree.nodes[static_cast<std::size_t>(node)];
    mutable_node.hdt_root = component.hdt_root;
    mutable_node.source_center = component.source_center;
    mutable_node.kind = component.kind;
    mutable_node.parent =
        workspace.parent[static_cast<std::size_t>(node)];
    tree.external_ports.push_back(component.external_ports);
  }
  for (int child : workspace.preorder) {
    const int parent = workspace.parent[static_cast<std::size_t>(child)];
    if (parent < 0) continue;
    const int edge_index =
        workspace.parent_edge[static_cast<std::size_t>(child)];
    const ActiveEdge& edge = active.edges[static_cast<std::size_t>(edge_index)];
    const int link_index = checked_hdt_int_size(
        tree.links.size(), "mutable tree exceeds link domain");
    tree.links.emplace_back();
    MutableActiveLink& link = tree.links.back();
    link.child = child;
    link.parent_port = endpoint_port(edge, parent);
    link.child_port = endpoint_port(edge, child);
    MutableActiveNode& parent_node =
        tree.nodes[static_cast<std::size_t>(parent)];
    if (parent_node.first_child < 0) {
      parent_node.first_child = link_index;
    } else {
      const int previous = parent_node.last_child;
      tree.links[static_cast<std::size_t>(previous)].next = link_index;
      tree.links[static_cast<std::size_t>(link_index)].previous = previous;
    }
    parent_node.last_child = link_index;
    tree.nodes[static_cast<std::size_t>(child)].parent_link = link_index;
  }

  auto set_output = [&](int active_index) {
    const MutableActiveNode& active_node =
        tree.nodes[static_cast<std::size_t>(active_index)];
    HdtNode& node = hdt.nodes[static_cast<std::size_t>(
        active_node.hdt_root)];
    node.output_port_count = 0;
    node.source_center_key = active_node.source_center;
    if (active_node.kind == ActiveKind::I) {
      node.output_kind = HdtOutputKind::I;
      return;
    }
    auto append = [&](const HdtPortKey& port) {
      if (node.output_port_count >= node.output_ports.size()) {
        throw std::logic_error("mutable HDT component exceeds two ports");
      }
      node.output_ports[node.output_port_count++] = port;
    };
    const ActivePortList& external_ports =
        tree.external_ports[static_cast<std::size_t>(active_index)];
    for (int i = 0; i < external_ports.count; ++i) {
      append(external_ports.values[static_cast<std::size_t>(i)]);
    }
    if (active_node.parent_link >= 0) {
      append(tree.links[static_cast<std::size_t>(active_node.parent_link)]
                 .child_port);
    }
    for (int link = active_node.first_child; link >= 0;
         link = tree.links[static_cast<std::size_t>(link)].next) {
      append(tree.links[static_cast<std::size_t>(link)].parent_port);
    }
  };
  for (int node = 0; node < static_cast<int>(tree.nodes.size()); ++node) {
    set_output(node);
  }
  return tree;
}

void mutable_unlink_child(MutableActiveTree& tree, int parent, int link) {
  MutableActiveLink& removed = tree.links[static_cast<std::size_t>(link)];
  const int previous = removed.previous;
  const int next = removed.next;
  if (previous < 0) {
    tree.nodes[static_cast<std::size_t>(parent)].first_child = next;
  } else {
    tree.links[static_cast<std::size_t>(previous)].next = next;
  }
  if (next < 0) {
    tree.nodes[static_cast<std::size_t>(parent)].last_child = previous;
  } else {
    tree.links[static_cast<std::size_t>(next)].previous = previous;
  }
  removed.previous = -1;
  removed.next = -1;
}

void set_mutable_hdt_output(Hdt& hdt,
                            const MutableActiveTree& tree,
                            int active_index) {
  const MutableActiveNode& active_node =
      tree.nodes[static_cast<std::size_t>(active_index)];
  HdtNode& node =
      hdt.nodes[static_cast<std::size_t>(active_node.hdt_root)];
  node.output_port_count = 0;
  node.source_center_key = active_node.source_center;
  if (active_node.kind == ActiveKind::I) {
    node.output_kind = HdtOutputKind::I;
    return;
  }
  auto append = [&](const HdtPortKey& port) {
    if (node.output_port_count >= node.output_ports.size()) {
      throw std::logic_error("mutable HDT component exceeds two ports");
    }
    node.output_ports[node.output_port_count++] = port;
  };
  const ActivePortList& external_ports =
      tree.external_ports[static_cast<std::size_t>(active_index)];
  for (int i = 0; i < external_ports.count; ++i) {
    append(external_ports.values[static_cast<std::size_t>(i)]);
  }
  if (active_node.parent_link >= 0) {
    append(tree.links[static_cast<std::size_t>(active_node.parent_link)]
               .child_port);
  }
  for (int link = active_node.first_child; link >= 0;
       link = tree.links[static_cast<std::size_t>(link)].next) {
    append(tree.links[static_cast<std::size_t>(link)].parent_port);
  }
}

void execute_mutable_round(MutableActiveTree& tree,
                           Hdt& hdt,
                           int round,
                           RoundWorkspace& workspace) {
  std::vector<int>& preorder = workspace.preorder;
  std::size_t active_end = 0;
  const std::size_t previous_size = preorder.size();
  for (std::size_t position = 0; position < previous_size; ++position) {
    const int component = preorder[position];
    MutableActiveNode& current =
        tree.nodes[static_cast<std::size_t>(component)];
    if (current.hdt_root < 0) continue;
    if (component != tree.root &&
        (current.parent < 0 ||
         tree.nodes[static_cast<std::size_t>(current.parent)].hdt_root < 0)) {
      throw std::logic_error("mutable HDT tree disconnected");
    }
    preorder[active_end++] = component;
    if (current.kind != ActiveKind::C || current.parent < 0 ||
        tree.nodes[static_cast<std::size_t>(current.parent)].kind !=
            ActiveKind::I) {
      continue;
    }
    bool down_open = current.first_child >= 0;
    const ActivePortList& current_external_ports =
        tree.external_ports[static_cast<std::size_t>(component)];
    for (int i = 0; i < current_external_ports.count && !down_open; ++i) {
      const HdtPortKey& port = current_external_ports.values[
          static_cast<std::size_t>(i)];
      if (port.kind == HdtPortKind::OrdinaryEdge &&
          port.orientation == HdtPortOrientation::Down) {
        down_open = true;
      }
    }
    if (down_open) continue;

    current.kind = ActiveKind::G;
    current.source_center =
        tree.nodes[static_cast<std::size_t>(current.parent)].source_center;
    HdtNode& root = hdt.nodes[static_cast<std::size_t>(
        current.hdt_root)];
    if (root.post_transform != HdtTransform::None) {
      throw std::logic_error("multiple fused transforms on HDT node");
    }
    root.post_transform = HdtTransform::CToG;
    replace_unary_payload_color_mask(
        hdt, root, HdtOutputKind::G,
        [&](std::uint8_t mask, const double* input, double* output) {
          full_bdp_c_to_g_color_mask_into(
              hdt.palette, mask, input, output);
        });
    root.output_kind = HdtOutputKind::G;
    root.source_center_key = current.source_center;
    const HdtPortKey center{HdtPortKind::Center,
                                  current.source_center,
                                  HdtPortOrientation::Center};
    root.output_port_count = 1;
    root.output_ports[0] = center;
    MutableActiveLink& parent_link =
        tree.links[static_cast<std::size_t>(current.parent_link)];
    parent_link.parent_port = center;
    parent_link.child_port = center;
    ++hdt.stats.c_to_g_ops;
  }
  preorder.resize(active_end);
  if (static_cast<int>(active_end) != tree.active_nodes) {
    throw std::logic_error("mutable HDT preorder lost an active node");
  }

  const int node_slots = checked_hdt_int_size(
      tree.nodes.size(), "mutable tree exceeds node domain");
  std::vector<int>& used_epoch = workspace.mutable_used_epoch;
  std::vector<PlannedGroup>& compositions = workspace.compositions;
  std::vector<PlannedGroup>& ig_compositions = workspace.ig_compositions;
  if (used_epoch.size() != static_cast<std::size_t>(node_slots)) {
    used_epoch.assign(static_cast<std::size_t>(node_slots), 0);
  }
  const auto is_used = [&](int node) {
    return used_epoch[static_cast<std::size_t>(node)] == round;
  };
  auto mark_used = [&](int node) {
    used_epoch[static_cast<std::size_t>(node)] = round;
  };
  compositions.clear();
  compositions.reserve(static_cast<std::size_t>(tree.active_nodes));
  ig_compositions.clear();
  ig_compositions.reserve(static_cast<std::size_t>(tree.active_nodes));

  for (int i_component : preorder) {
    const MutableActiveNode& i_node =
        tree.nodes[static_cast<std::size_t>(i_component)];
    if (i_node.kind != ActiveKind::I) continue;
    int pending_g = -1;
    int g_count = 0;
    int down_open = 0;
    const ActivePortList& i_external_ports =
        tree.external_ports[static_cast<std::size_t>(i_component)];
    for (int i = 0; i < i_external_ports.count; ++i) {
      const HdtPortKey& port = i_external_ports.values[
          static_cast<std::size_t>(i)];
      if (port.kind == HdtPortKind::OrdinaryEdge &&
          port.orientation == HdtPortOrientation::Down) {
        ++down_open;
      }
    }
    for (int link = i_node.first_child; link >= 0;
         link = tree.links[static_cast<std::size_t>(link)].next) {
      const int child = tree.links[static_cast<std::size_t>(link)].child;
      if (tree.nodes[static_cast<std::size_t>(child)].kind != ActiveKind::G) {
        ++down_open;
        continue;
      }
      ++g_count;
      if (pending_g < 0) {
        pending_g = child;
      } else {
        mark_used(pending_g);
        mark_used(child);
        compositions.push_back(
            PlannedGroup{PlannedOperation::GG, pending_g, child});
        pending_g = -1;
      }
    }
    if (g_count == 1 && down_open < 2 && !is_used(pending_g)) {
      mark_used(i_component);
      mark_used(pending_g);
      ig_compositions.push_back(
          PlannedGroup{PlannedOperation::IG, i_component, pending_g});
    }
  }
  compositions.insert(compositions.end(), ig_compositions.begin(),
                      ig_compositions.end());

  for (int top : preorder) {
    if (tree.nodes[static_cast<std::size_t>(top)].kind != ActiveKind::C ||
        is_used(top)) {
      continue;
    }
    const int parent = tree.nodes[static_cast<std::size_t>(top)].parent;
    if (parent >= 0 &&
        tree.nodes[static_cast<std::size_t>(parent)].kind == ActiveKind::C) {
      continue;
    }
    workspace.group_scratch.clear();
    int cursor = top;
    while (cursor >= 0 &&
           tree.nodes[static_cast<std::size_t>(cursor)].kind == ActiveKind::C) {
      workspace.group_scratch.push_back(cursor);
      int next_c = -1;
      for (int link = tree.nodes[static_cast<std::size_t>(cursor)].first_child;
           link >= 0; link = tree.links[static_cast<std::size_t>(link)].next) {
        const int child = tree.links[static_cast<std::size_t>(link)].child;
        if (tree.nodes[static_cast<std::size_t>(child)].kind == ActiveKind::C) {
          next_c = child;
          break;
        }
      }
      cursor = next_c;
    }
    for (int i = 0; i + 1 < static_cast<int>(workspace.group_scratch.size());
         i += 2) {
      const int higher = workspace.group_scratch[static_cast<std::size_t>(i)];
      const int lower =
          workspace.group_scratch[static_cast<std::size_t>(i + 1)];
      if (is_used(higher) || is_used(lower)) {
        continue;
      }
      mark_used(higher);
      mark_used(lower);
      compositions.push_back(
          PlannedGroup{PlannedOperation::CC, higher, lower});
    }
  }

  if (compositions.empty()) {
    throw std::logic_error("mutable SODA round made no progress");
  }

  for (const PlannedGroup& group : compositions) {
    MutableActiveNode& a_node =
        tree.nodes[static_cast<std::size_t>(group.a)];
    MutableActiveNode& b_node =
        tree.nodes[static_cast<std::size_t>(group.b)];
    const int a_hdt_root = a_node.hdt_root;
    const int b_hdt_root = b_node.hdt_root;
    const int a_source_center = a_node.source_center;
    const ActivePortList a_external_ports =
        tree.external_ports[static_cast<std::size_t>(group.a)];
    const ActivePortList b_external_ports =
        tree.external_ports[static_cast<std::size_t>(group.b)];
    ActiveKind result_kind = ActiveKind::I;
    int result_source_center = -1;
    ActivePortList result_external_ports;
    const int result_hdt_root = checked_hdt_int_size(
        hdt.nodes.size(), "HDT exceeds its node index domain");
    hdt.nodes.emplace_back();
    HdtNode& node = hdt.nodes.back();
    node.atom_count = hdt.nodes[static_cast<std::size_t>(a_hdt_root)].atom_count +
                      hdt.nodes[static_cast<std::size_t>(b_hdt_root)].atom_count;
    node.height = 1 + std::max(
        hdt.nodes[static_cast<std::size_t>(a_hdt_root)].height,
        hdt.nodes[static_cast<std::size_t>(b_hdt_root)].height);

    if (group.operation == PlannedOperation::GG) {
      if (a_node.parent < 0 || a_node.parent != b_node.parent ||
          a_source_center != b_node.source_center) {
        throw std::logic_error("mutable GG operands do not share a center");
      }
      result_kind = ActiveKind::G;
      result_source_center = a_source_center;
      evaluate_binary_payload_color_mask(
          hdt, node, HdtOutputKind::G, a_hdt_root, b_hdt_root,
          [&](std::uint8_t left_mask, std::uint8_t right_mask,
              const double* left, const double* right, double* output) {
            full_bdp_gg_color_mask_into(
                hdt.palette, left_mask, right_mask, left, right, output);
          });
      result_external_ports = a_external_ports;
      for (int i = 0; i < b_external_ports.count; ++i) {
        append_external_port(
            result_external_ports,
            b_external_ports.values[static_cast<std::size_t>(i)], true);
      }
      node.operation = HdtOperation::GG;
      node.left = a_hdt_root;
      node.right = b_hdt_root;
      ++hdt.stats.gg_ops;
      mutable_unlink_child(tree, a_node.parent, b_node.parent_link);
    } else if (group.operation == PlannedOperation::IG) {
      if (b_node.parent != group.a) {
        throw std::logic_error("mutable IG operands are not parent/child");
      }
      result_kind = ActiveKind::C;
      result_source_center = a_source_center;
      ActivePortList open_ports = a_external_ports;
      if (a_node.parent_link >= 0) {
        append_external_port(
            open_ports,
            tree.links[static_cast<std::size_t>(a_node.parent_link)]
                .child_port);
      }
      for (int link = a_node.first_child; link >= 0;
           link = tree.links[static_cast<std::size_t>(link)].next) {
        if (tree.links[static_cast<std::size_t>(link)].child == group.b) {
          continue;
        }
        append_external_port(
            open_ports,
            tree.links[static_cast<std::size_t>(link)].parent_port);
      }
      if (open_ports.count == 0) {
        evaluate_unary_payload_color_mask(
            hdt, node, HdtOutputKind::Closed, b_hdt_root,
            [&](std::uint8_t mask, const double* input, double* output) {
              full_bdp_ig_ff_color_mask_into(
                  hdt.palette, mask, input, output);
            });
      } else if (open_ports.count == 1 &&
                 open_ports.values[0].orientation ==
                     HdtPortOrientation::Up) {
        evaluate_unary_payload_color_mask(
            hdt, node, HdtOutputKind::C1U, b_hdt_root,
            [&](std::uint8_t mask, const double* input, double* output) {
              full_bdp_ig_tf_color_mask_into(
                  hdt.palette, mask, input, output);
            });
      } else if (open_ports.count == 2) {
        evaluate_unary_payload_color_mask(
            hdt, node, HdtOutputKind::C2, b_hdt_root,
            [&](std::uint8_t mask, const double* input, double* output) {
              full_bdp_ig_tt_color_mask_into(
                  hdt.palette, mask, input, output);
            });
      } else {
        const FullBdpDual value = ig_output(
            hdt.palette, hdt_node_dp(hdt, b_hdt_root), open_ports);
        store_node_payload(hdt, node, value,
                           output_kind_of(value.current));
      }
      result_external_ports = a_external_ports;
      node.operation = HdtOperation::IG;
      node.left = a_hdt_root;
      node.right = b_hdt_root;
      ++hdt.stats.ig_ops;
      mutable_unlink_child(tree, group.a, b_node.parent_link);
    } else {
      if (b_node.parent != group.a) {
        throw std::logic_error("mutable CC operands are not parent/child");
      }
      const MutableActiveLink& shared =
          tree.links[static_cast<std::size_t>(b_node.parent_link)];
      const bool a_higher =
          shared.parent_port.orientation == HdtPortOrientation::Down &&
          shared.child_port.orientation == HdtPortOrientation::Up;
      if (!a_higher) {
        throw std::logic_error("mutable CC ordinary-port orientation mismatch");
      }
      result_kind = ActiveKind::C;
      const HdtOutputKind cc_output =
          hdt.nodes[static_cast<std::size_t>(a_hdt_root)].output_kind ==
                      HdtOutputKind::C1D &&
                  hdt.nodes[static_cast<std::size_t>(b_hdt_root)].output_kind ==
                      HdtOutputKind::C1U
              ? HdtOutputKind::Closed
              : (hdt.nodes[static_cast<std::size_t>(a_hdt_root)].output_kind ==
                             HdtOutputKind::C1D
                     ? HdtOutputKind::C1D
                     : (hdt.nodes[static_cast<std::size_t>(b_hdt_root)]
                                    .output_kind ==
                                HdtOutputKind::C1U
                            ? HdtOutputKind::C1U
                            : HdtOutputKind::C2));
      const FullBdpInterface higher_interface = descriptor_for_kind(
          hdt.nodes[static_cast<std::size_t>(a_hdt_root)].output_kind).first;
      const FullBdpInterface lower_interface = descriptor_for_kind(
          hdt.nodes[static_cast<std::size_t>(b_hdt_root)].output_kind).first;
      evaluate_binary_payload_color_mask(
          hdt, node, cc_output, a_hdt_root, b_hdt_root,
          [&](std::uint8_t higher_mask, std::uint8_t lower_mask,
              const double* higher, const double* lower, double* output) {
            full_bdp_cc_color_mask_into(
                hdt.palette, higher_interface, lower_interface,
                higher_mask, lower_mask, higher, lower, output);
          });
      result_external_ports = a_external_ports;
      for (int i = 0; i < b_external_ports.count; ++i) {
        append_external_port(
            result_external_ports,
            b_external_ports.values[static_cast<std::size_t>(i)]);
      }
      result_source_center = a_source_center;
      node.operation = HdtOperation::CC;
      node.left = a_hdt_root;
      node.right = b_hdt_root;
      ++hdt.stats.cc_ops;

      const int removed_link = b_node.parent_link;
      MutableActiveLink& removed =
          tree.links[static_cast<std::size_t>(removed_link)];
      const int previous = removed.previous;
      const int after = removed.next;
      const int adopted = b_node.first_child;
      if (adopted < 0) {
        if (previous < 0) a_node.first_child = after;
        else tree.links[static_cast<std::size_t>(previous)].next = after;
        if (after < 0) a_node.last_child = previous;
        else tree.links[static_cast<std::size_t>(after)].previous = previous;
      } else {
        if (previous < 0) a_node.first_child = adopted;
        else tree.links[static_cast<std::size_t>(previous)].next = adopted;
        tree.links[static_cast<std::size_t>(adopted)].previous = previous;
        const int tail = b_node.last_child;
        for (int link = adopted; link >= 0;
             link = tree.links[static_cast<std::size_t>(link)].next) {
          const int child = tree.links[static_cast<std::size_t>(link)].child;
          tree.nodes[static_cast<std::size_t>(child)].parent = group.a;
        }
        tree.links[static_cast<std::size_t>(tail)].next = after;
        if (after < 0) a_node.last_child = tail;
        else tree.links[static_cast<std::size_t>(after)].previous = tail;
      }
      removed.previous = -1;
      removed.next = -1;
      b_node.first_child = -1;
      b_node.last_child = -1;
    }

    hdt.nodes[static_cast<std::size_t>(node.left)].parent = result_hdt_root;
    hdt.nodes[static_cast<std::size_t>(node.right)].parent = result_hdt_root;
    note_node_shape(hdt, node);
    a_node.hdt_root = result_hdt_root;
    a_node.source_center = result_source_center;
    a_node.kind = result_kind;
    tree.external_ports[static_cast<std::size_t>(group.a)] =
        result_external_ports;
    b_node.hdt_root = -1;
    b_node.parent = -1;
    b_node.parent_link = -1;
    --tree.active_nodes;
    set_mutable_hdt_output(hdt, tree, group.a);
  }

  const long long removed = static_cast<long long>(compositions.size());
  const long long before = tree.active_nodes + removed;
  if (9LL * removed < before) {
    throw std::logic_error(
        "mutable SODA round violated the 1-in-9 shrink bound");
  }
}

void finalize_hdt_metadata(Hdt& hdt,
                           int owner_count,
                           HdtOwnerIndexView owner_index) {
  const bool by_position = owner_index.owner_position != nullptr;
  if (by_position &&
      (owner_index.position_begin < 0 ||
       owner_index.position_end < owner_index.position_begin)) {
    throw std::invalid_argument("invalid owner-position index");
  }
  const int slot_count = by_position
                             ? owner_index.position_end -
                                   owner_index.position_begin
                             : owner_count;
  if (slot_count < 0) {
    throw std::invalid_argument("negative owner index domain");
  }
  hdt.owner_index_by_position = by_position;
  hdt.owner_position_begin =
      by_position ? owner_index.position_begin : 0;
  auto owner_slot = [&](int owner) {
    if (!by_position) return owner;
    if (owner < 0 ||
        static_cast<std::size_t>(owner) >=
            owner_index.owner_position->size()) {
      return -1;
    }
    const int position = (*owner_index.owner_position)[
        static_cast<std::size_t>(owner)];
    return position >= owner_index.position_begin &&
                   position < owner_index.position_end
               ? position - owner_index.position_begin
               : -1;
  };
  hdt.dirty.assign(hdt.nodes.size(), 0);
  hdt.recompute_stack.clear();
  hdt.recompute_stack.reserve(static_cast<std::size_t>(checked_hdt_int_mul_add(
      hdt.stats.max_update_height, 2, 4,
      "recompute stack capacity overflow")));
  hdt.owner_offsets.assign(
      checked_hdt_size_add(static_cast<std::size_t>(slot_count), 1u,
                           "owner-offset capacity overflow"),
      0);
  const int atom_count = checked_hdt_int_size(
      hdt.atoms.size(), "HDT exceeds its atom index domain");
  for (int node = 0; node < atom_count; ++node) {
    if (hdt.nodes[static_cast<std::size_t>(node)].operation !=
            HdtOperation::Atom ||
        hdt.nodes[static_cast<std::size_t>(node)].atom_descriptor != node) {
      throw std::logic_error("atom/node prefix invariant failed");
    }
    const int owner = hdt.nodes[static_cast<std::size_t>(node)].live_owner;
    const int slot = owner_slot(owner);
    if (slot >= 0 && slot < slot_count) {
      ++hdt.owner_offsets[static_cast<std::size_t>(slot) + 1u];
    }
  }
  for (int slot = 0; slot < slot_count; ++slot) {
    hdt.owner_offsets[static_cast<std::size_t>(slot) + 1u] +=
        hdt.owner_offsets[static_cast<std::size_t>(slot)];
  }
  hdt.owner_atom_nodes.assign(
      static_cast<std::size_t>(hdt.owner_offsets.back()), -1);
  for (int node = atom_count - 1; node >= 0; --node) {
    const int owner = hdt.nodes[static_cast<std::size_t>(node)].live_owner;
    const int slot = owner_slot(owner);
    if (slot < 0 || slot >= slot_count) continue;
    const int position = --hdt.owner_offsets[
        static_cast<std::size_t>(slot) + 1u];
    hdt.owner_atom_nodes[static_cast<std::size_t>(position)] = node;
  }
  for (int slot = 1; slot < slot_count; ++slot) {
    hdt.owner_offsets[static_cast<std::size_t>(slot)] =
        hdt.owner_offsets[static_cast<std::size_t>(slot) + 1u];
  }
  hdt.owner_offsets[static_cast<std::size_t>(slot_count)] =
      checked_hdt_int_size(hdt.owner_atom_nodes.size(),
                           "owner index exceeds integer domain");
}

void set_error(std::string* error, const std::string& message) {
  if (error != nullptr) *error = message;
}

}  // namespace

static void build_hdt_impl(Hdt& hdt,
                              const HdtComponentGraph& graph,
                              bool require_closed_root,
                              bool materialize_update_metadata,
                              bool background_alias_current,
                              HdtOwnerIndexView owner_index) {
  const long long scalar_fmas_before = full_bdp_thread_scalar_fmas();
#ifndef NDEBUG
  std::string graph_error;
  if (!validate_component_graph(graph, &graph_error)) {
    throw std::invalid_argument("invalid component graph: " + graph_error);
  }
#endif
  hdt.root = -1;
  hdt.root_anchor_atom_node = -1;
  hdt.nodes.clear();
  hdt.current_payload.clear();
  hdt.background_payload.clear();
  hdt.atoms.clear();
  hdt.atom_ports.clear();
  hdt.owner_offsets.clear();
  hdt.owner_atom_nodes.clear();
  hdt.recompute_stack.clear();
  hdt.owner_position_begin = 0;
  hdt.owner_index_by_position = false;
  hdt.dirty.clear();
  hdt.stats = HdtBuildStats{};
  hdt.palette = graph.palette;
  hdt.closed_root_required = require_closed_root;
  hdt.background_alias_current = background_alias_current;
  hdt.background_restricted = graph.background_restricted;
  hdt.stats.component_units = graph.component_units();
  if (graph.units.empty()) {
    // The empty quotient has A=H=R=0.  It is the total sparse-bank base case,
    // so both structural predicates are vacuously satisfied even though the
    // nonempty identity H=2A-1 is not applicable.
    hdt.stats.node_count_ok = true;
    hdt.stats.atom_bound_ok = true;
    hdt.stats.scalar_fmas =
        full_bdp_thread_scalar_fmas() - scalar_fmas_before;
    return;
  }

  const std::size_t component_units = graph.units.size();
  std::size_t initial_atoms = component_units;
  for (const HdtComponentUnit& unit : graph.units) {
    if (unit.kind == HdtComponentKind::IAnchor &&
        unit.local_background_g != nullptr) {
      initial_atoms = checked_hdt_size_add(
          initial_atoms, 1u, "HDT atom count overflow");
    }
  }
  const int initial_atom_count = checked_hdt_int_size(
      initial_atoms, "HDT exceeds its atom index domain");
  const int node_capacity = checked_hdt_int_mul_add(
      initial_atom_count, 2, -1,
      "HDT exceeds its node index domain");
  const std::size_t local_background_count = initial_atoms - component_units;
  const std::size_t atom_port_capacity = checked_hdt_size_add(
      checked_hdt_size_add(
          checked_hdt_size_mul(graph.edges.size(), 2u,
                               "HDT port capacity overflow"),
          graph.external_ports.size(), "HDT port capacity overflow"),
      local_background_count, "HDT port capacity overflow");
  (void)checked_hdt_int_size(
      atom_port_capacity, "HDT exceeds its port index domain");
  // The builder emits exactly A atoms and A-1 binary nodes. Payload arenas
  // grow independently by exact active-coordinate width.
  hdt.atoms.reserve(initial_atoms);
  hdt.nodes.reserve(static_cast<std::size_t>(node_capacity));
  hdt.atom_ports.reserve(atom_port_capacity);

  std::size_t build_scratch_peak_bytes = 0;
  ActiveGraph active = initial_active_graph(
      graph, initial_atoms, hdt, build_scratch_peak_bytes);
  hdt.stats.initial_atoms = checked_hdt_int_size(
      active.components.size(), "HDT exceeds its atom index domain");
  if (hdt.stats.initial_atoms != initial_atom_count) {
    throw std::logic_error("HDT atom count changed during expansion");
  }
  RoundWorkspace round_workspace;
  ActiveAdjacency adjacency;
  build_active_adjacency(active, adjacency,
                         round_workspace.adjacency_cursor);
  MutableActiveTree mutable_tree = make_mutable_active_tree(
      active, adjacency, hdt, round_workspace);

  int round = 0;
  while (mutable_tree.active_nodes > 1) {
    ++round;
    if (static_cast<long long>(round) >
        4LL * static_cast<long long>(hdt.stats.initial_atoms) + 16LL) {
      throw std::logic_error("SODA builder exceeded round bound");
    }
    execute_mutable_round(mutable_tree, hdt, round, round_workspace);
  }
  build_scratch_peak_bytes = std::max(
      build_scratch_peak_bytes,
      live_builder_scratch_bytes(active, adjacency, mutable_tree,
                                 round_workspace));
  hdt.stats.build_scratch_peak_bytes = build_scratch_peak_bytes;
  hdt.stats.rounds = round;
  hdt.root = mutable_tree.nodes[static_cast<std::size_t>(mutable_tree.root)]
                 .hdt_root;
  HdtNode& root = hdt.nodes[static_cast<std::size_t>(hdt.root)];
  if (require_closed_root && root.output_kind == HdtOutputKind::C1U) {
    root.post_transform = HdtTransform::CloseRoot;
    replace_unary_payload_color_mask(
        hdt, root, HdtOutputKind::Closed,
        [&](std::uint8_t mask, const double* input, double* output) {
          full_bdp_close_color_mask_into(
              hdt.palette, mask, input, output);
        });
    root.output_kind = HdtOutputKind::Closed;
    root.output_port_count = 0;
  }
  if (require_closed_root &&
      root.output_kind != HdtOutputKind::Closed) {
    throw std::logic_error("HDT root is not Closed");
  }
  if (materialize_update_metadata) {
    finalize_hdt_metadata(hdt, graph.owner_count, owner_index);
  }

  hdt.stats.reachable_nodes = checked_hdt_int_size(
      hdt.nodes.size(), "HDT exceeds its node index domain");
  hdt.stats.node_count_ok =
      static_cast<long long>(hdt.stats.reachable_nodes) ==
      2LL * static_cast<long long>(hdt.stats.initial_atoms) - 1LL;
  hdt.stats.atom_bound_ok =
      static_cast<long long>(hdt.stats.initial_atoms) <=
      2LL * static_cast<long long>(hdt.stats.component_units);
  if (!hdt.stats.node_count_ok || !hdt.stats.atom_bound_ok ||
      hdt.stats.max_recompute_degree > 2 ||
      hdt.stats.local_balance_violations != 0) {
    throw std::logic_error(
        "HDT construction violated a structural invariant");
  }
  if (materialize_update_metadata) {
#ifndef NDEBUG
    std::string error;
    if (!validate_hdt(graph, hdt, &error)) {
      throw std::logic_error("invalid HDT: " + error);
    }
#endif
  }
  hdt.stats.scalar_fmas =
      full_bdp_thread_scalar_fmas() - scalar_fmas_before;
}

void build_hdt_into(Hdt& hdt,
                              const HdtComponentGraph& graph,
                              bool require_closed_root,
                              HdtOwnerIndexView owner_index) {
  build_hdt_impl(hdt, graph, require_closed_root, true, false,
                           owner_index);
}

void build_frozen_fold_hdt_into(
    Hdt& hdt,
    const HdtComponentGraph& graph,
    bool require_closed_root) {
  build_hdt_impl(hdt, graph, require_closed_root, false, true, {});
}

Hdt build_hdt(const HdtComponentGraph& graph,
                                 bool require_closed_root,
                                 HdtOwnerIndexView owner_index) {
  Hdt hdt;
  build_hdt_into(hdt, graph, require_closed_root, owner_index);
  return hdt;
}

HdtState build_hdt_state(
    HdtComponentGraph graph,
    HdtOwnerIndexView owner_index) {
  HdtState state;
  state.components = std::move(graph);
  state.hdt = build_hdt(state.components, true, owner_index);
  return state;
}

FullBdpVec hdt_node_current(const Hdt& hdt, int node) {
  if (node < 0 || node >= static_cast<int>(hdt.nodes.size())) {
    throw std::out_of_range("HDT current payload node out of range");
  }
  const HdtNode& value =
      hdt.nodes[static_cast<std::size_t>(node)];
  return materialize_view(node_current_view(hdt, value), hdt.palette,
                          value.current_color_mask);
}

FullBdpVec hdt_node_background(const Hdt& hdt, int node) {
  if (node < 0 || node >= static_cast<int>(hdt.nodes.size())) {
    throw std::out_of_range(
        "HDT background payload node out of range");
  }
  const HdtNode& value =
      hdt.nodes[static_cast<std::size_t>(node)];
  return materialize_view(node_background_view(hdt, value), hdt.palette,
                          value.background_color_mask);
}

FullBdpDual hdt_node_dp(const Hdt& hdt, int node) {
  if (node < 0 || node >= static_cast<int>(hdt.nodes.size())) {
    throw std::out_of_range("HDT dual payload node out of range");
  }
  const HdtNode& value = hdt.nodes[static_cast<std::size_t>(node)];
  if (value.output_kind == HdtOutputKind::I) return FullBdpDual{};
  return FullBdpDual{materialize_view(node_current_view(hdt, value),
                                     hdt.palette,
                                     value.current_color_mask),
                     materialize_view(node_background_view(hdt, value),
                                     hdt.palette,
                                     value.background_color_mask),
                     value.current_equals_background};
}

bool validate_hdt(const HdtComponentGraph& graph,
                            const Hdt& hdt,
                            std::string* error) {
  if (graph.units.empty()) return hdt.root == -1 && hdt.nodes.empty();
  if (hdt.root < 0 || hdt.root >= static_cast<int>(hdt.nodes.size())) {
    set_error(error, "invalid root");
    return false;
  }
  if (!hdt.stats.node_count_ok) {
    set_error(error, "H != 2A-1");
    return false;
  }
  if (!hdt.stats.atom_bound_ok) {
    set_error(error, "A > 2R");
    return false;
  }
  if (hdt.stats.max_recompute_degree > 2) {
    set_error(error, "nonbinary recompute node");
    return false;
  }
  if (hdt.stats.local_balance_violations != 0) {
    set_error(error, "local balance violation");
    return false;
  }
  if (hdt.root_anchor_atom_node >= 0) {
    if (hdt.root_anchor_atom_node >= static_cast<int>(hdt.nodes.size()) ||
        hdt.nodes[static_cast<std::size_t>(hdt.root_anchor_atom_node)]
                .operation != HdtOperation::Atom) {
      set_error(error, "invalid stored root-anchor atom");
      return false;
    }
    const HdtExpandedAtom& atom = hdt.atoms[static_cast<std::size_t>(
        hdt.nodes[static_cast<std::size_t>(hdt.root_anchor_atom_node)]
            .atom_descriptor)];
    if (atom.component_kind != HdtComponentKind::IAnchor ||
        atom.component_unit != graph.root_unit) {
      set_error(error, "stored root-anchor atom does not name graph root");
      return false;
    }
  }
  const HdtNode& root =
      hdt.nodes[static_cast<std::size_t>(hdt.root)];
  if (hdt.closed_root_required) {
    if (root.output_kind != HdtOutputKind::Closed) {
      set_error(error, "root output is not Closed");
      return false;
    }
  } else if (root.output_kind == HdtOutputKind::I ||
             root.output_port_count != graph.external_ports.size()) {
    set_error(error, "open root output does not match external ports");
    return false;
  }
  if (root.parent != -1) {
    set_error(error, "root has a parent");
    return false;
  }
  std::vector<std::uint8_t> seen(hdt.nodes.size(), 0);
  std::vector<int> stack{hdt.root};
  int reached = 0;
  while (!stack.empty()) {
    const int node = stack.back();
    stack.pop_back();
    if (node < 0 || node >= static_cast<int>(hdt.nodes.size()) ||
        seen[static_cast<std::size_t>(node)] != 0) {
      set_error(error, "node reachability is not one-to-one");
      return false;
    }
    seen[static_cast<std::size_t>(node)] = 1;
    ++reached;
    const HdtNode& value = hdt.nodes[static_cast<std::size_t>(node)];
    if (value.output_kind == HdtOutputKind::I) {
      if (value.dp_offset != -1 || value.background_dp_offset != -1 ||
          value.dp_dimension != 0 || value.current_payload_class != 0xff ||
          value.background_payload_class != 0xff) {
        set_error(error, "I node owns a packed payload");
        return false;
      }
    } else {
      const FullBdpInterface interface =
          descriptor_for_kind(value.output_kind).first;
      if (value.dp_dimension != full_bdp_dimension(hdt.palette, interface)) {
        set_error(error, "packed payload dimension mismatch");
        return false;
      }
      const std::uint8_t expected_current_class = packed_payload_class_index(
          packed_payload_width(hdt.palette, interface,
                               value.current_color_mask));
      if (value.current_payload_class != expected_current_class ||
          value.dp_offset < 0) {
        set_error(error, "packed current payload class mismatch");
        return false;
      }
      const int current_width = packed_payload_width(expected_current_class);
      const auto& current_values = hdt.current_payload
                                       .classes[expected_current_class]
                                       .values;
      if (value.dp_offset % current_width != 0 ||
          static_cast<std::size_t>(value.dp_offset) + current_width >
              current_values.size()) {
        set_error(error, "packed current payload offset mismatch");
        return false;
      }
      if (value.background_dp_offset < 0) {
        if (value.background_payload_class != 0xff) {
          set_error(error, "aliased background has a payload class");
          return false;
        }
      } else {
        const std::uint8_t expected_background_class =
            packed_payload_class_index(packed_payload_width(
                hdt.palette, interface, value.background_color_mask));
        if (hdt.background_alias_current ||
            value.background_payload_class != expected_background_class) {
          set_error(error, "packed background payload class mismatch");
          return false;
        }
        const int background_width =
            packed_payload_width(expected_background_class);
        const auto& background_values =
            hdt.background_payload.classes[expected_background_class].values;
        if (value.background_dp_offset % background_width != 0 ||
            static_cast<std::size_t>(value.background_dp_offset) +
                    background_width >
                background_values.size()) {
          set_error(error, "packed background payload offset mismatch");
          return false;
        }
      }
    }
    if (value.operation == HdtOperation::Atom) {
      if (value.left >= 0 || value.right >= 0 || value.atom_descriptor < 0 ||
          value.atom_descriptor >= static_cast<int>(hdt.atoms.size())) {
        set_error(error, "invalid atom node");
        return false;
      }
    } else {
      if (value.left < 0 || value.right < 0 ||
          value.left >= static_cast<int>(hdt.nodes.size()) ||
          value.right >= static_cast<int>(hdt.nodes.size()) ||
          hdt.nodes[static_cast<std::size_t>(value.left)].parent != node ||
          hdt.nodes[static_cast<std::size_t>(value.right)].parent != node) {
        set_error(error, "invalid binary parent link");
        return false;
      }
      stack.push_back(value.right);
      stack.push_back(value.left);
    }
  }
  if (reached != static_cast<int>(hdt.nodes.size())) {
    set_error(error, "HDT contains unreachable nodes");
    return false;
  }
  return true;
}

int hdt_recolor_owner_ranges(
    HdtState& state,
    const std::vector<int>& owner_order,
    std::span<const std::pair<int, int>> owner_ranges,
    std::uint8_t color) {
  const long long scalar_fmas_before = full_bdp_thread_scalar_fmas();
  Hdt& hdt = state.hdt;
  if (color >= hdt.palette) {
    throw std::invalid_argument("recolor outside palette");
  }
  const int owner_order_count = checked_hdt_int_size(
      owner_order.size(), "recolor exceeds integer owner domain");
  const int component_count = state.components.component_units();
  const int owner_offset_count = checked_hdt_int_size(
      hdt.owner_offsets.size(),
      "owner index exceeds its integer domain");
  for (const auto& [begin_owner, end_owner] : owner_ranges) {
    if (begin_owner < 0 || end_owner < begin_owner ||
        end_owner > owner_order_count) {
      throw std::invalid_argument("invalid recolor owner range");
    }
  }
  auto recolor_component_unit = [&](int unit_index) {
    if (unit_index < 0 || unit_index >= component_count) {
      throw std::logic_error("live HDT atom lacks its component unit");
    }
    HdtComponentUnit& unit =
        state.components.units[static_cast<std::size_t>(unit_index)];
    if (unit.kind != HdtComponentKind::LiveLeaf) {
      throw std::logic_error("mutable HDT atom is not a live component");
    }
    unit.current_color = color;
    unit.current_color_mask = unit.weight == 0.0
                                  ? 0
                                  : static_cast<std::uint8_t>(1u << color);
    unit.current_equals_background =
        color == unit.background_color || unit.weight == 0.0;
  };
  int changed = 0;
  auto recolor_slot = [&](int slot) {
    if (slot < 0 || slot >= owner_offset_count - 1) {
      return;
    }
    const int begin = hdt.owner_offsets[static_cast<std::size_t>(slot)];
    const int end = hdt.owner_offsets[static_cast<std::size_t>(slot) + 1u];
    for (int p = begin; p < end; ++p) {
      const int node_index =
          hdt.owner_atom_nodes[static_cast<std::size_t>(p)];
      HdtNode& node =
          hdt.nodes[static_cast<std::size_t>(node_index)];
      HdtExpandedAtom& atom = hdt.atoms[static_cast<std::size_t>(
          node.atom_descriptor)];
      if (atom.current_color == color) continue;
      atom.current_color = color;
      atom.current_color_mask = atom.weight == 0.0
                                    ? 0
                                    : static_cast<std::uint8_t>(1u << color);
      recolor_component_unit(atom.component_unit);
      ++changed;
      for (int mark = node_index; mark >= 0 &&
                                      hdt.dirty[static_cast<std::size_t>(mark)] == 0;
           mark = hdt.nodes[static_cast<std::size_t>(mark)].parent) {
        hdt.dirty[static_cast<std::size_t>(mark)] = 1;
      }
    }
  };

  if (hdt.owner_index_by_position) {
    if (owner_offset_count == 0) {
      throw std::logic_error("owner-position index has no offsets");
    }
    const int indexed_end = checked_hdt_int_add(
        hdt.owner_position_begin, owner_offset_count - 1,
        "owner-position interval overflow");
    for (const auto& [begin_owner, end_owner] : owner_ranges) {
      const int begin_position =
          std::max(begin_owner, hdt.owner_position_begin);
      const int end_position = std::min(end_owner, indexed_end);
      for (int position = begin_position; position < end_position;
           ++position) {
        recolor_slot(position - hdt.owner_position_begin);
      }
    }
  } else {
    for (const auto& [begin_owner, end_owner] : owner_ranges) {
      for (int owner_index = begin_owner; owner_index < end_owner;
           ++owner_index) {
        recolor_slot(owner_order[static_cast<std::size_t>(owner_index)]);
      }
    }
  }
  int recomputed = 0;
  if (changed != 0) {
    std::vector<int>& stack = hdt.recompute_stack;
    stack.clear();
    stack.push_back(hdt.root);
    while (!stack.empty()) {
      const int token = stack.back();
      stack.pop_back();
      if (token < 0) {
        const int node = ~token;
        recompute_node(hdt, node);
        hdt.dirty[static_cast<std::size_t>(node)] = 0;
        ++recomputed;
        continue;
      }
      const int node = token;
      if (node < 0 || hdt.dirty[static_cast<std::size_t>(node)] == 0) continue;
      stack.push_back(~node);
      const HdtNode& value =
          hdt.nodes[static_cast<std::size_t>(node)];
      if (value.right >= 0 &&
          hdt.dirty[static_cast<std::size_t>(value.right)] != 0) {
        stack.push_back(value.right);
      }
      if (value.left >= 0 &&
          hdt.dirty[static_cast<std::size_t>(value.left)] != 0) {
        stack.push_back(value.left);
      }
    }
  }
  state.recolor_changed_leaves += changed;
  state.recolor_marked_nodes += recomputed;
  state.transition_evals += recomputed;
  state.scalar_fmas +=
      full_bdp_thread_scalar_fmas() - scalar_fmas_before;
  return changed;
}

int hdt_recolor_owner_range(HdtState& state,
                                  const std::vector<int>& owner_order,
                                  int begin_owner,
                                  int end_owner,
                                  std::uint8_t color) {
  const std::pair<int, int> range{begin_owner, end_owner};
  return hdt_recolor_owner_ranges(
      state, owner_order, std::span<const std::pair<int, int>>(&range, 1),
      color);
}

int hdt_recolor_owners(HdtState& state,
                             const std::vector<int>& owners,
                             std::uint8_t color) {
  return hdt_recolor_owner_range(
      state, owners, 0,
      checked_hdt_int_size(owners.size(),
                           "recolor exceeds integer owner domain"),
      color);
}

int hdt_recolor_owner(HdtState& state,
                            int owner,
                            std::uint8_t color) {
  return hdt_recolor_owners(state, std::vector<int>{owner}, color);
}

std::array<double, 3> hdt_root_counts(
    const HdtState& state) {
  if (state.hdt.root < 0) return {0.0, 0.0, 0.0};
  const HdtNode& root =
      state.hdt.nodes[static_cast<std::size_t>(state.hdt.root)];
  return full_bdp_color_mask_readout(
      state.hdt.palette, root.current_color_mask,
      node_current_data(state.hdt, root));
}

int hdt_live_atom_count(const HdtState& state) {
  return state.hdt.owner_offsets.empty() ? 0
                                         : state.hdt.owner_offsets.back();
}

}  // namespace wqfm::internal
