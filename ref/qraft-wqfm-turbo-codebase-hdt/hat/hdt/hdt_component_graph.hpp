#pragma once

#include "hat/hdt/hdt_boundary_dp.hpp"
#include "hat/hdt/hdt_checked.hpp"
#include "hat/hdt/owner_projection.hpp"

#include <array>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

namespace wqfm::internal {

enum class HdtPortKind : std::uint8_t {
  OrdinaryEdge,
  Center,
};

enum class HdtPortOrientation : std::uint8_t {
  Up,
  Down,
  Center,
};

struct HdtPortKey {
  int source_id = -1;
  HdtPortKind kind = HdtPortKind::OrdinaryEdge;
  HdtPortOrientation orientation = HdtPortOrientation::Up;

  constexpr HdtPortKey() = default;
  constexpr HdtPortKey(HdtPortKind kind_value,
                             int source_id_value,
                             HdtPortOrientation orientation_value)
      : source_id(source_id_value),
        kind(kind_value),
        orientation(orientation_value) {}
};

inline bool operator==(const HdtPortKey& left,
                       const HdtPortKey& right) {
  return left.kind == right.kind && left.source_id == right.source_id &&
         left.orientation == right.orientation;
}

enum class HdtComponentKind : std::uint8_t {
  LiveLeaf,
  IAnchor,
  ClosedOpaque,
  C1Opaque,
  C2Opaque,
  GOpaque,
};

struct HdtSourceIndex {
  std::vector<int> parent;
  std::vector<int> parent_edge;
  std::vector<int> child_position;
  std::vector<int> edge_parent;
  std::vector<int> edge_child;
  int root_vertex = -1;
  int edge_slots = 0;
};

struct HdtComponentUnit {
  HdtComponentKind kind = HdtComponentKind::IAnchor;
  std::uint64_t stable_id = 0;
  int owner = -1;
  int source_vertex_id = -1;
  int incidence_begin = 0;
  int incidence_count = 0;
  double weight = 0.0;
  std::uint8_t current_color = 0;
  std::uint8_t background_color = 0;
  std::uint8_t current_color_mask = 0;
  std::uint8_t background_color_mask = 0;
  int dp_offset = -1;
  int background_dp_offset = -1;
  std::uint8_t dp_dimension = 0;
  FullBdpOrientation dp_orientation = FullBdpOrientation::None;
  bool current_equals_background = false;
  std::shared_ptr<const FullBdpDual> local_background_g;
};

struct HdtIncidenceEdge {
  int a = -1;
  int b = -1;
  HdtPortKey a_port;
  HdtPortKey b_port;
};

struct HdtExternalPort {
  int unit = -1;
  HdtPortKey port;
};

struct HdtColorPolicy {
  std::uint8_t default_current = 1;
  std::uint8_t default_background = 0;
  std::vector<std::uint8_t> current_by_owner;
  std::vector<std::uint8_t> background_by_owner;
};

struct HdtComponentGraph {
  int palette = 0;
  int root_unit = -1;
  std::shared_ptr<const HdtSourceIndex> source;
  std::vector<HdtComponentUnit> units;
  std::vector<double> current_payload;
  std::vector<double> background_payload;
  std::vector<HdtIncidenceEdge> edges;
  std::vector<HdtExternalPort> external_ports;
  std::vector<int> incidence_edges;
  int owner_count = 0;
  long long construction_scalar_fmas = 0;
  bool background_restricted = false;

  int component_units() const {
    return checked_hdt_int_size(
        units.size(), "component graph exceeds its integer domain");
  }
};

HdtComponentGraph build_base_component_graph(
    const GeneIndex& gene,
    const OwnerProjection& owners,
    int palette,
    const HdtColorPolicy& colors);

// Builds the same base quotient after first suppressing source leaves that are
// absent from the owner projection. This avoids constructing and immediately
// Contracting zero-mass atoms when a physical child bank is reused.
HdtComponentGraph build_induced_base_component_graph(
    const GeneIndex& gene,
    const OwnerProjection& owners,
    int palette,
    const HdtColorPolicy& colors);

// Returns the stable induced source used by
// build_induced_base_component_graph. `changed` is false when every physical
// leaf is live, in which case callers should continue using `gene`.
GeneIndex induce_gene_for_owner_projection(
    const GeneIndex& gene,
    const OwnerProjection& owners,
    bool& changed);

void rebuild_component_graph_indices(HdtComponentGraph& graph);

bool validate_component_graph(const HdtComponentGraph& graph,
                              std::string* error = nullptr);

void hdt_component_store_dp(HdtComponentGraph& graph,
                            HdtComponentUnit& unit,
                            const FullBdpDual& value);
FullBdpDual hdt_component_dp(const HdtComponentGraph& graph,
                             int unit);
FullBdpVec hdt_component_background(
    const HdtComponentGraph& graph,
    int unit);

std::vector<int> hdt_component_neighbors(
    const HdtComponentGraph& graph,
    int unit);

}  // namespace wqfm::internal
