#pragma once

#include "hat/hdt/hdt_component_graph.hpp"

#include <array>
#include <cstdint>
#include <span>
#include <string>
#include <utility>
#include <vector>

namespace wqfm::internal {

enum class HdtOperation : std::uint8_t {
  Atom,
  GG,
  IG,
  CC,
};

enum class HdtTransform : std::uint8_t {
  None,
  CToG,
  CloseRoot,
};

enum class HdtOutputKind : std::uint8_t {
  I,
  Closed,
  C1U,
  C1D,
  C2,
  G,
};

struct HdtExpandedAtom {
  HdtComponentKind component_kind = HdtComponentKind::IAnchor;
  int component_unit = -1;
  int owner = -1;
  int source_vertex_id = -1;
  int port_begin = 0;
  int port_count = 0;
  double weight = 0.0;
  std::uint8_t current_color = 0;
  std::uint8_t background_color = 0;
  std::uint8_t current_color_mask = 0;
  std::uint8_t background_color_mask = 0;
};

struct HdtNode {
  int left = -1;
  int right = -1;
  int parent = -1;
  int atom_descriptor = -1;
  std::array<HdtPortKey, 2> output_ports{};
  int source_center_key = -1;
  int dp_offset = -1;
  int background_dp_offset = -1;
  int atom_count = 1;
  int height = 0;
  int live_owner = -1;
  HdtOperation operation = HdtOperation::Atom;
  HdtTransform post_transform = HdtTransform::None;
  HdtOutputKind output_kind = HdtOutputKind::I;
  std::uint8_t output_port_count = 0;
  std::uint8_t dp_dimension = 0;
  bool current_equals_background = false;
  std::uint8_t current_color_mask = 0;
  std::uint8_t background_color_mask = 0;
  std::uint8_t current_payload_class = 0xff;
  std::uint8_t background_payload_class = 0xff;
};

struct HdtBuildStats {
  int component_units = 0;
  int initial_atoms = 0;
  int reachable_nodes = 0;
  int max_recompute_degree = 0;
  int local_balance_violations = 0;
  int max_update_height = 0;
  int rounds = 0;
  long long gg_ops = 0;
  long long ig_ops = 0;
  long long cc_ops = 0;
  long long c_to_g_ops = 0;
  long long scalar_fmas = 0;
  std::size_t build_scratch_peak_bytes = 0;
  bool node_count_ok = false;
  bool atom_bound_ok = false;
};

struct HdtOwnerIndexView {
  const std::vector<int>* owner_position = nullptr;
  int position_begin = 0;
  int position_end = 0;
};

struct HdtPackedPayloadClass {
  std::vector<double> values;
  int free_head = -1;
};

struct HdtPackedPayloadArena {
  // Union of all reachable active-coordinate widths for K=2 and K=4.
  std::array<HdtPackedPayloadClass, 17> classes;

  void clear() {
    for (HdtPackedPayloadClass& payload_class : classes) {
      payload_class.values.clear();
      payload_class.free_head = -1;
    }
  }

  std::size_t capacity_bytes() const {
    std::size_t bytes = 0;
    for (const HdtPackedPayloadClass& payload_class : classes) {
      checked_hdt_accumulate_bytes(
          bytes, payload_class.values.capacity(), sizeof(double),
          "packed payload byte count exceeds size_t");
    }
    return bytes;
  }
};

struct Hdt {
  int palette = 0;
  int root = -1;
  int root_anchor_atom_node = -1;
  std::vector<HdtNode> nodes;
  HdtPackedPayloadArena current_payload;
  HdtPackedPayloadArena background_payload;
  std::vector<HdtExpandedAtom> atoms;
  std::vector<HdtPortKey> atom_ports;
  std::vector<int> owner_offsets;
  std::vector<int> owner_atom_nodes;
  std::vector<int> recompute_stack;
  int owner_position_begin = 0;
  bool owner_index_by_position = false;
  std::vector<std::uint8_t> dirty;
  HdtBuildStats stats;
  bool closed_root_required = true;
  bool background_alias_current = false;
  bool background_restricted = false;
};

struct HdtState {
  HdtState() = default;
  HdtState(const HdtState&) = delete;
  HdtState& operator=(const HdtState&) = delete;
  HdtState(HdtState&&) noexcept = default;
  HdtState& operator=(HdtState&&) noexcept = default;

  HdtComponentGraph components;
  Hdt hdt;
  long long recolor_changed_leaves = 0;
  long long recolor_marked_nodes = 0;
  long long transition_evals = 0;
  long long scalar_fmas = 0;
  bool initial_hdt_accounted = false;
};

Hdt build_hdt(const HdtComponentGraph& graph,
              bool require_closed_root = true,
              HdtOwnerIndexView owner_index = {});

void build_hdt_into(Hdt& output,
                    const HdtComponentGraph& graph,
                    bool require_closed_root = true,
                    HdtOwnerIndexView owner_index = {});

void build_frozen_fold_hdt_into(
    Hdt& output,
    const HdtComponentGraph& graph,
    bool require_closed_root = true);

HdtState build_hdt_state(HdtComponentGraph graph,
                         HdtOwnerIndexView owner_index = {});

bool validate_hdt(const HdtComponentGraph& graph,
                  const Hdt& hdt,
                  std::string* error = nullptr);

int hdt_recolor_owners(HdtState& state,
                       const std::vector<int>& owners,
                       std::uint8_t color);

int hdt_recolor_owner_range(HdtState& state,
                            const std::vector<int>& owner_order,
                            int begin,
                            int end,
                            std::uint8_t color);

int hdt_recolor_owner_ranges(
    HdtState& state,
    const std::vector<int>& owner_order,
    std::span<const std::pair<int, int>> owner_ranges,
    std::uint8_t color);

int hdt_recolor_owner(HdtState& state,
                      int owner,
                      std::uint8_t color);

std::array<double, 3> hdt_root_counts(const HdtState& state);

FullBdpDual hdt_node_dp(const Hdt& hdt, int node);
FullBdpVec hdt_node_current(const Hdt& hdt, int node);
FullBdpVec hdt_node_background(const Hdt& hdt, int node);

int hdt_live_atom_count(const HdtState& state);

}  // namespace wqfm::internal
