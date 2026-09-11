#pragma once

#include "hat/hdt/hdt_state.hpp"

#include <cstdint>
#include <vector>

namespace wqfm::internal {

enum class HdtContractMode : std::uint8_t {
  LightCopy,
  HeavyReplace,
  InitialCompact,
};

struct HdtContractStats {
  HdtContractMode mode = HdtContractMode::LightCopy;
  int input_units = 0;
  int output_units = 0;
  int selected_live_atoms = 0;
  int anchors = 0;
  int corridors_emitted = 0;
  int background_units_absorbed = 0;
  int polytomy_arms_folded = 0;
  int size_bound = 0;
  bool size_bound_ok = false;
  long long fragment_hdt_nodes = 0;
  long long gg_ops = 0;
  long long ig_ops = 0;
  long long cc_ops = 0;
  long long c_to_g_ops = 0;
  long long mark_ns = 0;
  long long background_fold_ns = 0;
  long long corridor_fold_ns = 0;
  long long emit_ns = 0;
  long long hdt_rebuild_ns = 0;
  long long scalar_fmas = 0;
  std::size_t persistent_scratch_bytes = 0;
  std::size_t ephemeral_scratch_peak_bytes = 0;
  std::size_t scratch_peak_bytes = 0;
};

struct HdtContractResult {
  HdtState state;
  HdtContractStats stats;
};

struct HdtInitialStateResult {
  HdtState state;
  HdtContractStats contract;
  bool contracted = false;
};

struct HdtExtractStats {
  int source_hdt_nodes = 0;
  int marked_hdt_nodes = 0;
  int hanging_components = 0;
  int intermediate_units = 0;
  int selected_live_atoms = 0;
};

struct HdtExtractContractResult {
  HdtState state;
  HdtExtractStats extract;
  HdtContractStats contract;
};

struct HdtExtractWorkspace {
  std::uint64_t epoch = 0;
  std::vector<std::uint64_t> hdt_mark_epoch;
  std::vector<std::uint64_t> ordinary_epoch;
  std::vector<int> ordinary_unit;
  std::vector<HdtPortKey> ordinary_port;
  std::vector<std::uint64_t> center_epoch;
  std::vector<int> center_anchor;
  std::size_t peak_bytes = 0;
};

// A production selector has two independent views of the same owner set:
// direct iteration over represented owners, and O(1) membership through query
// Euler positions.  The owning dense form is retained for standalone callers
// and tests, but traversal code uses the borrowed interval form so constructing
// a heavy selector is O(1).
struct HdtLiveSelector {
  const std::vector<int>* owner_order = nullptr;
  const std::vector<int>* owner_position = nullptr;
  int owner_begin = 0;
  int owner_end = 0;
  int membership_begin = 0;
  int membership_end = 0;
  int exact_selected_atom_count = -1;
  std::vector<int> owned_owners;
  std::vector<std::uint8_t> dense_membership;

  int selected_owner_count() const;
  int selected_owner_at(int index) const;
  bool contains_owner(int owner) const;
  std::uint8_t operator[](std::size_t owner) const {
    return contains_owner(static_cast<int>(owner)) ? 1u : 0u;
  }
};

HdtContractResult hdt_contract(
    const HdtState& input,
    const HdtLiveSelector& selector,
    HdtContractMode mode = HdtContractMode::LightCopy);

HdtInitialStateResult hdt_build_initial_state(HdtComponentGraph graph);

HdtLiveSelector hdt_owner_selector(int owner_count,
                                   const std::vector<int>& owners);

HdtLiveSelector hdt_interval_owner_selector(
    const std::vector<int>& represented_owner_order,
    const std::vector<int>& owner_position,
    int owner_begin,
    int owner_end,
    int membership_begin,
    int membership_end,
    int exact_selected_atom_count);

int hdt_selected_atom_count(const HdtState& state,
                            const HdtLiveSelector& selector);

int hdt_contract_size_bound(int selected_live_atoms);

HdtExtractContractResult hdt_extract_contract(
    const HdtState& parent,
    const HdtLiveSelector& selector,
    HdtExtractWorkspace& workspace);

HdtExtractContractResult hdt_extract_contract(
    const HdtState& parent,
    const HdtLiveSelector& selector);

}  // namespace wqfm::internal
