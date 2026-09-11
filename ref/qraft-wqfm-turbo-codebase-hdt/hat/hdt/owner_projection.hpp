#pragma once

#include "hat/hat_internal.hpp"

#include <cstdint>
#include <utility>
#include <vector>

namespace wqfm::internal {

enum class OwnerKeyMode : std::uint8_t {
  Missing = 0,
  Dense = 1,
  Sorted = 2,
  Direct = 3,
};

struct OwnerKeyMap {
  OwnerKeyMode mode = OwnerKeyMode::Missing;
  std::vector<int> dense;
  std::vector<std::pair<int, int>> sorted;
  int direct_size = 0;

  int find(int key) const;
  int inferred_owner_count() const;
};

struct OwnerProjection {
  OwnerKeyMap taxon_keys;
  OwnerKeyMap identity_keys;
  std::vector<int> taxon_to_owner;
  std::vector<int> identity_to_owner;
  std::vector<std::uint8_t> live_owner;
  std::vector<int> live_taxa;
  std::vector<int> live_identities;
  int owner_count = 0;
};

OwnerProjection make_owner_projection(
    const SolverContext& ctx,
    const HatTaxonIndex& local_by_taxon,
    const std::vector<int>* source_identity_taxa);

int owner_for_taxon(const OwnerProjection& owners, int taxon);
int owner_for_identity(const OwnerProjection& owners, int identity);
int owner_for_source_leaf(int taxon,
                          int identity,
                          const OwnerProjection& owners);
int owner_for_gene_leaf(const GeneIndex::Node& node,
                        const OwnerProjection& owners);

}  // namespace wqfm::internal
