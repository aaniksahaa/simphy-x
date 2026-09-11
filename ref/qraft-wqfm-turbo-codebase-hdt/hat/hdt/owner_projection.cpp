#include "hat/hdt/owner_projection.hpp"

#include <algorithm>
#include <limits>

namespace wqfm::internal {

int OwnerKeyMap::find(int key) const {
  if (key < 0) return -1;
  switch (mode) {
    case OwnerKeyMode::Missing:
      return -1;
    case OwnerKeyMode::Dense:
      return key < static_cast<int>(dense.size())
                 ? dense[static_cast<std::size_t>(key)]
                 : -1;
    case OwnerKeyMode::Sorted: {
      const auto it = std::lower_bound(
          sorted.begin(), sorted.end(), key,
          [](const std::pair<int, int>& entry, int value) {
            return entry.first < value;
          });
      return it != sorted.end() && it->first == key ? it->second : -1;
    }
    case OwnerKeyMode::Direct:
      return key < direct_size ? key : -1;
  }
  return -1;
}

int OwnerKeyMap::inferred_owner_count() const {
  if (mode == OwnerKeyMode::Direct) return std::max(0, direct_size);
  int count = 0;
  const auto note_owner = [&](int owner) {
    if (owner >= 0) count = std::max(count, owner + 1);
  };
  if (mode == OwnerKeyMode::Dense) {
    for (int owner : dense) note_owner(owner);
  } else if (mode == OwnerKeyMode::Sorted) {
    for (const auto& entry : sorted) note_owner(entry.second);
  }
  return count;
}

int owner_for_taxon(const OwnerProjection& owners, int taxon) {
  if (owners.taxon_keys.mode != OwnerKeyMode::Missing) {
    return owners.taxon_keys.find(taxon);
  }
  return taxon >= 0 && taxon < static_cast<int>(owners.taxon_to_owner.size())
             ? owners.taxon_to_owner[static_cast<std::size_t>(taxon)]
             : -1;
}

int owner_for_identity(const OwnerProjection& owners, int identity) {
  if (owners.identity_keys.mode != OwnerKeyMode::Missing) {
    return owners.identity_keys.find(identity);
  }
  return identity >= 0 &&
                 identity < static_cast<int>(owners.identity_to_owner.size())
             ? owners.identity_to_owner[static_cast<std::size_t>(identity)]
             : -1;
}

int owner_for_source_leaf(int taxon,
                          int identity,
                          const OwnerProjection& owners) {
  const int owner = identity >= 0 ? owner_for_identity(owners, identity)
                                  : owner_for_taxon(owners, taxon);
  int owner_count = owners.owner_count;
  if (owner_count <= 0) {
    owner_count = std::max(owners.taxon_keys.inferred_owner_count(),
                           owners.identity_keys.inferred_owner_count());
    const auto note_owner = [&](int mapped) {
      if (mapped >= 0) owner_count = std::max(owner_count, mapped + 1);
    };
    for (int mapped : owners.taxon_to_owner) note_owner(mapped);
    for (int mapped : owners.identity_to_owner) note_owner(mapped);
  }
  if (owner < 0 || owner >= owner_count) return -1;
  if (!owners.live_owner.empty() &&
      (owner >= static_cast<int>(owners.live_owner.size()) ||
       owners.live_owner[static_cast<std::size_t>(owner)] == 0)) {
    return -1;
  }
  return owner;
}

int owner_for_gene_leaf(const GeneIndex::Node& node,
                        const OwnerProjection& owners) {
  return owner_for_source_leaf(node.taxon, node.identity, owners);
}

OwnerProjection make_owner_projection(
    const SolverContext& ctx,
    const HatTaxonIndex& local_by_taxon,
    const std::vector<int>* source_identity_taxa) {
  OwnerProjection owners;
  owners.owner_count = static_cast<int>(local_by_taxon.taxa.size());
  owners.live_owner.assign(static_cast<std::size_t>(owners.owner_count), 1);

  const auto linear_domain = [&](std::size_t domain) {
    return owners.owner_count > 0 &&
           domain <= static_cast<std::size_t>(owners.owner_count) * 2u;
  };
  std::size_t taxon_domain = context_taxa(ctx).size();
  if (!local_by_taxon.by_taxon.empty() &&
      local_by_taxon.by_taxon.back().first >= 0) {
    taxon_domain = std::max(
        taxon_domain,
        static_cast<std::size_t>(local_by_taxon.by_taxon.back().first) + 1u);
  }
  if (!local_by_taxon.by_taxon.empty()) {
    if (linear_domain(taxon_domain)) {
      owners.taxon_keys.mode = OwnerKeyMode::Dense;
      owners.taxon_keys.dense.assign(taxon_domain, -1);
      for (const auto& [taxon, owner] : local_by_taxon.by_taxon) {
        if (taxon >= 0 && taxon < static_cast<int>(taxon_domain)) {
          owners.taxon_keys.dense[static_cast<std::size_t>(taxon)] = owner;
        }
      }
    } else {
      owners.taxon_keys.mode = OwnerKeyMode::Sorted;
      owners.taxon_keys.sorted = local_by_taxon.by_taxon;
    }
  }

  auto build_identity_map = [&](std::size_t domain, const auto& taxon_at) {
    OwnerKeyMap& map = owners.identity_keys;
    if (linear_domain(domain)) {
      map.mode = OwnerKeyMode::Dense;
      map.dense.assign(domain, -1);
      for (int identity = 0; identity < static_cast<int>(domain); ++identity) {
        map.dense[static_cast<std::size_t>(identity)] =
            owner_for_taxon(owners, taxon_at(identity));
      }
      return;
    }
    map.mode = OwnerKeyMode::Sorted;
    map.sorted.reserve(std::min(domain, local_by_taxon.taxa.size()));
    for (int identity = 0; identity < static_cast<int>(domain); ++identity) {
      const int owner = owner_for_taxon(owners, taxon_at(identity));
      if (owner >= 0) map.sorted.push_back({identity, owner});
    }
  };

  if (source_identity_taxa != nullptr) {
    build_identity_map(source_identity_taxa->size(), [&](int identity) {
      return (*source_identity_taxa)[static_cast<std::size_t>(identity)];
    });
  } else if (ctx.active_gene_identity_ids != nullptr) {
    build_identity_map(ctx.active_gene_identity_ids->size(), [&](int identity) {
      const int64_t id =
          (*ctx.active_gene_identity_ids)[static_cast<std::size_t>(identity)];
      return id >= 0 && id <= std::numeric_limits<int>::max()
                 ? static_cast<int>(id)
                 : -1;
    });
  } else if (ctx.active_genes_match_identities) {
    owners.identity_keys.mode = OwnerKeyMode::Direct;
    owners.identity_keys.direct_size = owners.owner_count;
  }
  return owners;
}

}  // namespace wqfm::internal
