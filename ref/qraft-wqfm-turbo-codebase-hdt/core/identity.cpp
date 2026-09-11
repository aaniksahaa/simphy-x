#include "internal.hpp"

#include <sstream>

namespace wqfm::internal {
namespace {

constexpr int kMaterializeSupportTaxonLimit = 4096;

/// Returns whether an identity support bitset should be materialized.
bool should_materialize_support(int taxon_count) {
  return taxon_count >= 0 && taxon_count <= kMaterializeSupportTaxonLimit;
}

/// Materializes a support bitset for small taxon universes.
void materialize_support_if_small(Identity& identity) {
  if (!should_materialize_support(identity.universe_size)) {
    identity.support = Bitset();
    return;
  }
  identity.support = Bitset(identity.universe_size);
  for (int taxon : identity.support_taxa) {
    if (taxon >= 0 && taxon < identity.universe_size) {
      identity.support.set(taxon);
    }
  }
}

/// Assigns taxon ownership and weights for one identity.
void assign_nested_weights_for_identity(const Identity& identity,
                                        int owner,
                                        std::vector<int>& taxon_to_identity,
                                        std::vector<double>& taxon_weight) {
  if (identity.support_taxa.empty()) {
    return;
  }
  const double equal_weight =
      1.0 / static_cast<double>(identity.support_taxa.size());
  for (int taxon : identity.support_taxa) {
    if (taxon < 0 || taxon >= static_cast<int>(taxon_weight.size())) {
      throw std::invalid_argument("identity support taxon outside universe");
    }
    taxon_to_identity[taxon] = owner;
    taxon_weight[taxon] = equal_weight;
  }
}

}  // namespace

/// Creates one real singleton identity per taxon.
std::vector<Identity> make_initial_identities(const TaxonRegistry& taxa) {
  std::vector<Identity> identities;
  identities.reserve(taxa.size());
  for (int taxon = 0; taxon < taxa.size(); ++taxon) {
    Identity identity;
    identity.id = taxon;
    identity.universe_size = taxa.size();
    identity.support_taxa = {taxon};
    identity.stable_key = taxa.label_for(taxon);
    materialize_support_if_small(identity);
    identities.push_back(std::move(identity));
  }
  return identities;
}

/// Returns the universe size represented by one identity.
int identity_universe_size(const Identity& identity) {
  if (identity.universe_size > 0) {
    return identity.universe_size;
  }
  if (identity.support.size() > 0) {
    return identity.support.size();
  }
  int max_taxon = -1;
  for (int taxon : identity.support_taxa) {
    max_taxon = std::max(max_taxon, taxon);
  }
  return max_taxon + 1;
}

/// Returns the maximum universe size represented by identities.
int identity_universe_size(const std::vector<Identity>& identities) {
  int result = 0;
  for (const Identity& identity : identities) {
    result = std::max(result, identity_universe_size(identity));
  }
  return result;
}

/// Tests whether an identity contains a taxon.
bool identity_contains_taxon(const Identity& identity, int taxon) {
  if (taxon < 0) {
    return false;
  }
  if (taxon < identity.support.size()) {
    return identity.support.test(taxon);
  }
  return std::find(identity.support_taxa.begin(), identity.support_taxa.end(), taxon) !=
         identity.support_taxa.end();
}

/// Builds per-taxon ownership and mass for the current subproblem.
SubproblemWeights build_subproblem_weights(const std::vector<Identity>& identities,
                                           int taxon_count) {
  SubproblemWeights weights;
  weights.taxon_to_identity.assign(taxon_count, -1);
  weights.taxon_weight.assign(taxon_count, 0.0);
  std::vector<char> seen(taxon_count, false);
  for (int i = 0; i < static_cast<int>(identities.size()); ++i) {
    const Identity& identity = identities[i];
    for (int taxon : identity.support_taxa) {
      if (taxon < 0 || taxon >= taxon_count) {
        throw std::invalid_argument("identity support taxon outside universe");
      }
      if (seen[taxon]) {
        throw std::invalid_argument("current identity supports must form a partition");
      }
      seen[taxon] = true;
    }
    assign_nested_weights_for_identity(identity, i, weights.taxon_to_identity,
                                       weights.taxon_weight);
  }
  for (int taxon = 0; taxon < taxon_count; ++taxon) {
    if (seen[taxon] && weights.taxon_to_identity[taxon] < 0) {
      throw std::runtime_error("nested weight builder failed to assign an active taxon");
    }
  }
  return weights;
}

}  // namespace wqfm::internal
