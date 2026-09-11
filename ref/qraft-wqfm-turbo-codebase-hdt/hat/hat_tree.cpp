#include "hat/hat_internal.hpp"

#include <algorithm>
#include <stdexcept>

namespace wqfm::internal {

HatTaxonIndex make_hat_taxon_index(const std::vector<int>& taxa,
                                   int taxon_count) {
  HatTaxonIndex index;
  std::vector<std::pair<int, int>> first_positions;
  first_positions.reserve(taxa.size());
  for (int position = 0; position < static_cast<int>(taxa.size()); ++position) {
    const int taxon = taxa[static_cast<std::size_t>(position)];
    if (taxon >= 0 && taxon < taxon_count) {
      first_positions.push_back({taxon, position});
    }
  }
  std::sort(first_positions.begin(), first_positions.end());
  first_positions.erase(
      std::unique(first_positions.begin(), first_positions.end(),
                  [](const auto& a, const auto& b) {
                    return a.first == b.first;
                  }),
      first_positions.end());
  std::sort(first_positions.begin(), first_positions.end(),
            [](const auto& a, const auto& b) {
              return a.second < b.second;
            });

  index.taxa.reserve(first_positions.size());
  index.by_taxon.reserve(first_positions.size());
  for (const auto& [taxon, unused_position] : first_positions) {
    (void)unused_position;
    const int local = static_cast<int>(index.taxa.size());
    index.taxa.push_back(taxon);
    index.by_taxon.push_back({taxon, local});
  }
  std::sort(index.by_taxon.begin(), index.by_taxon.end());
  benchmark_counter_add(benchmark_stats().hat_active_local_taxon_entries,
                        static_cast<long long>(taxa.size()));
  benchmark_counter_add(benchmark_stats().hat_active_local_index_entries,
                        static_cast<long long>(index.by_taxon.size()));
  return index;
}

int hat_active_index(const HatActiveSet& active, int taxon) {
  return active.active_by_taxon.find(taxon);
}

HatActiveSet make_hat_active_set(const std::vector<Identity>& identities,
                                 int taxon_count) {
  HatActiveSet active;
  std::vector<int> requested_taxa;
  for (const Identity& identity : identities) {
    for (int taxon : identity.support_taxa) {
      if (taxon >= 0 && taxon < taxon_count) requested_taxa.push_back(taxon);
    }
  }
  active.active_by_taxon =
      make_hat_taxon_index(requested_taxa, taxon_count);
  active.taxa = active.active_by_taxon.taxa;
  return active;
}

void note_hat_source_bank_scoring_use(const std::vector<GeneIndex>& genes) {
  long long source_nodes = 0;
  for (const GeneIndex& gene : genes) {
    source_nodes += static_cast<long long>(gene.nodes.size());
  }
  // Keep the established telemetry keys stable while recording immutable-bank
  // scoring use instead of removed localized-bank construction.
  benchmark_counter_add(benchmark_stats().hat_localized_scoring_calls);
  benchmark_counter_add(benchmark_stats().hat_localized_scoring_gene_nodes,
                        source_nodes);
}

std::vector<int> hat_taxa_from_identities(
    const std::vector<Identity>& identities) {
  std::vector<int> taxa;
  for (const Identity& identity : identities) {
    taxa.insert(taxa.end(), identity.support_taxa.begin(),
                identity.support_taxa.end());
  }
  std::sort(taxa.begin(), taxa.end());
  taxa.erase(std::unique(taxa.begin(), taxa.end()), taxa.end());
  return taxa;
}

SolverTree hat_base_tree(const std::vector<Identity>& identities) {
  SolverTree tree;
  if (identities.empty()) return tree;
  if (identities.size() == 1) {
    tree.nodes.push_back(SolverTree::Node{{}, identities[0].id, true});
    return tree;
  }
  if (identities.size() == 2) {
    tree.nodes.push_back(SolverTree::Node{{1}, identities[0].id, true});
    tree.nodes.push_back(SolverTree::Node{{0}, identities[1].id, true});
    return tree;
  }
  tree.nodes.push_back(SolverTree::Node{{}, -1, true});
  for (const Identity& identity : identities) {
    const int leaf = static_cast<int>(tree.nodes.size());
    tree.nodes.push_back(SolverTree::Node{{0}, identity.id, true});
    tree.nodes[0].neighbors.push_back(leaf);
  }
  return tree;
}

std::vector<int> hat_active_neighbors(const SolverTree& tree, int node) {
  std::vector<int> out;
  if (node < 0 || node >= static_cast<int>(tree.nodes.size()) ||
      !tree.nodes[static_cast<std::size_t>(node)].active) {
    return out;
  }
  for (int nb : tree.nodes[static_cast<std::size_t>(node)].neighbors) {
    if (nb >= 0 && nb < static_cast<int>(tree.nodes.size()) &&
        tree.nodes[static_cast<std::size_t>(nb)].active) {
      out.push_back(nb);
    }
  }
  return out;
}

bool hat_has_edge(const SolverTree& tree, int a, int b) {
  if (a < 0 || b < 0 || a >= static_cast<int>(tree.nodes.size()) ||
      b >= static_cast<int>(tree.nodes.size())) {
    return false;
  }
  const auto& nbs = tree.nodes[static_cast<std::size_t>(a)].neighbors;
  return tree.nodes[static_cast<std::size_t>(a)].active &&
         tree.nodes[static_cast<std::size_t>(b)].active &&
         std::find(nbs.begin(), nbs.end(), b) != nbs.end();
}

void hat_add_edge(SolverTree& tree, int a, int b) {
  if (a < 0 || b < 0 || a >= static_cast<int>(tree.nodes.size()) ||
      b >= static_cast<int>(tree.nodes.size())) {
    return;
  }
  auto& an = tree.nodes[static_cast<std::size_t>(a)].neighbors;
  if (std::find(an.begin(), an.end(), b) == an.end()) an.push_back(b);
  auto& bn = tree.nodes[static_cast<std::size_t>(b)].neighbors;
  if (std::find(bn.begin(), bn.end(), a) == bn.end()) bn.push_back(a);
}

void hat_remove_edge(SolverTree& tree, int a, int b) {
  if (a < 0 || b < 0 || a >= static_cast<int>(tree.nodes.size()) ||
      b >= static_cast<int>(tree.nodes.size())) {
    return;
  }
  auto& an = tree.nodes[static_cast<std::size_t>(a)].neighbors;
  an.erase(std::remove(an.begin(), an.end(), b), an.end());
  auto& bn = tree.nodes[static_cast<std::size_t>(b)].neighbors;
  bn.erase(std::remove(bn.begin(), bn.end(), a), bn.end());
}

namespace {

void validate_attachment_edge(const SolverTree& tree,
                              const AttachmentEdge& edge) {
  if (edge.a >= 0 && edge.b >= 0) {
    if (hat_has_edge(tree, edge.a, edge.b) &&
        hat_has_edge(tree, edge.b, edge.a)) {
      return;
    }
    throw std::logic_error("invalid HAT attachment edge");
  }

  int active_count = 0;
  int singleton = -1;
  for (int a = 0; a < static_cast<int>(tree.nodes.size()); ++a) {
    if (tree.nodes[static_cast<std::size_t>(a)].active) {
      ++active_count;
      singleton = a;
    }
  }

  if (active_count == 0) {
    if (edge.a == -1 && edge.b == -1) return;
    throw std::logic_error("nonempty attachment endpoint for an empty HAT tree");
  }
  if (active_count == 1) {
    if (edge.a == singleton && edge.b == -1) return;
    throw std::logic_error("invalid singleton HAT attachment endpoint");
  }
  throw std::logic_error("invalid HAT attachment edge sentinel");
}

int subdivide_attachment_edge(SolverTree& tree, const AttachmentEdge& edge) {
  if (edge.a < 0) {
    if (edge.b != -1) {
      throw std::logic_error("invalid empty-tree HAT attachment endpoint");
    }
    const int node = static_cast<int>(tree.nodes.size());
    tree.nodes.push_back(SolverTree::Node{{}, -1, true});
    return node;
  }
  if (edge.b == -1) {
    // A singleton's sole vertex is already its attachment point. Introducing
    // a connector here leaves an active unlabeled degree-two vertex. Besides
    // being non-binary, that vertex is not a one-X median, so later range
    // accumulation assigns the two incident edges incomplete scores and can
    // make a spurious attachment win during a singleton-completion chain.
    return edge.a;
  }
  if (!hat_has_edge(tree, edge.a, edge.b)) {
    throw std::logic_error("HAT attachment edge disappeared during glue");
  }
  hat_remove_edge(tree, edge.a, edge.b);
  const int node = static_cast<int>(tree.nodes.size());
  tree.nodes.push_back(SolverTree::Node{{}, -1, true});
  hat_add_edge(tree, node, edge.a);
  hat_add_edge(tree, node, edge.b);
  return node;
}

}  // namespace

SolverTree hat_glue(SolverTree left,
                    SolverTree right,
                    AttachmentEdge left_edge,
                    AttachmentEdge right_edge) {
  validate_attachment_edge(left, left_edge);
  validate_attachment_edge(right, right_edge);
  const int right_offset = static_cast<int>(left.nodes.size());
  for (const SolverTree::Node& node : right.nodes) {
    SolverTree::Node copy;
    copy.identity_id = node.identity_id;
    copy.active = node.active;
    copy.neighbors.reserve(node.neighbors.size());
    for (int nb : node.neighbors) copy.neighbors.push_back(nb + right_offset);
    left.nodes.push_back(std::move(copy));
  }
  if (right_edge.a >= 0) right_edge.a += right_offset;
  if (right_edge.b >= 0) right_edge.b += right_offset;
  const int left_connector = subdivide_attachment_edge(left, left_edge);
  const int right_connector = subdivide_attachment_edge(left, right_edge);
  hat_add_edge(left, left_connector, right_connector);
  return left;
}

}  // namespace wqfm::internal
