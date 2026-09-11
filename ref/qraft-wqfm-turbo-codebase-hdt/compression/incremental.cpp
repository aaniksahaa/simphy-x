#include "internal.hpp"

#include <limits>

namespace wqfm::internal {
namespace {

GeneIndex canonical_empty_gene_index() {
  GeneIndex result;
  result.root = -1;
  result.reference_root = -1;
  result.root_traversal = traversal(result, result.root);
  result.reference_traversal = result.root_traversal;
  return result;
}

/// Mutable tree adjacency with stable incidence order and O(1) edge unlink.
class MutableTreeAdjacency {
 public:
  explicit MutableTreeAdjacency(const std::vector<GeneIndex::Node>& nodes,
                                std::size_t extra_directed_edges = 0)
      : nodes_(nodes.size()) {
    constexpr std::size_t kMaxIndexedEdges =
        static_cast<std::size_t>(std::numeric_limits<int>::max());
    if (nodes.size() > kMaxIndexedEdges) {
      throw std::length_error("gene adjacency exceeds its integer node domain");
    }
    std::size_t directed_edges = 0;
    for (const GeneIndex::Node& node : nodes) {
      if (node.neighbors.size() > kMaxIndexedEdges - directed_edges) {
        throw std::length_error("gene adjacency exceeds its half-edge domain");
      }
      directed_edges += node.neighbors.size();
    }
    if (extra_directed_edges > kMaxIndexedEdges - directed_edges) {
      throw std::length_error("gene adjacency exceeds its half-edge domain");
    }
    edges_.reserve(directed_edges + extra_directed_edges);
    for (int node = 0; node < static_cast<int>(nodes.size()); ++node) {
      for (int neighbor : nodes[static_cast<std::size_t>(node)].neighbors) {
        append_directed(node, neighbor);
      }
    }
    pair_original_edges();
  }

  void add_node() { nodes_.push_back(NodeLinks{}); }

  int degree(int node) const {
    return nodes_[static_cast<std::size_t>(node)].degree;
  }

  int first_edge(int node) const {
    return nodes_[static_cast<std::size_t>(node)].first;
  }

  int next_edge(int edge) const {
    return edges_[static_cast<std::size_t>(edge)].next;
  }

  int from(int edge) const {
    return edges_[static_cast<std::size_t>(edge)].from;
  }

  int neighbor(int edge) const {
    return edges_[static_cast<std::size_t>(edge)].to;
  }

  template <typename Fn>
  void for_each_edge(int node, Fn&& fn) const {
    for (int edge = first_edge(node); edge >= 0; edge = next_edge(edge)) {
      fn(edge);
    }
  }

  int add_edge(int a, int b) {
    const int ab = append_directed(a, b);
    const int ba = append_directed(b, a);
    edges_[static_cast<std::size_t>(ab)].reverse = ba;
    edges_[static_cast<std::size_t>(ba)].reverse = ab;
    return ab;
  }

  void remove_edge(int edge) {
    HalfEdge& forward = edges_[static_cast<std::size_t>(edge)];
    if (!forward.active) return;
    const int reverse = forward.reverse;
    if (reverse < 0 ||
        reverse >= static_cast<int>(edges_.size()) ||
        !edges_[static_cast<std::size_t>(reverse)].active) {
      throw std::invalid_argument("gene adjacency must be a symmetric tree");
    }
    unlink_directed(edge);
    unlink_directed(reverse);
  }

 private:
  struct NodeLinks {
    int first = -1;
    int last = -1;
    int degree = 0;
  };

  struct HalfEdge {
    int from = -1;
    int to = -1;
    int reverse = -1;
    int previous = -1;
    int next = -1;
    bool active = true;
  };

  int append_directed(int from, int to) {
    if (edges_.size() >=
        static_cast<std::size_t>(std::numeric_limits<int>::max())) {
      throw std::length_error("gene adjacency exceeds its half-edge domain");
    }
    const int edge = static_cast<int>(edges_.size());
    NodeLinks& links = nodes_[static_cast<std::size_t>(from)];
    edges_.push_back(HalfEdge{from, to, -1, links.last, -1, true});
    if (links.last >= 0) {
      edges_[static_cast<std::size_t>(links.last)].next = edge;
    } else {
      links.first = edge;
    }
    links.last = edge;
    ++links.degree;
    return edge;
  }

  void unlink_directed(int edge) {
    HalfEdge& current = edges_[static_cast<std::size_t>(edge)];
    NodeLinks& links = nodes_[static_cast<std::size_t>(current.from)];
    if (current.previous >= 0) {
      edges_[static_cast<std::size_t>(current.previous)].next = current.next;
    } else {
      links.first = current.next;
    }
    if (current.next >= 0) {
      edges_[static_cast<std::size_t>(current.next)].previous = current.previous;
    } else {
      links.last = current.previous;
    }
    current.active = false;
    --links.degree;
  }

  void pair_original_edges() {
    const int n = static_cast<int>(nodes_.size());
    std::vector<unsigned char> discovered(static_cast<std::size_t>(n), 0);
    std::vector<int> stack;
    stack.reserve(static_cast<std::size_t>(n));
    int component_count = 0;
    for (int start = 0; start < n; ++start) {
      if (discovered[static_cast<std::size_t>(start)] != 0) continue;
      ++component_count;
      discovered[static_cast<std::size_t>(start)] = 1;
      stack.push_back(start);
      while (!stack.empty()) {
        const int node = stack.back();
        stack.pop_back();
        for (int edge = first_edge(node); edge >= 0; edge = next_edge(edge)) {
          HalfEdge& forward = edges_[static_cast<std::size_t>(edge)];
          if (forward.reverse >= 0) continue;
          const int child = forward.to;
          if (child < 0 || child >= n) {
            throw std::invalid_argument("gene adjacency contains an invalid node");
          }
          if (discovered[static_cast<std::size_t>(child)] != 0) continue;

          int reverse = -1;
          for (int candidate = first_edge(child); candidate >= 0;
               candidate = next_edge(candidate)) {
            const HalfEdge& back = edges_[static_cast<std::size_t>(candidate)];
            if (back.to == node && back.reverse < 0) {
              reverse = candidate;
              break;
            }
          }
          if (reverse < 0) {
            throw std::invalid_argument("gene adjacency must be symmetric");
          }
          forward.reverse = reverse;
          edges_[static_cast<std::size_t>(reverse)].reverse = edge;
          discovered[static_cast<std::size_t>(child)] = 1;
          stack.push_back(child);
        }
      }
    }
    for (const HalfEdge& edge : edges_) {
      if (edge.reverse < 0) {
        throw std::invalid_argument("gene adjacency must be acyclic and symmetric");
      }
    }
    if (component_count > 1) {
      throw std::invalid_argument("gene adjacency must be connected");
    }
  }

  std::vector<NodeLinks> nodes_;
  std::vector<HalfEdge> edges_;
};

}  // namespace

/// Removes unlabeled degree-two gene nodes while preserving leaf taxon mapping.
void suppress_gene_degree_two(GeneIndex& gene, int taxon_count) {
  if (gene.nodes.size() >
      static_cast<std::size_t>(std::numeric_limits<int>::max())) {
    throw std::length_error("gene index exceeds its integer node domain");
  }
  const int n = static_cast<int>(gene.nodes.size());
  if (n == 0) return;
  MutableTreeAdjacency adjacency(
      gene.nodes, static_cast<std::size_t>(n) * 2);
  std::vector<unsigned char> active(static_cast<std::size_t>(n), 1);

  for (int i = 0; i < n; ++i) {
    if (active[static_cast<std::size_t>(i)] == 0 ||
        gene.nodes[static_cast<std::size_t>(i)].taxon >= 0 ||
        gene.nodes[static_cast<std::size_t>(i)].identity >= 0 ||
        adjacency.degree(i) != 2) {
      continue;
    }
    const int first = adjacency.first_edge(i);
    const int second = adjacency.next_edge(first);
    const int a = adjacency.neighbor(first);
    const int b = adjacency.neighbor(second);
    adjacency.remove_edge(first);
    adjacency.remove_edge(second);
    active[static_cast<std::size_t>(i)] = 0;
    adjacency.add_edge(a, b);
  }

  std::vector<int> remap(n, -1);
  int next = 0;
  for (int i = 0; i < n; ++i) {
    if (active[static_cast<std::size_t>(i)] != 0) remap[i] = next++;
  }
  if (next == n) return;

  std::vector<GeneIndex::Node> compact(next);
  for (int i = 0; i < n; ++i) {
    const int ni = remap[i];
    if (ni < 0) continue;
    compact[ni] = gene.nodes[i];
    compact[ni].neighbors.clear();
    adjacency.for_each_edge(i, [&](int edge) {
      const int nb = adjacency.neighbor(edge);
      if (nb >= 0 && nb < n && remap[nb] >= 0) {
        compact[ni].neighbors.push_back(remap[nb]);
      }
    });
  }
  gene.nodes = std::move(compact);
  gene.taxon_to_node.assign(taxon_count, -1);
  for (int i = 0; i < static_cast<int>(gene.nodes.size()); ++i) {
    const int taxon = gene.nodes[i].taxon;
    if (taxon >= 0 && taxon < taxon_count) gene.taxon_to_node[taxon] = i;
  }
  gene.root = -1;
  for (int i = 0; i < static_cast<int>(gene.nodes.size()); ++i) {
    if (gene.nodes[i].taxon < 0 && gene.nodes[i].identity < 0 &&
        gene.nodes[i].neighbors.size() >= 2) {
      gene.root = i;
      break;
    }
  }
  if (gene.root < 0 && !gene.nodes.empty()) gene.root = 0;
}

/// Converts a public Newick tree into taxon-indexed internal storage.
GeneIndex index_gene(const GeneTree& gene, const TaxonRegistry& taxa) {
  if (gene.nodes.size() >
      static_cast<std::size_t>(std::numeric_limits<int>::max())) {
    throw std::length_error("gene tree exceeds its integer node domain");
  }
  if (gene.nodes.empty()) return canonical_empty_gene_index();
  if (gene.root < 0 ||
      gene.root >= static_cast<int>(gene.nodes.size())) {
    throw std::invalid_argument("gene tree has an invalid root");
  }
  for (int i = 0; i < static_cast<int>(gene.nodes.size()); ++i) {
    const Tree::Node& node = gene.nodes[static_cast<std::size_t>(i)];
    const bool labeled = !node.label.empty();
    const std::size_t degree = node.neighbors.size();
    if (labeled && degree > 1u) {
      throw std::invalid_argument("gene tree has a labeled internal node");
    }
    // Newick permits a unary, unlabeled parse root (for example `(A);`).
    // Every other degree-zero/one vertex must be a physical labeled leaf.
    if (!labeled && degree <= 1u && !(i == gene.root && degree == 1u)) {
      throw std::invalid_argument("gene tree has an unlabeled leaf");
    }
  }
  GeneIndex out;
  out.nodes.resize(gene.nodes.size());
  out.root = gene.root;
  out.taxon_to_node.assign(taxa.size(), -1);
  for (int i = 0; i < static_cast<int>(gene.nodes.size()); ++i) {
    out.nodes[i].neighbors = gene.nodes[i].neighbors;
    if (!gene.nodes[i].label.empty()) {
      const int taxon = taxa.id_for(gene.nodes[i].label);
      out.nodes[i].taxon = taxon;
      out.taxon_to_node[taxon] = i;
    }
  }
  suppress_gene_degree_two(out, taxa.size());
  out.root_traversal = traversal(out, out.root);
  // Root every nonempty taxon tree at its smallest present taxon. This makes
  // the represented side of an unrooted split independent of the arbitrary
  // Newick parse root, including when global taxon zero is absent.
  out.reference_root = out.root;
  for (int taxon = 0; taxon < static_cast<int>(out.taxon_to_node.size());
       ++taxon) {
    if (out.taxon_to_node[static_cast<std::size_t>(taxon)] >= 0) {
      out.reference_root =
          out.taxon_to_node[static_cast<std::size_t>(taxon)];
      break;
    }
  }
  out.reference_traversal = out.reference_root == out.root
                                ? out.root_traversal
                                : traversal(out, out.reference_root);
  return out;
}

GeneSplitLeafProjection make_gene_split_leaf_projection(
    const GeneIndex& gene,
    int taxon_count,
    const std::vector<std::array<uint64_t, 4>>& taxon_fingerprints) {
  if (taxon_count < 0 ||
      taxon_fingerprints.size() < static_cast<std::size_t>(taxon_count)) {
    throw std::invalid_argument("invalid taxon fingerprint domain");
  }
  if (gene.nodes.size() >
      static_cast<std::size_t>(std::numeric_limits<int>::max())) {
    throw std::length_error("gene split projection exceeds its integer domain");
  }

  GeneSplitLeafProjection out;
  const int node_count = static_cast<int>(gene.nodes.size());
  out.unique_size.assign(static_cast<std::size_t>(node_count), 0);
  out.unique_fingerprint.assign(static_cast<std::size_t>(node_count), {});
  out.duplicate_balance.assign(static_cast<std::size_t>(node_count), {});
  out.reference_marker.assign(static_cast<std::size_t>(node_count), 0);
  if (node_count == 0) return out;

  std::vector<int> occurrence_count(static_cast<std::size_t>(taxon_count), 0);
  std::vector<int> representative(static_cast<std::size_t>(taxon_count), -1);
  for (int node = 0; node < node_count; ++node) {
    const int taxon = gene.nodes[static_cast<std::size_t>(node)].taxon;
    if (taxon < 0) continue;
    if (taxon >= taxon_count) {
      throw std::invalid_argument("gene taxon lies outside the split universe");
    }
    ++occurrence_count[static_cast<std::size_t>(taxon)];
    representative[static_cast<std::size_t>(taxon)] = node;
  }

  for (int taxon = 0; taxon < taxon_count; ++taxon) {
    if (occurrence_count[static_cast<std::size_t>(taxon)] == 0) continue;
    if (out.reference_taxon < 0) out.reference_taxon = taxon;
    ++out.distinct_taxa;
  }
  if (out.reference_taxon < 0) return out;

  const int canonical_reference =
      representative[static_cast<std::size_t>(out.reference_taxon)];
  if (gene.reference_root != canonical_reference ||
      gene.reference_traversal.preorder.size() != gene.nodes.size()) {
    throw std::logic_error("gene split traversal is not canonically rooted");
  }

  for (int node = 0; node < node_count; ++node) {
    const int taxon = gene.nodes[static_cast<std::size_t>(node)].taxon;
    if (taxon < 0) continue;
    const int copies = occurrence_count[static_cast<std::size_t>(taxon)];
    const bool is_representative =
        representative[static_cast<std::size_t>(taxon)] == node;
    if (is_representative) {
      out.unique_size[static_cast<std::size_t>(node)] = 1;
      out.unique_fingerprint[static_cast<std::size_t>(node)] =
          taxon_fingerprints[static_cast<std::size_t>(taxon)];
      if (taxon == out.reference_taxon) {
        out.reference_marker[static_cast<std::size_t>(node)] = 1;
      }
    }
    if (copies <= 1) continue;
    for (int word = 0; word < 4; ++word) {
      const uint64_t fingerprint =
          taxon_fingerprints[static_cast<std::size_t>(taxon)]
                            [static_cast<std::size_t>(word)];
      out.duplicate_balance[static_cast<std::size_t>(node)]
                           [static_cast<std::size_t>(word)] =
          is_representative
              ? uint64_t{0} -
                    static_cast<uint64_t>(copies - 1) * fingerprint
              : fingerprint;
    }
  }
  return out;
}

/// Builds rooted traversal arrays for a gene index.
GeneTraversal traversal(const GeneIndex& gene, int root) {
  GeneTraversal tr;
  const int n = static_cast<int>(gene.nodes.size());
  tr.parent.assign(n, -1);
  tr.tin.assign(n, 0);
  tr.tout.assign(n, 0);
  if (n == 0) return tr;
  root = root >= 0 ? root : gene.root;
  std::vector<std::pair<int, int>> stack{{root, 0}};
  tr.parent[root] = root;
  while (!stack.empty()) {
    const int node = stack.back().first;
    int& next = stack.back().second;
    if (next == 0) {
      tr.tin[node] = static_cast<int>(tr.preorder.size());
      tr.preorder.push_back(node);
    }
    if (next < static_cast<int>(gene.nodes[node].neighbors.size())) {
      const int nb = gene.nodes[node].neighbors[next++];
      if (nb == tr.parent[node]) continue;
      tr.parent[nb] = node;
      stack.push_back({nb, 0});
      continue;
    }
    tr.tout[node] = static_cast<int>(tr.preorder.size());
    tr.postorder.push_back(node);
    stack.pop_back();
  }
  return tr;
}

/// Returns the identity owner represented at a gene node.
int leaf_owner(const GeneIndex& gene,
               const SubproblemWeights& weights,
               int node) {
  const auto& current = gene.nodes[node];
  if (current.identity >= 0) {
    return current.identity;
  }
  const int taxon = current.taxon;
  if (taxon >= 0 && taxon < static_cast<int>(weights.taxon_to_identity.size())) {
    return weights.taxon_to_identity[taxon];
  }
  return -1;
}

/// Returns the identity mass represented at a gene node.
double leaf_mass(const GeneIndex& gene,
                 const SubproblemWeights& weights,
                 int node) {
  const auto& current = gene.nodes[node];
  if (current.identity >= 0) {
    return current.mass;
  }
  const int taxon = current.taxon;
  if (taxon >= 0 && taxon < static_cast<int>(weights.taxon_weight.size())) {
    return weights.taxon_weight[taxon];
  }
  return 0.0;
}

/// Returns the real-identity count represented at a gene node.
double leaf_real_count(const GeneIndex& gene,
                       const std::vector<Identity>& identities,
                       const SubproblemWeights& weights,
                       int node) {
  const auto& current = gene.nodes[node];
  if (current.identity >= 0) {
    return current.real_count;
  }
  const int owner = leaf_owner(gene, weights, node);
  return owner >= 0 && owner < static_cast<int>(identities.size()) ? 1.0
                                                                   : 0.0;
}

/// Builds the exact induced tree for a caller-supplied physical-leaf map.
template <typename InputValidator, typename LeafMapper>
GeneIndex compress_gene_with_leaf_mapper(const GeneIndex& original,
                                         int identity_count,
                                         InputValidator&& validate_input,
                                         LeafMapper&& map_leaf) {
  (void)identity_count;
  /// Mutable node used while pruning inactive leaves and suppressing corridors.
  struct MNode {
    int identity = -1;
    double mass = 0.0;
    double real_count = 0.0;
    bool is_leaf = false;
    bool deleted = false;
    bool queued = false;
  };

  if (original.nodes.size() >
      static_cast<std::size_t>(std::numeric_limits<int>::max())) {
    throw std::length_error("gene compression exceeds its integer node domain");
  }
  const int n = static_cast<int>(original.nodes.size());
  MutableTreeAdjacency adjacency(
      original.nodes, static_cast<std::size_t>(n) * 2 + 8);
  validate_input();
  if (n == 0) return canonical_empty_gene_index();
  std::vector<MNode> nodes(n);
  int mapped_leaf_count = 0;
  for (int i = 0; i < n; ++i) {
    int identity = -1;
    double mass = 0.0;
    double real_count = 0.0;
    if (map_leaf(original.nodes[i], identity, mass, real_count)) {
      nodes[i].is_leaf = true;
      nodes[i].identity = identity;
      nodes[i].mass = mass;
      nodes[i].real_count = real_count;
      ++mapped_leaf_count;
    }
  }
  if (mapped_leaf_count == 0) return canonical_empty_gene_index();

  std::vector<int> queue;
  queue.reserve(nodes.size());
  auto schedule = [&](int node) {
    if (node < 0 || node >= static_cast<int>(nodes.size())) return;
    MNode& current = nodes[static_cast<std::size_t>(node)];
    if (current.deleted || current.is_leaf || current.queued) return;
    current.queued = true;
    queue.push_back(node);
  };

  for (int i = 0; i < static_cast<int>(nodes.size()); ++i) schedule(i);
  for (std::size_t qpos = 0; qpos < queue.size(); ++qpos) {
    const int node = queue[qpos];
    MNode& current = nodes[static_cast<std::size_t>(node)];
    current.queued = false;
    if (current.deleted || current.is_leaf) continue;

    const int degree = adjacency.degree(node);
    if (degree == 2) {
      const int first = adjacency.first_edge(node);
      const int second = adjacency.next_edge(first);
      const int a = adjacency.neighbor(first);
      const int b = adjacency.neighbor(second);
      adjacency.remove_edge(first);
      adjacency.remove_edge(second);
      current.deleted = true;
      adjacency.add_edge(a, b);
    } else if (degree <= 1) {
      int neighbor = -1;
      if (degree == 1) {
        const int edge = adjacency.first_edge(node);
        neighbor = adjacency.neighbor(edge);
        adjacency.remove_edge(edge);
      }
      current.deleted = true;
      schedule(neighbor);
    }
  }

  std::vector<int> old_to_new(nodes.size(), -1);
  int active = 0;
  for (int i = 0; i < static_cast<int>(nodes.size()); ++i) {
    if (!nodes[i].deleted) old_to_new[i] = active++;
  }
  if (active == 0) return canonical_empty_gene_index();

  GeneIndex result;
  result.nodes.resize(active);
  result.root = 0;
  for (int i = 0; i < static_cast<int>(nodes.size()); ++i) {
    const int ni = old_to_new[i];
    if (ni < 0) continue;
    adjacency.for_each_edge(i, [&](int edge) {
      const int nb = adjacency.neighbor(edge);
      const int nn = nb >= 0 && nb < static_cast<int>(old_to_new.size())
                         ? old_to_new[nb]
                         : -1;
      if (nn >= 0) result.nodes[ni].neighbors.push_back(nn);
    });
    if (nodes[i].is_leaf) {
      result.nodes[ni].identity = nodes[i].identity;
      result.nodes[ni].mass = nodes[i].mass;
      result.nodes[ni].real_count = nodes[i].real_count;
    }
  }
  for (int i = 0; i < static_cast<int>(result.nodes.size()); ++i) {
    if (result.nodes[i].identity < 0 && result.nodes[i].neighbors.size() >= 2) {
      result.root = i;
      break;
    }
  }
  result.reference_root = result.root;
  result.root_traversal = traversal(result, result.root);
  result.reference_traversal = result.root_traversal;
  return result;
}

/// Prunes inactive leaves for the current subproblem without coalescing owners.
GeneIndex compress_gene_for_subproblem(const GeneIndex& original,
                                       const std::vector<Identity>& identities,
                                       const SubproblemWeights& weights) {
  return compress_gene_with_leaf_mapper(
      original, static_cast<int>(identities.size()),
      [&]() {
        if (weights.taxon_to_identity.size() != weights.taxon_weight.size()) {
          throw std::invalid_argument(
              "subproblem weight vectors must have equal length");
        }
      },
      [&](const GeneIndex::Node& node, int& identity, double& mass,
          double& real_count) {
        const int taxon = node.taxon;
        if (taxon < 0) return false;
        if (taxon >= static_cast<int>(weights.taxon_to_identity.size())) {
          return false;
        }
        identity = weights.taxon_to_identity[static_cast<std::size_t>(taxon)];
        if (identity < 0 || identity >= static_cast<int>(identities.size())) {
          return false;
        }
        mass = weights.taxon_weight[static_cast<std::size_t>(taxon)];
        real_count = 1.0;
        return true;
      });
}

/// Prunes a previously reduced gene under an identity remapping.
GeneIndex compress_gene_for_identity_remap(
    const GeneIndex& original,
    const std::vector<Identity>& identities,
    const std::vector<int>& old_identity_to_new,
    const std::vector<double>& old_identity_scale) {
  return compress_gene_with_leaf_mapper(
      original, static_cast<int>(identities.size()),
      []() {},
      [&](const GeneIndex::Node& node, int& identity, double& mass,
          double& real_count) {
        const int old_identity = node.identity;
        if (old_identity < 0 ||
            old_identity >= static_cast<int>(old_identity_to_new.size())) {
          return false;
        }
        identity = old_identity_to_new[old_identity];
        if (identity < 0 || identity >= static_cast<int>(identities.size())) {
          return false;
        }
        const double scale =
            old_identity < static_cast<int>(old_identity_scale.size())
                ? old_identity_scale[old_identity]
                : 1.0;
        mass = node.mass * scale;
        real_count = node.real_count;
        return true;
      });
}

}  // namespace wqfm::internal
