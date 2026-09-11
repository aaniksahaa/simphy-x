#include "internal.hpp"

#include <cctype>
#include <sstream>

namespace wqfm {

/// Returns every non-empty node label in the tree.
std::vector<std::string> Tree::leaf_labels() const {
  std::vector<std::string> labels;
  for (const auto &node : nodes) {
    if (!node.label.empty()) {
      labels.push_back(node.label);
    }
  }
  return labels;
}

namespace {

/// Escapes one label for Newick output.
std::string escape_newick_label(const std::string &label) {
  bool needs_quote = label.empty();
  for (char ch : label) {
    if (std::isspace(static_cast<unsigned char>(ch)) || ch == '(' ||
        ch == ')' || ch == ',' || ch == ':' || ch == ';' || ch == '\'' ||
        ch == '[' || ch == ']') {
      needs_quote = true;
      break;
    }
  }
  if (!needs_quote) {
    return label;
  }
  std::string out = "'";
  for (char ch : label) {
    if (ch == '\'') {
      out += "''";
    } else {
      out += ch;
    }
  }
  out += "'";
  return out;
}

/// Serializes a rooted view of an unrooted tree.
std::string write_subtree(const SpeciesTree &tree, int node, int parent) {
  const auto &current = tree.nodes[node];
  std::vector<int> children;
  for (int neighbor : current.neighbors) {
    if (neighbor != parent) {
      children.push_back(neighbor);
    }
  }
  if (!current.label.empty() && children.empty()) {
    return escape_newick_label(current.label);
  }
  std::vector<std::string> parts;
  parts.reserve(children.size());
  for (int child : children) {
    parts.push_back(write_subtree(tree, child, node));
  }
  std::sort(parts.begin(), parts.end());
  std::ostringstream out;
  out << "(";
  for (std::size_t i = 0; i < parts.size(); ++i) {
    if (i != 0) {
      out << ",";
    }
    out << parts[i];
  }
  out << ")";
  if (!current.label.empty()) {
    out << escape_newick_label(current.label);
  }
  return out.str();
}

} // namespace

/// Serializes a species tree in deterministic Newick order.
std::string to_newick(const SpeciesTree &tree) {
  if (tree.nodes.empty()) {
    return ";";
  }
  if (tree.nodes.size() == 1) {
    return escape_newick_label(tree.nodes[0].label) + ";";
  }
  if (tree.nodes.size() == 2 && !tree.nodes[0].label.empty() &&
      !tree.nodes[1].label.empty()) {
    std::vector<std::string> labels = {
        escape_newick_label(tree.nodes[0].label),
        escape_newick_label(tree.nodes[1].label)};
    std::sort(labels.begin(), labels.end());
    return "(" + labels[0] + "," + labels[1] + ");";
  }
  int root = -1;
  for (int i = 0; i < static_cast<int>(tree.nodes.size()); ++i) {
    if (tree.nodes[i].neighbors.size() > 1 || tree.nodes[i].label.empty()) {
      root = i;
      break;
    }
  }
  if (root == -1) {
    root = 0;
  }
  return write_subtree(tree, root, -1) + ";";
}

} // namespace wqfm

namespace wqfm::internal {

/// Builds the taxon registry from input gene trees.
TaxonRegistry build_registry(const std::vector<GeneTree> &genes) {
  std::vector<std::string> labels;
  for (const GeneTree &gene : genes) {
    auto leaf_labels = gene.leaf_labels();
    labels.insert(labels.end(), leaf_labels.begin(), leaf_labels.end());
  }
  if (labels.empty()) {
    throw std::invalid_argument("no taxa found in input gene trees");
  }
  return TaxonRegistry(std::move(labels));
}

/// Builds parent/order/interval arrays for an indexed tree.
RootedTraversal build_rooted_traversal(const IndexedTree &tree) {
  RootedTraversal out;
  const int n = static_cast<int>(tree.nodes.size());
  out.parent.assign(n, -1);
  out.tin.assign(n, 0);
  out.tout.assign(n, 0);
  if (n == 0) {
    return out;
  }
  const int root = tree.root >= 0 ? tree.root : 0;
  std::vector<std::pair<int, int>> stack{{root, 0}};
  out.parent[root] = root;
  while (!stack.empty()) {
    const int node = stack.back().first;
    int &next = stack.back().second;
    if (next == 0) {
      out.tin[node] = static_cast<int>(out.preorder.size());
      out.preorder.push_back(node);
    }
    if (next < static_cast<int>(tree.nodes[node].neighbors.size())) {
      const int nb = tree.nodes[node].neighbors[next++];
      if (nb == out.parent[node]) {
        continue;
      }
      out.parent[nb] = node;
      stack.push_back({nb, 0});
      continue;
    }
    out.tout[node] = static_cast<int>(out.preorder.size());
    out.postorder.push_back(node);
    stack.pop_back();
  }
  return out;
}

/// Mixes three integers into a deterministic 64-bit hash.
uint64_t stable_hash_u64(uint64_t a, uint64_t b, uint64_t c) {
  uint64_t h = 1469598103934665603ULL;
  auto mix = [&](uint64_t x) {
    h ^= x + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
    h *= 1099511628211ULL;
  };
  mix(a);
  mix(b);
  mix(c);
  return h;
}

} // namespace wqfm::internal
