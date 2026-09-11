#include "internal.hpp"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace wqfm::internal {
namespace {

/// Fills traversal intervals for the consensus tree.
void fill_consensus_traversal(ConsensusTree& consensus) {
  consensus.preorder.clear();
  consensus.tin.assign(consensus.nodes.size(), 0);
  consensus.tout.assign(consensus.nodes.size(), 0);
  if (consensus.root < 0) return;
  std::vector<std::pair<int, int>> stack{{consensus.root, 0}};
  while (!stack.empty()) {
    const int node = stack.back().first;
    int& next = stack.back().second;
    if (next == 0) {
      consensus.tin[static_cast<std::size_t>(node)] =
          static_cast<int>(consensus.preorder.size());
      consensus.preorder.push_back(node);
    }
    if (next < static_cast<int>(
                   consensus.nodes[static_cast<std::size_t>(node)]
                       .children.size())) {
      stack.push_back(
          {consensus.nodes[static_cast<std::size_t>(node)]
               .children[static_cast<std::size_t>(next++)],
           0});
      continue;
    }
    consensus.tout[static_cast<std::size_t>(node)] =
        static_cast<int>(consensus.preorder.size());
    stack.pop_back();
  }
}

/// Escapes a consensus leaf label for Newick output.
std::string escape_consensus_label(const std::string& label) {
  bool quote = label.empty();
  for (char ch : label) {
    if (std::isspace(static_cast<unsigned char>(ch)) || ch == '(' ||
        ch == ')' || ch == ',' || ch == ':' || ch == ';' || ch == '\'' ||
        ch == '[' || ch == ']') {
      quote = true;
      break;
    }
  }
  if (!quote) return label;
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

/// Serializes a consensus subtree as Newick.
std::string consensus_subtree_newick(const ConsensusTree& consensus,
                                     const TaxonRegistry& taxa,
                                     int node) {
  if (consensus.nodes[static_cast<std::size_t>(node)].taxon >= 0) {
    return escape_consensus_label(
        taxa.label_for(consensus.nodes[static_cast<std::size_t>(node)].taxon));
  }
  std::vector<std::string> parts;
  parts.reserve(
      consensus.nodes[static_cast<std::size_t>(node)].children.size());
  for (int child :
       consensus.nodes[static_cast<std::size_t>(node)].children) {
    parts.push_back(consensus_subtree_newick(consensus, taxa, child));
  }
  std::sort(parts.begin(), parts.end());
  std::ostringstream out;
  out << "(";
  for (int i = 0; i < static_cast<int>(parts.size()); ++i) {
    if (i != 0) out << ",";
    out << parts[static_cast<std::size_t>(i)];
  }
  out << ")";
  return out.str();
}

}  // namespace

/// Converts a parsed Newick tree into the rooted consensus representation.
ConsensusTree consensus_from_tree(const Tree& tree, const TaxonRegistry& taxa) {
  if (tree.nodes.empty()) {
    throw std::invalid_argument("external consensus tree is empty");
  }

  ConsensusTree out;
  out.nodes.resize(tree.nodes.size());
  out.root = tree.root >= 0 ? tree.root : 0;
  out.taxon_to_leaf.assign(taxa.size(), -1);

  std::vector<int> parent(tree.nodes.size(), -1);
  std::vector<int> stack{out.root};
  parent[static_cast<std::size_t>(out.root)] = out.root;
  while (!stack.empty()) {
    const int node = stack.back();
    stack.pop_back();
    out.nodes[static_cast<std::size_t>(node)].parent =
        node == out.root ? -1 : parent[static_cast<std::size_t>(node)];

    int child_count = 0;
    for (int nb : tree.nodes[static_cast<std::size_t>(node)].neighbors) {
      if (nb == parent[static_cast<std::size_t>(node)]) continue;
      if (parent[static_cast<std::size_t>(nb)] != -1) {
        throw std::invalid_argument("external consensus tree is not acyclic");
      }
      parent[static_cast<std::size_t>(nb)] = node;
      out.nodes[static_cast<std::size_t>(node)].children.push_back(nb);
      stack.push_back(nb);
      ++child_count;
    }

    if (!tree.nodes[static_cast<std::size_t>(node)].label.empty() &&
        child_count == 0) {
      const int taxon =
          taxa.id_for(tree.nodes[static_cast<std::size_t>(node)].label);
      if (out.taxon_to_leaf[static_cast<std::size_t>(taxon)] != -1) {
        throw std::invalid_argument(
            "duplicate taxon in external consensus: " +
            tree.nodes[static_cast<std::size_t>(node)].label);
      }
      out.nodes[static_cast<std::size_t>(node)].taxon = taxon;
      out.taxon_to_leaf[static_cast<std::size_t>(taxon)] = node;
    }
  }
  for (int i = 0; i < static_cast<int>(parent.size()); ++i) {
    if (parent[static_cast<std::size_t>(i)] == -1) {
      throw std::invalid_argument("external consensus tree is disconnected");
    }
  }
  for (int taxon = 0; taxon < taxa.size(); ++taxon) {
    if (out.taxon_to_leaf[static_cast<std::size_t>(taxon)] < 0) {
      throw std::invalid_argument("external consensus is missing taxon: " +
                                  taxa.label_for(taxon));
    }
  }

  fill_consensus_traversal(out);
  return out;
}

/// Writes the consensus tree in Newick form when requested by Options.
void write_consensus_dump(const Options& options, const SolverContext& ctx) {
  const ConsensusTree& consensus = context_consensus(ctx);
  if (options.dump_consensus.empty() || consensus.root < 0) return;
  std::ofstream out(options.dump_consensus);
  if (!out) return;
  out << consensus_subtree_newick(consensus, context_taxa(ctx), consensus.root)
      << ";\n";
}

}  // namespace wqfm::internal
