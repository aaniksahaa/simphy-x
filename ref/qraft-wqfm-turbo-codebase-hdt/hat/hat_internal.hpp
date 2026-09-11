#pragma once

#include "internal.hpp"

#include <algorithm>
#include <array>
#include <utility>

namespace wqfm::internal {

struct HatTaxonIndex {
  std::vector<int> taxa;
  std::vector<std::pair<int, int>> by_taxon;

  int find(int taxon) const {
    const auto it = std::lower_bound(
        by_taxon.begin(), by_taxon.end(), taxon,
        [](const std::pair<int, int>& entry, int value) {
          return entry.first < value;
        });
    return it != by_taxon.end() && it->first == taxon ? it->second : -1;
  }
};

struct HatActiveSet {
  std::vector<int> taxa;
  HatTaxonIndex active_by_taxon;
};

struct HatSplitChoice {
  Candidate candidate;
  int candidate_index = -1;
  double score = 0.0;
  int representative_order = 0;
  int candidate_count = 0;
  int left_count = 0;
  int right_count = 0;
};

struct HatActiveSplitSelectionProbe {
  double score = 0.0;
  int left_count = 0;
  int representative_order = 0;
  int edge_child_node = 0;
  std::array<int, 4> leaf_intervals{};
};

struct HatFourColorCounts {
  double q01_23 = 0.0;
  double q02_13 = 0.0;
  double q03_12 = 0.0;
};

inline int hat_hdt_workspace_worker_limit(const Options& options) {
  (void)options;
  return current_inner_parallel_worker_limit();
}

// Give every input gene its own score row. The scheduler may execute rows in
// any order, while the final gene-order reduction remains deterministic and
// exposes all available gene-level parallelism without a fixed worker ceiling.
inline int hat_score_gene_partition_count(int gene_count) {
  return std::max(0, gene_count);
}

struct HatQueryTree {
  struct Node {
    int first_adj = -1;
    int degree = 0;
    int taxon_active = -1;
    int source_node = -1;
  };

  std::vector<Node> nodes;
  std::vector<int> adjacency;
  std::vector<int> parent;
  std::vector<int> child_first;
  std::vector<int> child_count;
  std::vector<int> child_nodes;
  std::vector<int> preorder;
  std::vector<int> postorder;
  std::vector<int> tin;
  std::vector<int> tout;
  std::vector<int> subtree_size;
  std::vector<int> subtree_leaf_first;
  std::vector<int> subtree_leaf_count;
  std::vector<int> leaf_active_order;
  std::vector<int> heavy_child;
  std::vector<int> active_to_node;
  int root = -1;
  int leaf_count = 0;
};

struct HatMedianBranchColoring {
  std::array<int, 3> branch_nodes{-1, -1, -1};
  std::array<std::uint8_t, 3> branch_colors{0, 1, 2};
  int temporary_color_two_node = -1;
};

inline int hat_query_leaf_count(const HatQueryTree& query, int node) {
  return node >= 0 && node < static_cast<int>(query.subtree_leaf_count.size())
             ? query.subtree_leaf_count[static_cast<std::size_t>(node)]
             : 0;
}

inline HatMedianBranchColoring hat_binary_median_branch_coloring(
    const HatQueryTree& query,
    int node) {
  HatMedianBranchColoring coloring;
  if (node < 0 || node >= static_cast<int>(query.nodes.size()) ||
      node >= static_cast<int>(query.child_first.size())) {
    return coloring;
  }
  const int child_begin = query.child_first[static_cast<std::size_t>(node)];
  if (node == query.root) {
    if (query.child_count[static_cast<std::size_t>(node)] != 3) {
      return coloring;
    }
    coloring.branch_nodes = std::array<int, 3>{
        query.child_nodes[static_cast<std::size_t>(child_begin)],
        query.child_nodes[static_cast<std::size_t>(child_begin + 1)],
        query.child_nodes[static_cast<std::size_t>(child_begin + 2)]};
    coloring.branch_colors = std::array<std::uint8_t, 3>{0, 1, 2};
    coloring.temporary_color_two_node = coloring.branch_nodes[2];
    return coloring;
  }
  if (query.child_count[static_cast<std::size_t>(node)] != 2) {
    return coloring;
  }
  const int child0 = query.child_nodes[static_cast<std::size_t>(child_begin)];
  const int child1 =
      query.child_nodes[static_cast<std::size_t>(child_begin + 1)];
  int color_one_child = child0;
  int color_two_child = child1;
  if (hat_query_leaf_count(query, child0) <
      hat_query_leaf_count(query, child1)) {
    color_one_child = child1;
    color_two_child = child0;
  }
  coloring.branch_nodes = std::array<int, 3>{
      query.parent[static_cast<std::size_t>(node)], color_one_child,
      color_two_child};
  coloring.branch_colors = std::array<std::uint8_t, 3>{0, 1, 2};
  coloring.temporary_color_two_node = color_two_child;
  return coloring;
}

struct HatNniMove {
  int edge_a = -1;
  int edge_b = -1;
  int swap_a_neighbor = -1;
  int swap_b_neighbor = -1;
  int alternative = 0;
  double gain = 0.0;
};

HatActiveSet make_hat_active_set(const std::vector<Identity>& identities,
                                 int taxon_count);

HatTaxonIndex make_hat_taxon_index(const std::vector<int>& taxa,
                                   int taxon_count);

int hat_active_index(const HatActiveSet& active, int taxon);

void note_hat_source_bank_scoring_use(const std::vector<GeneIndex>& genes);

std::vector<int> hat_taxa_from_identities(
    const std::vector<Identity>& identities);

HatQueryTree make_hat_query_tree_from_consensus_induced(
    const ConsensusTree& consensus,
    const HatActiveSet& active);

HatQueryTree make_hat_query_tree_from_solver_tree(
    const SolverTree& tree,
    const std::vector<int>& active_taxa,
    int taxon_count);

std::vector<int> query_leaf_owner_order(
    const HatQueryTree& query,
    const std::vector<int>& query_active_to_owner);

SolverTree hat_base_tree(const std::vector<Identity>& identities);

std::vector<int> hat_active_neighbors(const SolverTree& tree, int node);

bool hat_has_edge(const SolverTree& tree, int a, int b);

void hat_add_edge(SolverTree& tree, int a, int b);

void hat_remove_edge(SolverTree& tree, int a, int b);

SolverTree hat_glue(SolverTree left,
                    SolverTree right,
                    AttachmentEdge left_edge,
                    AttachmentEdge right_edge);

// Selects a structural-child failure by logical side after both tasks join.
void rethrow_hat_parallel_child_failure(
    const std::exception_ptr& first,
    const std::exception_ptr& second);

HatSplitChoice hat_choose_active_split_hdt(
    const SolverContext& ctx,
    const std::vector<Identity>& identities,
    uint64_t path,
    int depth);

int hat_select_active_split_candidate_for_tests(
    const std::vector<HatActiveSplitSelectionProbe>& probes,
    int active_count,
    int min_side);

AttachmentEdge hat_best_one_x_edge_hdt(
    const SolverContext& ctx,
    const SolverTree& child_tree,
    const std::vector<int>& active_taxa,
    const std::vector<int>& context_taxa,
    const std::vector<int>* identity_taxa);

std::pair<AttachmentEdge, AttachmentEdge> hat_best_one_x_edges_hdt(
    const SolverContext& ctx,
    const SolverTree& left_tree,
    const SolverTree& right_tree,
    const std::vector<int>& left_taxa,
    const std::vector<int>& right_taxa,
    const std::vector<int>* identity_taxa);

std::vector<double> hat_score_active_split_edges_hdt(
    const SolverContext& ctx,
    const HatActiveSet& active,
    const HatQueryTree& query,
    const std::vector<int>& edge_child_nodes);

HatNniMove hat_best_nni_move_hdt(const SolverContext& ctx,
                                 const SolverTree& tree,
                                 const std::vector<int>& active_taxa);

void run_hat_nni_repair(SolverContext& ctx,
                        SolverTree& tree,
                        const std::vector<int>& active_taxa);

}  // namespace wqfm::internal
