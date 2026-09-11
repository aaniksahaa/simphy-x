#include "hat/hat_internal.hpp"
#include "hat/hat_score_selection.hpp"
#include "hat/hdt/hdt_traversal.hpp"
#include "candidates/split_records.hpp"
#include "hat/hdt/owner_projection.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iterator>
#include <limits>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>

namespace wqfm::internal {
namespace {

constexpr double kHatNniGainEps = 1e-9;
std::size_t checked_hat_nni_size_product(std::size_t first,
                                         std::size_t second,
                                         const char* what) {
  if (first != 0 && second > std::numeric_limits<std::size_t>::max() / first) {
    throw std::length_error(what);
  }
  return first * second;
}

int checked_hat_nni_int_size(std::size_t size, const char* what) {
  if (size > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
    throw std::length_error(what);
  }
  return static_cast<int>(size);
}

std::size_t checked_hat_nni_doubled_capacity(std::size_t size,
                                              bool add_one,
                                              const char* what) {
  const std::size_t extra = add_one ? 1u : 0u;
  if (size > (std::numeric_limits<std::size_t>::max() - extra) / 2u) {
    throw std::length_error(what);
  }
  return size * 2u + extra;
}

long long checked_hat_nni_doubled_counter(std::size_t size,
                                           const char* what) {
  if (size > static_cast<std::size_t>(
                 std::numeric_limits<long long>::max() / 2ll)) {
    throw std::length_error(what);
  }
  return static_cast<long long>(size) * 2ll;
}

struct HatNniEdgeProposal {
  int p_query = -1;
  int q_query = -1;
  int q_split_query = -1;
  std::array<std::uint8_t, 2> p_colors{0, 1};
  std::array<std::uint8_t, 2> q_colors{2, 3};
  HatNniMove swap_p0_q0;
  HatNniMove swap_p0_q1;
};

using HatNniSplitKey = SplitKey256;
using HatNniTaxonFingerprints = SplitTaxonFingerprints256;

struct HatNniSplitSupportCache {
  bool enabled = false;
  int taxon_count = 0;
  HatNniTaxonFingerprints fingerprints;
  std::vector<SplitAggregate256> splits;
};

enum class HatNniGateDecision { Allow, StopQuartet, StopSplit };

uint64_t hat_nni_node_edge_key(int a, int b) {
  const uint32_t lo = static_cast<uint32_t>(std::min(a, b));
  const uint32_t hi = static_cast<uint32_t>(std::max(a, b));
  return (static_cast<uint64_t>(lo) << 32) | static_cast<uint64_t>(hi);
}

bool hat_nni_gain_meets_threshold(double gain, double min_gain) {
  return gain > kHatNniGainEps && gain + kHatNniGainEps >= min_gain;
}

double hat_nni_min_gain_threshold(const SolverContext& ctx) {
  const double per_gene = std::max(0.0, ctx.options.hat_nni_min_gain_per_gene);
  const int gene_count = checked_hat_nni_int_size(
      scoring_genes(ctx).size(), "NNI gene count exceeds its integer domain");
  return std::max(kHatNniGainEps,
                  per_gene * static_cast<double>(std::max(1, gene_count)));
}

void hat_atomic_double_add(std::atomic<double>& counter, double delta) {
  double current = counter.load(std::memory_order_relaxed);
  while (!counter.compare_exchange_weak(current, current + delta,
                                        std::memory_order_relaxed,
                                        std::memory_order_relaxed)) {
  }
}

void hat_atomic_double_max(std::atomic<double>& counter, double value) {
  double current = counter.load(std::memory_order_relaxed);
  while (current < value &&
         !counter.compare_exchange_weak(current, value,
                                        std::memory_order_relaxed,
                                        std::memory_order_relaxed)) {
  }
}

void hat_note_accepted_nni_gain(double gain) {
  benchmark_stats().hat_nni_last_accepted_gain.store(
      gain, std::memory_order_relaxed);
  hat_atomic_double_max(benchmark_stats().hat_nni_max_accepted_gain, gain);
}

std::array<uint64_t, 4> hat_nni_sub_words(
    const std::array<uint64_t, 4>& a,
    const std::array<uint64_t, 4>& b) {
  std::array<uint64_t, 4> out{};
  for (int i = 0; i < 4; ++i) out[static_cast<std::size_t>(i)] =
      a[static_cast<std::size_t>(i)] - b[static_cast<std::size_t>(i)];
  return out;
}

HatNniSplitSupportCache build_hat_nni_split_support_cache(
    const SolverContext& ctx) {
  const auto started = std::chrono::steady_clock::now();
  HatNniSplitSupportCache cache;
  if (!ctx.options.hat_nni_prefix_gate) return cache;
  cache.enabled = true;
  cache.taxon_count = checked_hat_nni_int_size(
      context_taxa(ctx).size(), "NNI taxon count exceeds its integer domain");
  cache.fingerprints =
      build_split_taxon_fingerprints256(cache.taxon_count, ctx.options);
  const std::vector<GeneIndex>& genes = scoring_genes(ctx);
  const int gene_count = checked_hat_nni_int_size(
      genes.size(), "NNI gene count exceeds its integer domain");
  const int chunks =
      std::max(1, std::min(gene_count,
                           hat_hdt_workspace_worker_limit(ctx.options)));
  std::vector<std::vector<SplitOccurrence256>> local_records(
      static_cast<std::size_t>(chunks));

  parallel_for_inner(0, chunks, [&](int chunk) {
    const int begin = static_cast<int>(
        static_cast<int64_t>(gene_count) * chunk / chunks);
    const int end = static_cast<int>(
        static_cast<int64_t>(gene_count) * (chunk + 1) / chunks);
    auto& local = local_records[static_cast<std::size_t>(chunk)];
    std::size_t capacity = 0;
    for (int g = begin; g < end; ++g) {
      const std::size_t nodes = genes[static_cast<std::size_t>(g)].nodes.size();
      if (nodes > std::numeric_limits<std::size_t>::max() - capacity) {
        throw std::length_error(
            "NNI split records exceed the addressable size domain");
      }
      capacity += nodes;
    }
    local.reserve(capacity);
    std::vector<uint32_t> subtree_size;
    std::vector<uint8_t> contains_ref;
    std::vector<std::array<uint64_t, 4>> subtree_sum;
    std::vector<std::array<uint64_t, 4>> duplicate_balance;

    for (int g = begin; g < end; ++g) {
      const GeneIndex& gene = genes[static_cast<std::size_t>(g)];
      if (gene.nodes.empty()) continue;
      if (gene.nodes.size() >
          static_cast<std::size_t>(std::numeric_limits<int>::max())) {
        throw std::length_error("NNI split cache exceeds its integer node domain");
      }
      const GeneTraversal& tr = gene.reference_traversal;
      const int node_count = static_cast<int>(gene.nodes.size());
      const GeneSplitLeafProjection projection =
          make_gene_split_leaf_projection(gene, cache.taxon_count,
                                          cache.fingerprints.word);
      subtree_size.assign(static_cast<std::size_t>(node_count), 0);
      contains_ref.assign(static_cast<std::size_t>(node_count), 0);
      subtree_sum.assign(static_cast<std::size_t>(node_count),
                         std::array<uint64_t, 4>{});
      duplicate_balance.assign(static_cast<std::size_t>(node_count),
                               std::array<uint64_t, 4>{});
      for (int node : tr.postorder) {
        subtree_size[static_cast<std::size_t>(node)] =
            projection.unique_size[static_cast<std::size_t>(node)];
        contains_ref[static_cast<std::size_t>(node)] =
            projection.reference_marker[static_cast<std::size_t>(node)];
        subtree_sum[static_cast<std::size_t>(node)] =
            projection.unique_fingerprint[static_cast<std::size_t>(node)];
        duplicate_balance[static_cast<std::size_t>(node)] =
            projection.duplicate_balance[static_cast<std::size_t>(node)];
        for (int nb : gene.nodes[static_cast<std::size_t>(node)].neighbors) {
          if (tr.parent[static_cast<std::size_t>(nb)] != node) continue;
          subtree_size[static_cast<std::size_t>(node)] +=
              subtree_size[static_cast<std::size_t>(nb)];
          contains_ref[static_cast<std::size_t>(node)] =
              static_cast<uint8_t>(
                  contains_ref[static_cast<std::size_t>(node)] |
                  contains_ref[static_cast<std::size_t>(nb)]);
          for (int word = 0; word < 4; ++word) {
            subtree_sum[static_cast<std::size_t>(node)]
                       [static_cast<std::size_t>(word)] +=
                subtree_sum[static_cast<std::size_t>(nb)]
                           [static_cast<std::size_t>(word)];
            duplicate_balance[static_cast<std::size_t>(node)]
                             [static_cast<std::size_t>(word)] +=
                duplicate_balance[static_cast<std::size_t>(nb)]
                                 [static_cast<std::size_t>(word)];
          }
        }
      }

      uint32_t edge_id = 0;
      for (int child : tr.preorder) {
        if (child == gene.reference_root) continue;
        if (contains_ref[static_cast<std::size_t>(child)] != 0) {
          ++edge_id;
          continue;
        }
        if (duplicate_balance[static_cast<std::size_t>(child)] !=
            std::array<uint64_t, 4>{}) {
          ++edge_id;
          continue;
        }
        const uint32_t size = subtree_size[static_cast<std::size_t>(child)];
        if (projection.distinct_taxa < 4 || size < 2 ||
            size > projection.distinct_taxa - 2) {
          ++edge_id;
          continue;
        }
        HatNniSplitKey key;
        key.size = size;
        key.h = subtree_sum[static_cast<std::size_t>(child)];
        local.push_back(SplitOccurrence256{
            key, static_cast<uint32_t>(g), edge_id,
            static_cast<uint32_t>(child), 0});
        ++edge_id;
      }
    }
  }, 1);

  std::size_t observed = 0;
  for (const auto& local : local_records) {
    if (local.size() > std::numeric_limits<std::size_t>::max() - observed) {
      throw std::length_error("NNI split cache exceeds the addressable size domain");
    }
    observed += local.size();
  }
  std::vector<SplitOccurrence256> occurrences;
  occurrences.reserve(observed);
  for (auto& local : local_records) {
    occurrences.insert(occurrences.end(),
                       std::make_move_iterator(local.begin()),
                       std::make_move_iterator(local.end()));
    std::vector<SplitOccurrence256>().swap(local);
  }
  SplitRadixStats radix_stats;
  cache.splits = aggregate_split_occurrences256(occurrences, radix_stats);
  benchmark_counter_max(benchmark_stats().hat_nni_split_support_unique_splits,
                        static_cast<long long>(cache.splits.size()));
  benchmark_timer_add(benchmark_stats().hat_nni_split_support_ns,
                      std::chrono::duration_cast<std::chrono::nanoseconds>(
                          std::chrono::steady_clock::now() - started));
  return cache;
}

std::string hat_nni_pass_tree_name(int pass) {
  std::ostringstream name;
  name << "pass_" << std::setw(6) << std::setfill('0') << pass << ".tre";
  return name.str();
}

std::vector<int> snapshot_active_neighbors(const SolverTree& tree, int node) {
  std::vector<int> neighbors;
  if (node < 0 || node >= static_cast<int>(tree.nodes.size()) ||
      !tree.nodes[static_cast<std::size_t>(node)].active) {
    return neighbors;
  }
  for (int nb : tree.nodes[static_cast<std::size_t>(node)].neighbors) {
    if (nb >= 0 && nb < static_cast<int>(tree.nodes.size()) &&
        tree.nodes[static_cast<std::size_t>(nb)].active) {
      neighbors.push_back(nb);
    }
  }
  return neighbors;
}

void snapshot_remove_edge(SolverTree& tree, int a, int b) {
  auto erase_one = [](std::vector<int>& values, int value) {
    values.erase(std::remove(values.begin(), values.end(), value),
                 values.end());
  };
  erase_one(tree.nodes[static_cast<std::size_t>(a)].neighbors, b);
  erase_one(tree.nodes[static_cast<std::size_t>(b)].neighbors, a);
}

void snapshot_add_edge(SolverTree& tree, int a, int b) {
  auto has = [](const std::vector<int>& values, int value) {
    return std::find(values.begin(), values.end(), value) != values.end();
  };
  if (!has(tree.nodes[static_cast<std::size_t>(a)].neighbors, b)) {
    tree.nodes[static_cast<std::size_t>(a)].neighbors.push_back(b);
  }
  if (!has(tree.nodes[static_cast<std::size_t>(b)].neighbors, a)) {
    tree.nodes[static_cast<std::size_t>(b)].neighbors.push_back(a);
  }
}

SpeciesTree snapshot_solver_tree_to_species(const SolverTree& tree,
                                            const TaxonRegistry& taxa) {
  SolverTree compact = tree;
  bool changed = true;
  while (changed) {
    changed = false;
    for (int i = 0; i < static_cast<int>(compact.nodes.size()); ++i) {
      if (!compact.nodes[static_cast<std::size_t>(i)].active) continue;
      const std::vector<int> neighbors = snapshot_active_neighbors(compact, i);
      if (neighbors.size() != 2) continue;
      snapshot_remove_edge(compact, i, neighbors[0]);
      snapshot_remove_edge(compact, i, neighbors[1]);
      snapshot_add_edge(compact, neighbors[0], neighbors[1]);
      compact.nodes[static_cast<std::size_t>(i)].active = false;
      compact.nodes[static_cast<std::size_t>(i)].neighbors.clear();
      changed = true;
    }
  }

  SpeciesTree species;
  std::vector<int> map(compact.nodes.size(), -1);
  for (int i = 0; i < static_cast<int>(compact.nodes.size()); ++i) {
    const SolverTree::Node& source = compact.nodes[static_cast<std::size_t>(i)];
    if (!source.active || source.identity_id >= static_cast<int64_t>(taxa.size())) {
      continue;
    }
    map[static_cast<std::size_t>(i)] = static_cast<int>(species.nodes.size());
    SpeciesTree::Node node;
    if (source.identity_id >= 0) {
      node.label = taxa.label_for(static_cast<int>(source.identity_id));
    }
    species.nodes.push_back(std::move(node));
  }
  for (int i = 0; i < static_cast<int>(compact.nodes.size()); ++i) {
    if (map[static_cast<std::size_t>(i)] < 0) continue;
    for (int nb : compact.nodes[static_cast<std::size_t>(i)].neighbors) {
      if (nb >= 0 && nb < static_cast<int>(map.size()) &&
          map[static_cast<std::size_t>(nb)] >= 0) {
        species.nodes[static_cast<std::size_t>(map[static_cast<std::size_t>(i)])]
            .neighbors.push_back(map[static_cast<std::size_t>(nb)]);
      }
    }
  }
  return species;
}

class HatNniPassRecorder {
 public:
  HatNniPassRecorder(const SolverContext& ctx, const std::string& mode)
      : ctx_(ctx), mode_(mode) {
    if (ctx.options.hat_nni_pass_tree_dir.empty()) return;
    dir_ = ctx.options.hat_nni_pass_tree_dir;
    std::filesystem::create_directories(dir_);
    tsv_.open(dir_ / "passes.tsv", std::ios::trunc);
    if (!tsv_) {
      throw std::runtime_error("failed to open HAT NNI pass log: " +
                               (dir_ / "passes.tsv").string());
    }
    tsv_ << "pass\tmode\tpass_moves\taccepted_moves\tpass_gain"
         << "\tcumulative_gain\tsplit_support_delta"
         << "\tcumulative_split_support_delta"
         << "\tedges_scored\tmin_gain\ttree\n";
  }

  bool enabled() const { return tsv_.is_open(); }

  void record(const SolverTree& tree,
              int pass,
              long long pass_moves,
              long long accepted_moves,
              double pass_gain,
              double cumulative_gain,
              int split_support_delta,
              long long cumulative_split_support_delta,
              long long edges_scored,
              double min_gain) {
    if (!enabled()) return;
    const std::filesystem::path tree_path = dir_ / hat_nni_pass_tree_name(pass);
    std::ofstream tree_out(tree_path);
    if (!tree_out) {
      throw std::runtime_error("failed to open HAT NNI pass tree: " +
                               tree_path.string());
    }
    tree_out << to_newick(snapshot_solver_tree_to_species(tree,
                                                          context_taxa(ctx_)))
             << "\n";
    tsv_ << pass << '\t' << mode_ << '\t' << pass_moves << '\t'
         << accepted_moves << '\t' << std::setprecision(17) << pass_gain
         << '\t' << cumulative_gain << '\t' << split_support_delta << '\t'
         << cumulative_split_support_delta << '\t' << edges_scored << '\t'
         << min_gain << '\t' << tree_path.string() << '\n';
    tsv_.flush();
  }

 private:
  const SolverContext& ctx_;
  std::string mode_;
  std::filesystem::path dir_;
  std::ofstream tsv_;
};

std::vector<int> query_children(const HatQueryTree& query, int node) {
  std::vector<int> children;
  if (node < 0 || node >= static_cast<int>(query.nodes.size())) {
    return children;
  }
  const int first = query.nodes[static_cast<std::size_t>(node)].first_adj;
  const int degree = query.nodes[static_cast<std::size_t>(node)].degree;
  children.reserve(static_cast<std::size_t>(degree));
  for (int i = 0; i < degree; ++i) {
    const int nb = query.adjacency[static_cast<std::size_t>(first + i)];
    if (query.parent[static_cast<std::size_t>(nb)] == node) {
      children.push_back(nb);
    }
  }
  return children;
}

bool query_internal_binary(const HatQueryTree& query, int node) {
  return node >= 0 && node < static_cast<int>(query.nodes.size()) &&
         query.nodes[static_cast<std::size_t>(node)].taxon_active < 0 &&
         query.nodes[static_cast<std::size_t>(node)].degree == 3;
}

std::array<int, 3> sorted_pair_pair(int a, int b, int c, int d) {
  std::array<int, 2> left{std::min(a, b), std::max(a, b)};
  std::array<int, 2> right{std::min(c, d), std::max(c, d)};
  if (right < left) std::swap(left, right);
  return std::array<int, 3>{left[0], left[1], right[1]};
}

double pair_pair_score(const HatFourColorCounts& counts,
                       int a,
                       int b,
                       int c,
                       int d) {
  const std::array<int, 3> key = sorted_pair_pair(a, b, c, d);
  if (key == std::array<int, 3>{0, 1, 3}) return counts.q01_23;
  if (key == std::array<int, 3>{0, 2, 3}) return counts.q02_13;
  if (key == std::array<int, 3>{0, 3, 2}) return counts.q03_12;
  return 0.0;
}

int query_source_node(const HatQueryTree& query, int node) {
  if (node < 0 || node >= static_cast<int>(query.nodes.size())) return -1;
  return query.nodes[static_cast<std::size_t>(node)].source_node;
}

std::vector<HatNniEdgeProposal> collect_nni_proposals(
    const SolverTree& tree,
    const HatQueryTree& query,
    const std::vector<uint64_t>* allowed_central_edges = nullptr) {
  std::vector<HatNniEdgeProposal> proposals;
  if (query.root < 0 || query.nodes.empty()) return proposals;
  for (int q : query.preorder) {
    if (q == query.root || !query_internal_binary(query, q)) continue;
    const int p = query.parent[static_cast<std::size_t>(q)];
    if (!query_internal_binary(query, p)) continue;

    const int edge_a = query_source_node(query, p);
    const int edge_b = query_source_node(query, q);
    if (edge_a < 0 || edge_b < 0 || !hat_has_edge(tree, edge_a, edge_b)) {
      continue;
    }
    if (allowed_central_edges != nullptr &&
        !std::binary_search(allowed_central_edges->begin(),
                            allowed_central_edges->end(),
                            hat_nni_node_edge_key(edge_a, edge_b))) {
      continue;
    }

    const HatMedianBranchColoring p_coloring =
        hat_binary_median_branch_coloring(query, p);
    const std::array<int, 3> p_branches = p_coloring.branch_nodes;
    const std::array<std::uint8_t, 3> p_branch_colors =
        p_coloring.branch_colors;
    int q_branch_index = -1;
    for (int i = 0; i < 3; ++i) {
      if (p_branches[static_cast<std::size_t>(i)] == q) {
        q_branch_index = i;
        break;
      }
    }
    if (q_branch_index < 0) continue;

    std::array<int, 2> p_side_query{-1, -1};
    std::array<std::uint8_t, 2> p_colors{0, 1};
    int p_side_count = 0;
    for (int i = 0; i < 3; ++i) {
      if (i == q_branch_index) continue;
      p_side_query[static_cast<std::size_t>(p_side_count)] =
          p_branches[static_cast<std::size_t>(i)];
      p_colors[static_cast<std::size_t>(p_side_count)] =
          p_branch_colors[static_cast<std::size_t>(i)];
      ++p_side_count;
    }
    if (p_side_count != 2) continue;

    const std::vector<int> q_children = query_children(query, q);
    if (q_children.size() != 2) continue;
    const int q0 = q_children[0];
    const int q1 = q_children[1];
    const int q0_size = hat_query_leaf_count(query, q0);
    const int q1_size = hat_query_leaf_count(query, q1);
    const int q_split = q0_size <= q1_size ? q0 : q1;
    const int q_base = q_split == q0 ? q1 : q0;

    const int p0_source = query_source_node(query, p_side_query[0]);
    const int p1_source = query_source_node(query, p_side_query[1]);
    const int q0_source = query_source_node(query, q_base);
    const int q1_source = query_source_node(query, q_split);
    if (p0_source < 0 || p1_source < 0 || q0_source < 0 || q1_source < 0 ||
        !hat_has_edge(tree, edge_a, p0_source) ||
        !hat_has_edge(tree, edge_a, p1_source) ||
        !hat_has_edge(tree, edge_b, q0_source) ||
        !hat_has_edge(tree, edge_b, q1_source)) {
      continue;
    }
    HatNniEdgeProposal proposal;
    proposal.p_query = p;
    proposal.q_query = q;
    proposal.q_split_query = q_split;
    proposal.p_colors = p_colors;
    proposal.q_colors = std::array<std::uint8_t, 2>{
        p_branch_colors[static_cast<std::size_t>(q_branch_index)], 3};

    proposal.swap_p0_q0.edge_a = edge_a;
    proposal.swap_p0_q0.edge_b = edge_b;
    proposal.swap_p0_q0.swap_a_neighbor = p0_source;
    proposal.swap_p0_q0.swap_b_neighbor = q0_source;
    proposal.swap_p0_q0.alternative = 1;

    proposal.swap_p0_q1.edge_a = edge_a;
    proposal.swap_p0_q1.edge_b = edge_b;
    proposal.swap_p0_q1.swap_a_neighbor = p0_source;
    proposal.swap_p0_q1.swap_b_neighbor = q1_source;
    proposal.swap_p0_q1.alternative = 2;

    proposals.push_back(proposal);
  }
  return proposals;
}

struct HatNniProposalScores {
  double current = 0.0;
  double swap_p0_q0 = 0.0;
  double swap_p0_q1 = 0.0;
};

struct HatScoredNniData {
  std::vector<HatNniEdgeProposal> proposals;
  std::vector<HatNniProposalScores> scores;
};

struct HatNniOwnerRange {
  int begin = 0;
  int end = 0;
};

struct HatNniPrefixGate {
  bool enabled = false;
  bool initialized = false;
  bool middle_gain_regime = false;
  int pass_cap = std::numeric_limits<int>::max();
  int gene_count = 1;
  double split_drop_threshold = -1.0;
  double max_pass_gain = 0.0;
};

HatNniPrefixGate make_hat_nni_prefix_gate(const SolverContext& ctx) {
  HatNniPrefixGate gate;
  gate.enabled = ctx.options.hat_nni_prefix_gate;
  gate.gene_count =
      std::max(1, static_cast<int>(scoring_genes(ctx).size()));
  gate.pass_cap = ctx.options.hat_nni_unbounded
                      ? std::numeric_limits<int>::max()
                      : std::max(0, ctx.options.hat_nni_passes);
  const int absolute_split_drop = ctx.options.hat_nni_gate_split_drop_v1;
  const double per_gene_split_drop =
      ctx.options.hat_nni_gate_split_drop_per_gene_v1;
  gate.split_drop_threshold =
      absolute_split_drop >= 0
          ? static_cast<double>(absolute_split_drop)
          : std::max(0.0, per_gene_split_drop) *
                static_cast<double>(gate.gene_count);
  return gate;
}

HatNniGateDecision hat_nni_prefix_gate_decide(const SolverContext& ctx,
                                              HatNniPrefixGate& gate,
                                              int accepted_passes,
                                              double pass_gain,
                                              int split_support_delta) {
  if (!gate.enabled) return HatNniGateDecision::Allow;
  const double gain_per_gene =
      pass_gain / static_cast<double>(std::max(1, gate.gene_count));
  if (!gate.initialized) {
    gate.initialized = true;
    if (gain_per_gene <
        ctx.options.hat_nni_gate_low_first_gain_per_gene) {
      gate.pass_cap = 0;
      return HatNniGateDecision::StopQuartet;
    }
    if (gain_per_gene >=
        ctx.options.hat_nni_gate_high_first_gain_per_gene) {
      gate.pass_cap =
          std::min(gate.pass_cap,
                   std::max(0, ctx.options.hat_nni_gate_high_gain_pass_cap));
    } else {
      gate.middle_gain_regime = true;
    }
  }
  if (accepted_passes >= gate.pass_cap) {
    return HatNniGateDecision::StopQuartet;
  }
  if (gate.middle_gain_regime) {
    const int next_pass = accepted_passes + 1;
    const double middle_min_gain =
        std::max(0.0, ctx.options.hat_nni_gate_middle_min_gain_per_gene) *
        static_cast<double>(std::max(1, gate.gene_count));
    if (pass_gain + kHatNniGainEps < middle_min_gain) {
      return HatNniGateDecision::StopQuartet;
    }
    if (next_pass > ctx.options.hat_nni_gate_middle_min_passes &&
        gate.max_pass_gain > 0.0 &&
        pass_gain < ctx.options.hat_nni_gate_middle_rel_max_gain *
                        gate.max_pass_gain) {
      return HatNniGateDecision::StopQuartet;
    }
  }
  if (gate.split_drop_threshold >= 0 &&
      split_support_delta < -gate.split_drop_threshold) {
    return HatNniGateDecision::StopSplit;
  }
  return HatNniGateDecision::Allow;
}

void hat_nni_prefix_gate_note_accepted(HatNniPrefixGate& gate,
                                       double pass_gain,
                                       int split_support_delta) {
  if (!gate.enabled) return;
  gate.max_pass_gain = std::max(gate.max_pass_gain, pass_gain);
  benchmark_counter_add(benchmark_stats().hat_nni_split_support_delta,
                        split_support_delta);
}

void hat_nni_prefix_gate_note_stop(HatNniGateDecision decision) {
  if (decision == HatNniGateDecision::Allow) return;
  benchmark_counter_add(benchmark_stats().hat_nni_prefix_gate_stops);
  if (decision == HatNniGateDecision::StopQuartet) {
    benchmark_counter_add(benchmark_stats().hat_nni_prefix_gate_quartet_stops);
  } else if (decision == HatNniGateDecision::StopSplit) {
    benchmark_counter_add(benchmark_stats().hat_nni_prefix_gate_split_stops);
  }
}

std::vector<HatNniProposalScores> score_nni_proposals_hdt(
    const SolverContext& ctx,
    const HatQueryTree& query,
    const std::vector<HatNniEdgeProposal>& proposals,
    const std::vector<int>& active_taxa,
    int taxon_count) {
  std::vector<HatNniProposalScores> totals(proposals.size());
  if (proposals.empty()) return totals;
  const int proposal_count = checked_hat_nni_int_size(
      proposals.size(), "NNI proposal count exceeds its integer domain");
  (void)checked_hat_nni_int_size(
      query.nodes.size(), "NNI query exceeds its integer node domain");
  (void)checked_hat_nni_int_size(
      active_taxa.size(), "NNI active-taxon set exceeds its integer domain");

  std::vector<int> proposal_by_q(query.nodes.size(), -1);
  for (int proposal_index = 0; proposal_index < proposal_count;
       ++proposal_index) {
    const int q = proposals[static_cast<std::size_t>(proposal_index)].q_query;
    if (q >= 0 && q < static_cast<int>(proposal_by_q.size())) {
      proposal_by_q[static_cast<std::size_t>(q)] = proposal_index;
    }
  }

  std::vector<int> query_active_to_owner(query.leaf_count, -1);
  for (int active = 0; active < query.leaf_count; ++active) {
    if (active < static_cast<int>(query.active_to_node.size()) &&
        query.active_to_node[static_cast<std::size_t>(active)] >= 0 &&
        active < static_cast<int>(active_taxa.size())) {
      query_active_to_owner[static_cast<std::size_t>(active)] =
          active_taxa[static_cast<std::size_t>(active)];
    }
  }
  const std::vector<int> leaf_owner_order =
      query_leaf_owner_order(query, query_active_to_owner);
  std::vector<HatNniOwnerRange> q_split_ranges(proposals.size());
  long long qsplit_owner_slots = 0;
  long long max_qsplit_leaves = 0;
  for (int proposal_index = 0; proposal_index < proposal_count;
       ++proposal_index) {
    const int q_split =
        proposals[static_cast<std::size_t>(proposal_index)].q_split_query;
    if (q_split < 0 ||
        q_split >= static_cast<int>(query.subtree_leaf_first.size()) ||
        query.leaf_active_order.empty()) {
      continue;
    }
    const int raw_begin =
        query.subtree_leaf_first[static_cast<std::size_t>(q_split)];
    const int raw_end =
        raw_begin + query.subtree_leaf_count[static_cast<std::size_t>(q_split)];
    const int begin =
        std::clamp(raw_begin, 0,
                   static_cast<int>(leaf_owner_order.size()));
    const int end =
        std::clamp(raw_end, begin,
                   static_cast<int>(leaf_owner_order.size()));
    q_split_ranges[static_cast<std::size_t>(proposal_index)] =
        HatNniOwnerRange{begin, end};
    const long long leaves = std::max(0, end - begin);
    qsplit_owner_slots += leaves;
    max_qsplit_leaves = std::max(max_qsplit_leaves, leaves);
  }
  benchmark_counter_max(benchmark_stats().hat_nni_peak_proposals,
                        static_cast<long long>(proposals.size()));
  benchmark_counter_max(benchmark_stats().hat_nni_max_qsplit_leaves,
                        max_qsplit_leaves);
  benchmark_counter_add(benchmark_stats().hat_nni_qsplit_owner_bytes_avoided,
                        qsplit_owner_slots *
                            static_cast<long long>(sizeof(int)));

  OwnerProjection owners;
  owners.taxon_to_owner.resize(static_cast<std::size_t>(taxon_count));
  std::iota(owners.taxon_to_owner.begin(), owners.taxon_to_owner.end(), 0);
  owners.owner_count = taxon_count;
  owners.live_owner.assign(static_cast<std::size_t>(taxon_count), 0);
  owners.live_taxa.reserve(active_taxa.size());
  for (int taxon : active_taxa) {
    if (taxon >= 0 && taxon < taxon_count) {
      owners.live_owner[static_cast<std::size_t>(taxon)] = 1;
      owners.live_taxa.push_back(taxon);
    }
  }

  const std::vector<GeneIndex>& genes = scoring_genes(ctx);
  const int gene_count = checked_hat_nni_int_size(
      genes.size(), "NNI gene count exceeds its integer domain");
  const int chunks = hat_score_gene_partition_count(gene_count);
  const std::size_t partial_count = checked_hat_nni_size_product(
      static_cast<std::size_t>(chunks), proposals.size(),
      "NNI partial score table exceeds the addressable size domain");
  std::vector<HatNniProposalScores> partial(partial_count);
  const std::size_t partial_bytes = checked_hat_nni_size_product(
      partial_count, sizeof(HatNniProposalScores),
      "NNI partial score bytes exceed the addressable size domain");
  if (partial_bytes >
      static_cast<std::size_t>(std::numeric_limits<long long>::max())) {
    throw std::length_error(
        "NNI partial score bytes exceed the telemetry domain");
  }
  benchmark_counter_max(benchmark_stats().hat_score_partial_peak_bytes,
                        static_cast<long long>(partial_bytes));
  std::vector<long long> local_queries(static_cast<std::size_t>(chunks), 0);
  std::vector<long long> local_recolors(static_cast<std::size_t>(chunks), 0);
  std::vector<long long> local_recomputed(static_cast<std::size_t>(chunks), 0);
  std::vector<HdtTraversalStats> local_hdt(
      static_cast<std::size_t>(chunks));

  parallel_for_inner(0, chunks, [&](int chunk) {
    const int begin = chunk;
    const int end = begin + 1;
    HatNniProposalScores* row =
        partial.data() + static_cast<std::size_t>(chunk) * proposals.size();
    long long queries = 0;
    HdtTraversalStats hdt_stats;
    HdtExtractWorkspace hdt_workspace;
    // NNI can score large root-scale queries after normalized splits. Scope the
    // fixed-vector workspace to this chunk so root-scale buffers do not persist in
    // Taskflow worker thread-locals across passes.
    for (int g = begin; g < end; ++g) {
      const GeneIndex& gene = genes[static_cast<std::size_t>(g)];
      HdtColorPolicy colors;
      colors.default_current = 0;
      colors.default_background = 0;
      colors.current_by_owner.assign(
          static_cast<std::size_t>(owners.owner_count), 0);
      colors.background_by_owner.assign(
          static_cast<std::size_t>(owners.owner_count), 0);
      for (int owner : leaf_owner_order) {
        if (owner < 0 || owner >= owners.owner_count) continue;
        colors.current_by_owner[static_cast<std::size_t>(owner)] = 1;
      }
      HdtInitialStateResult initial = hdt_build_initial_state(
          build_base_component_graph(gene, owners, 4, colors));
      hdt_account_initial_state(initial, hdt_stats);
      hdt_traverse_one_x_medians(
          std::move(initial.state), query, leaf_owner_order,
          hdt_workspace, hdt_stats,
          [&](int median, const std::array<int, 3>& branch_nodes,
              const std::array<std::uint8_t, 3>&,
              HdtState& current) {
            for (int branch_node : branch_nodes) {
              if (branch_node < 0 ||
                  branch_node >= static_cast<int>(proposal_by_q.size()) ||
                  query.parent[static_cast<std::size_t>(branch_node)] !=
                      median) {
                continue;
              }
              const int proposal_index =
                  proposal_by_q[static_cast<std::size_t>(branch_node)];
              if (proposal_index < 0) continue;
              const HatNniEdgeProposal& proposal =
                  proposals[static_cast<std::size_t>(proposal_index)];
              const HatNniOwnerRange range =
                  q_split_ranges[static_cast<std::size_t>(proposal_index)];
              hdt_recolor_owner_range(
                  current, leaf_owner_order, range.begin, range.end, 3);
              const std::array<double, 3> values =
                  hdt_root_counts(current);
              const HatFourColorCounts counts{values[0], values[1], values[2]};
              HatNniProposalScores& out =
                  row[static_cast<std::size_t>(proposal_index)];
              const int p0 = proposal.p_colors[0];
              const int p1 = proposal.p_colors[1];
              const int q0 = proposal.q_colors[0];
              const int q1 = proposal.q_colors[1];
              out.current += pair_pair_score(counts, p0, p1, q0, q1);
              out.swap_p0_q0 += pair_pair_score(counts, p1, q0, p0, q1);
              out.swap_p0_q1 += pair_pair_score(counts, p1, q1, p0, q0);
              hdt_recolor_owner_range(
                  current, leaf_owner_order, range.begin, range.end,
                  static_cast<std::uint8_t>(proposal.q_colors[0]));
              ++queries;
            }
          });
    }
    local_queries[static_cast<std::size_t>(chunk)] = queries;
    local_recolors[static_cast<std::size_t>(chunk)] =
        hdt_stats.recolor_changed_leaves;
    local_recomputed[static_cast<std::size_t>(chunk)] =
        hdt_stats.recolor_marked_nodes;
    local_hdt[static_cast<std::size_t>(chunk)] = hdt_stats;
  }, 1);

  for (int chunk = 0; chunk < chunks; ++chunk) {
    const HatNniProposalScores* row =
        partial.data() + static_cast<std::size_t>(chunk) * proposals.size();
    for (int proposal_index = 0;
         proposal_index < static_cast<int>(proposals.size());
         ++proposal_index) {
      const HatNniProposalScores& in =
          row[static_cast<std::size_t>(proposal_index)];
      HatNniProposalScores& out =
          totals[static_cast<std::size_t>(proposal_index)];
      out.current += in.current;
      out.swap_p0_q0 += in.swap_p0_q0;
      out.swap_p0_q1 += in.swap_p0_q1;
    }
    benchmark_counter_add(benchmark_stats().hat_hdt_queries,
                          local_queries[static_cast<std::size_t>(chunk)]);
    benchmark_counter_add(benchmark_stats().hat_hdt_recolor_ops,
                          local_recolors[static_cast<std::size_t>(chunk)]);
    benchmark_counter_add(benchmark_stats().hat_hdt_nodes_recomputed,
                          local_recomputed[static_cast<std::size_t>(chunk)]);
    const HdtTraversalStats& traversal =
        local_hdt[static_cast<std::size_t>(chunk)];
    benchmark_merge_hdt_traversal_stats(traversal);
  }
  benchmark_counter_add(benchmark_stats().hat_hdt_banks_built,
                        static_cast<long long>(gene_count));
  benchmark_counter_add(benchmark_stats().attachment_four_color_passes,
                        static_cast<long long>(proposals.size()));
  return totals;
}

bool apply_nni_move(SolverTree& tree, const HatNniMove& move) {
  if (!hat_has_edge(tree, move.edge_a, move.edge_b) ||
      !hat_has_edge(tree, move.edge_a, move.swap_a_neighbor) ||
      !hat_has_edge(tree, move.edge_b, move.swap_b_neighbor)) {
    return false;
  }
  hat_remove_edge(tree, move.edge_a, move.swap_a_neighbor);
  hat_remove_edge(tree, move.edge_b, move.swap_b_neighbor);
  hat_add_edge(tree, move.edge_a, move.swap_b_neighbor);
  hat_add_edge(tree, move.edge_b, move.swap_a_neighbor);
  return true;
}

HatScoredNniData score_hat_nni_edges(
    const SolverContext& ctx,
    const SolverTree& tree,
    const std::vector<int>& active_taxa,
    const std::vector<uint64_t>* allowed_central_edges = nullptr) {
  const int taxon_count = static_cast<int>(context_taxa(ctx).size());
  const HatQueryTree query =
      make_hat_query_tree_from_solver_tree(tree, active_taxa, taxon_count);
  HatScoredNniData data;
  data.proposals = collect_nni_proposals(tree, query, allowed_central_edges);
  const long long scored_moves = checked_hat_nni_doubled_counter(
      data.proposals.size(), "NNI proposal counter exceeds its integer domain");
  benchmark_counter_add(benchmark_stats().hat_nni_edges_scored,
                        scored_moves);
  if (allowed_central_edges == nullptr) {
    benchmark_counter_add(benchmark_stats().hat_nni_full_score_passes);
  } else {
    benchmark_counter_add(benchmark_stats().hat_nni_local_score_passes);
    benchmark_counter_add(benchmark_stats().hat_nni_local_edges_scored,
                          scored_moves);
  }
  data.scores =
      score_nni_proposals_hdt(ctx, query, data.proposals, active_taxa,
                              taxon_count);
  return data;
}

struct HatCachedNniMove {
  HatNniMove move;
  uint64_t central_edge_key = 0;
};

void add_cached_positive_moves(const HatScoredNniData& scored,
                               std::vector<HatCachedNniMove>& cache,
                               double min_gain) {
  for (int i = 0; i < static_cast<int>(scored.proposals.size()); ++i) {
    const HatNniEdgeProposal& proposal =
        scored.proposals[static_cast<std::size_t>(i)];
    const HatNniProposalScores& scores =
        scored.scores[static_cast<std::size_t>(i)];
    const uint64_t central_key =
        hat_nni_node_edge_key(proposal.swap_p0_q0.edge_a,
                              proposal.swap_p0_q0.edge_b);
    const double gain0 = scores.swap_p0_q0 - scores.current;
    if (hat_nni_gain_meets_threshold(gain0, min_gain)) {
      HatNniMove move = proposal.swap_p0_q0;
      move.gain = gain0;
      cache.push_back(HatCachedNniMove{move, central_key});
    } else if (gain0 > kHatNniGainEps) {
      benchmark_counter_add(
          benchmark_stats().hat_nni_moves_rejected_by_threshold);
    }
    const double gain1 = scores.swap_p0_q1 - scores.current;
    if (hat_nni_gain_meets_threshold(gain1, min_gain)) {
      HatNniMove move = proposal.swap_p0_q1;
      move.gain = gain1;
      cache.push_back(HatCachedNniMove{move, central_key});
    } else if (gain1 > kHatNniGainEps) {
      benchmark_counter_add(
          benchmark_stats().hat_nni_moves_rejected_by_threshold);
    }
  }
}

bool nni_move_applicable(const SolverTree& tree, const HatNniMove& move) {
  return hat_has_edge(tree, move.edge_a, move.edge_b) &&
         hat_has_edge(tree, move.edge_a, move.swap_a_neighbor) &&
         hat_has_edge(tree, move.edge_b, move.swap_b_neighbor);
}

std::array<int, 4> nni_move_stable_key(const HatNniMove& move) {
  return std::array<int, 4>{move.edge_a, move.edge_b, move.swap_a_neighbor,
                            move.swap_b_neighbor};
}

HatNniMove select_nni_moves_by_gain(
    const std::vector<HatNniMove>& candidates) {
  const auto selection = hat_select_score_maximum(
      candidates.begin(), candidates.end(), kHatNniGainEps,
      [](const HatNniMove& move) { return move.gain; },
      nni_move_stable_key);
  return selection.selected == candidates.end() ? HatNniMove{}
                                                 : *selection.selected;
}

HatNniMove best_cached_nni_move(const SolverTree& tree,
                                const std::vector<HatCachedNniMove>& cache) {
  const auto selection = hat_select_score_maximum(
      cache.begin(), cache.end(), kHatNniGainEps,
      [](const HatCachedNniMove& cached) { return cached.move.gain; },
      [](const HatCachedNniMove& cached) {
        return nni_move_stable_key(cached.move);
      },
      [&](const HatCachedNniMove& cached) {
        return cached.move.gain > kHatNniGainEps &&
               nni_move_applicable(tree, cached.move);
      });
  return selection.selected == cache.end() ? HatNniMove{}
                                            : selection.selected->move;
}

std::vector<int> nni_touched_nodes(const HatNniMove& move) {
  std::vector<int> seeds;
  seeds.reserve(4);
  const std::array<int, 4> nodes{move.edge_a, move.edge_b,
                                 move.swap_a_neighbor,
                                 move.swap_b_neighbor};
  for (int node : nodes) {
    if (node < 0) continue;
    if (std::find(seeds.begin(), seeds.end(), node) == seeds.end()) {
      seeds.push_back(node);
    }
  }
  return seeds;
}

bool unordered_edge_equal(int a, int b, int c, int d) {
  return (a == c && b == d) || (a == d && b == c);
}

template <typename NeighborFn>
std::vector<int> component_taxa_with_neighbors(const SolverTree& tree,
                                               int start,
                                               int blocked,
                                               int taxon_count,
                                               NeighborFn&& neighbors) {
  std::vector<int> taxa;
  if (start < 0 || start >= static_cast<int>(tree.nodes.size())) return taxa;
  std::vector<uint8_t> seen(tree.nodes.size(), 0);
  std::vector<int> stack;
  std::vector<int> nbs;
  stack.push_back(start);
  if (blocked >= 0 && blocked < static_cast<int>(seen.size())) {
    seen[static_cast<std::size_t>(blocked)] = 1;
  }
  while (!stack.empty()) {
    const int node = stack.back();
    stack.pop_back();
    if (node < 0 || node >= static_cast<int>(tree.nodes.size()) ||
        seen[static_cast<std::size_t>(node)] ||
        !tree.nodes[static_cast<std::size_t>(node)].active) {
      continue;
    }
    seen[static_cast<std::size_t>(node)] = 1;
    const int64_t id = tree.nodes[static_cast<std::size_t>(node)].identity_id;
    if (id >= 0 && id < taxon_count) {
      taxa.push_back(static_cast<int>(id));
    }
    nbs.clear();
    nbs.reserve(tree.nodes[static_cast<std::size_t>(node)].neighbors.size() + 2);
    neighbors(node, nbs);
    for (int nb : nbs) {
      if (nb >= 0 && nb < static_cast<int>(seen.size()) &&
          !seen[static_cast<std::size_t>(nb)]) {
        stack.push_back(nb);
      }
    }
  }
  return taxa;
}

bool split_key_from_taxa(const HatNniSplitSupportCache& support,
                         const std::vector<int>& taxa,
                         HatNniSplitKey& key) {
  if (!support.enabled || support.taxon_count < 4) return false;
  const int side_size = static_cast<int>(taxa.size());
  const int other_size = support.taxon_count - side_size;
  if (side_size <= 1 || other_size <= 1) return false;
  std::array<uint64_t, 4> sum{};
  bool contains_ref = false;
  for (int taxon : taxa) {
    if (taxon < 0 || taxon >= support.taxon_count) continue;
    contains_ref = contains_ref || taxon == 0;
    for (int word = 0; word < 4; ++word) {
      sum[static_cast<std::size_t>(word)] +=
          support.fingerprints.word[static_cast<std::size_t>(taxon)]
                                   [static_cast<std::size_t>(word)];
    }
  }
  key.size = static_cast<uint32_t>(side_size);
  key.h = sum;
  if (contains_ref) {
    key.size = static_cast<uint32_t>(other_size);
    key.h = hat_nni_sub_words(support.fingerprints.all, sum);
  } else if (side_size > other_size) {
    key.size = static_cast<uint32_t>(other_size);
    key.h = hat_nni_sub_words(support.fingerprints.all, sum);
  }
  return key.size >= 2 &&
         key.size <= static_cast<uint32_t>(support.taxon_count - 2);
}

int split_support_count(const HatNniSplitSupportCache& support,
                        const HatNniSplitKey& key) {
  const auto it = std::lower_bound(
      support.splits.begin(), support.splits.end(), key,
      [](const SplitAggregate256& entry, const HatNniSplitKey& sought) {
        return split_key_less(entry.key, sought);
      });
  return it == support.splits.end() || !(it->key == key)
             ? 0
             : static_cast<int>(it->count);
}

int split_support_delta_from_keys(
    const HatNniSplitSupportCache& support,
    const std::vector<HatNniSplitKey>& old_keys,
    const std::vector<HatNniSplitKey>& new_keys) {
  if (!support.enabled) return 0;
  std::vector<HatNniSplitKey> old_unique = old_keys;
  std::vector<HatNniSplitKey> new_unique = new_keys;
  const auto less = [](const HatNniSplitKey& left,
                       const HatNniSplitKey& right) {
    return split_key_less(left, right);
  };
  const auto deduplicate = [&](std::vector<HatNniSplitKey>& keys) {
    std::sort(keys.begin(), keys.end(), less);
    keys.erase(std::unique(keys.begin(), keys.end()), keys.end());
  };
  deduplicate(old_unique);
  deduplicate(new_unique);
  int delta = 0;
  for (const HatNniSplitKey& key : new_unique) {
    if (!std::binary_search(old_unique.begin(), old_unique.end(), key, less)) {
      delta += split_support_count(support, key);
    }
  }
  for (const HatNniSplitKey& key : old_unique) {
    if (!std::binary_search(new_unique.begin(), new_unique.end(), key, less)) {
      delta -= split_support_count(support, key);
    }
  }
  return delta;
}

bool old_central_split_key(const HatNniSplitSupportCache& support,
                           const SolverTree& tree,
                           int left,
                           int right,
                           HatNniSplitKey& key) {
  const auto neighbors = [&](int node, std::vector<int>& out) {
    for (int nb : tree.nodes[static_cast<std::size_t>(node)].neighbors) {
      out.push_back(nb);
    }
  };
  return split_key_from_taxa(
      support,
      component_taxa_with_neighbors(tree, right, left, support.taxon_count,
                                    neighbors),
      key);
}

int split_support_delta_for_move(const HatNniSplitSupportCache& support,
                                 const SolverTree& tree,
                                 const HatNniMove& move) {
  if (!support.enabled) return 0;
  HatNniSplitKey old_key;
  if (!old_central_split_key(support, tree, move.edge_a, move.edge_b,
                             old_key)) {
    return 0;
  }
  const auto neighbors = [&](int node, std::vector<int>& out) {
    for (int nb : tree.nodes[static_cast<std::size_t>(node)].neighbors) {
      if (unordered_edge_equal(node, nb, move.edge_a, move.swap_a_neighbor) ||
          unordered_edge_equal(node, nb, move.edge_b, move.swap_b_neighbor)) {
        continue;
      }
      out.push_back(nb);
    }
    if (node == move.edge_a) {
      out.push_back(move.swap_b_neighbor);
    } else if (node == move.swap_b_neighbor) {
      out.push_back(move.edge_a);
    }
    if (node == move.edge_b) {
      out.push_back(move.swap_a_neighbor);
    } else if (node == move.swap_a_neighbor) {
      out.push_back(move.edge_b);
    }
  };
  HatNniSplitKey new_key;
  if (!split_key_from_taxa(
          support,
          component_taxa_with_neighbors(tree, move.edge_b, move.edge_a,
                                        support.taxon_count, neighbors),
          new_key)) {
    return 0;
  }
  return split_support_delta_from_keys(support, {old_key}, {new_key});
}

std::vector<uint64_t> collect_affected_nni_edge_keys(
    const SolverTree& tree,
    const std::vector<int>& seeds,
    int radius) {
  std::vector<uint64_t> affected;
  if (tree.nodes.empty() || seeds.empty()) return affected;

  std::vector<int> distance(tree.nodes.size(), -1);
  std::vector<int> queue;
  queue.reserve(tree.nodes.size());
  for (int seed : seeds) {
    if (seed < 0 || seed >= static_cast<int>(tree.nodes.size()) ||
        !tree.nodes[static_cast<std::size_t>(seed)].active ||
        distance[static_cast<std::size_t>(seed)] >= 0) {
      continue;
    }
    distance[static_cast<std::size_t>(seed)] = 0;
    queue.push_back(seed);
  }

  for (std::size_t head = 0; head < queue.size(); ++head) {
    const int node = queue[head];
    const int node_distance = distance[static_cast<std::size_t>(node)];
    for (int nb : hat_active_neighbors(tree, node)) {
      affected.push_back(hat_nni_node_edge_key(node, nb));
      if (node_distance < radius &&
          distance[static_cast<std::size_t>(nb)] < 0) {
        distance[static_cast<std::size_t>(nb)] = node_distance + 1;
        queue.push_back(nb);
      }
    }
  }
  std::sort(affected.begin(), affected.end());
  affected.erase(std::unique(affected.begin(), affected.end()), affected.end());
  return affected;
}

void erase_cached_nni_edges(std::vector<HatCachedNniMove>& cache,
                            const SolverTree& tree,
                            const std::vector<uint64_t>& affected) {
  cache.erase(
      std::remove_if(
          cache.begin(), cache.end(),
          [&](const HatCachedNniMove& cached) {
            return std::binary_search(affected.begin(), affected.end(),
                                      cached.central_edge_key) ||
                   !nni_move_applicable(tree, cached.move);
          }),
      cache.end());
}

void merge_edge_key_sets(std::vector<uint64_t>& into,
                         const std::vector<uint64_t>& from) {
  std::vector<uint64_t> merged;
  merged.reserve(into.size() + from.size());
  std::set_union(into.begin(), into.end(), from.begin(), from.end(),
                 std::back_inserter(merged));
  into = std::move(merged);
}

}  // namespace

void hat_select_nni_paths_for_tests(
    const SolverTree& tree,
    const std::vector<HatNniMove>& candidates,
    double min_gain,
    HatNniMove& cached_result,
    HatNniMove& full_rescore_result) {
  std::vector<HatCachedNniMove> cache;
  std::vector<HatNniMove> full_candidates;
  cache.reserve(candidates.size());
  full_candidates.reserve(candidates.size());
  for (const HatNniMove& move : candidates) {
    if (!hat_nni_gain_meets_threshold(move.gain, min_gain)) continue;
    cache.push_back(HatCachedNniMove{
        move, hat_nni_node_edge_key(move.edge_a, move.edge_b)});
    full_candidates.push_back(move);
  }
  cached_result = best_cached_nni_move(tree, cache);
  full_rescore_result = select_nni_moves_by_gain(full_candidates);
}

HatNniMove hat_best_nni_move_hdt(const SolverContext& ctx,
                                 const SolverTree& tree,
                                 const std::vector<int>& active_taxa) {
  const double min_gain = hat_nni_min_gain_threshold(ctx);
  const HatScoredNniData scored = score_hat_nni_edges(ctx, tree, active_taxa);
  std::vector<HatNniMove> candidates;
  candidates.reserve(checked_hat_nni_doubled_capacity(
      scored.proposals.size(), false,
      "NNI candidate set exceeds the addressable size domain"));
  auto consider = [&](HatNniMove move, double gain) {
    if (!hat_nni_gain_meets_threshold(gain, min_gain)) {
      if (gain > kHatNniGainEps) {
        benchmark_counter_add(
            benchmark_stats().hat_nni_moves_rejected_by_threshold);
      }
      return;
    }
    move.gain = gain;
    candidates.push_back(move);
  };

  for (int i = 0; i < static_cast<int>(scored.proposals.size()); ++i) {
    const HatNniEdgeProposal& proposal =
        scored.proposals[static_cast<std::size_t>(i)];
    const HatNniProposalScores& scores =
        scored.scores[static_cast<std::size_t>(i)];
    consider(proposal.swap_p0_q0, scores.swap_p0_q0 - scores.current);
    consider(proposal.swap_p0_q1, scores.swap_p0_q1 - scores.current);
  }
  return select_nni_moves_by_gain(candidates);
}

void run_hat_nni_repair(SolverContext& ctx,
                        SolverTree& tree,
                        const std::vector<int>& active_taxa) {
  const bool unbounded = ctx.options.hat_nni_unbounded;
  const int pass_cap = std::max(0, ctx.options.hat_nni_passes);
  if (!unbounded && pass_cap <= 0) return;
  const double min_gain = hat_nni_min_gain_threshold(ctx);
  benchmark_stats().hat_nni_gain_threshold.store(
      min_gain, std::memory_order_relaxed);
  const auto started = std::chrono::steady_clock::now();
  progress_note_stage("hat-nni");
  const HatNniSplitSupportCache split_support =
      build_hat_nni_split_support_cache(ctx);

  const auto progress_edges = []() {
    return benchmark_stats().hat_nni_edges_scored.load(
        std::memory_order_relaxed);
  };

  if (ctx.options.hat_nni_local_rescore) {
    progress_note_hat_nni("v1-local-init", 0, pass_cap, unbounded,
                          progress_edges(), 0, 0.0, min_gain);
    HatNniPassRecorder recorder(ctx, "v1-local");
    recorder.record(tree, 0, 0, 0, 0.0, 0.0, 0, 0, progress_edges(),
                    min_gain);
    HatScoredNniData scored = score_hat_nni_edges(ctx, tree, active_taxa);
    std::vector<HatCachedNniMove> cache;
    cache.reserve(checked_hat_nni_doubled_capacity(
        scored.proposals.size(), false,
        "NNI move cache exceeds the addressable size domain"));
    add_cached_positive_moves(scored, cache, min_gain);

    int pass = 0;
    long long accepted_moves = 0;
    long long accepted_split_delta = 0;
    double accepted_gain = 0.0;
    HatNniPrefixGate gate = make_hat_nni_prefix_gate(ctx);
    while (unbounded || pass < pass_cap) {
      HatNniMove move = best_cached_nni_move(tree, cache);
      if (move.edge_a < 0 ||
          !hat_nni_gain_meets_threshold(move.gain, min_gain)) {
        break;
      }
      const int split_delta =
          split_support_delta_for_move(split_support, tree, move);
      const HatNniGateDecision decision =
          hat_nni_prefix_gate_decide(ctx, gate, pass, move.gain, split_delta);
      if (decision != HatNniGateDecision::Allow) {
        hat_nni_prefix_gate_note_stop(decision);
        break;
      }
      const std::vector<int> seeds = nni_touched_nodes(move);
      std::vector<uint64_t> affected =
          collect_affected_nni_edge_keys(tree, seeds, 2);
      if (!apply_nni_move(tree, move)) {
        const uint64_t failed_key =
            hat_nni_node_edge_key(move.edge_a, move.edge_b);
        const auto position =
            std::lower_bound(affected.begin(), affected.end(), failed_key);
        if (position == affected.end() || *position != failed_key) {
          affected.insert(position, failed_key);
        }
        erase_cached_nni_edges(cache, tree, affected);
        continue;
      }

      ++pass;
      ++accepted_moves;
      accepted_gain += move.gain;
      accepted_split_delta += split_delta;
      benchmark_counter_add(benchmark_stats().hat_nni_passes_executed);
      benchmark_counter_add(benchmark_stats().hat_nni_moves_accepted);
      hat_atomic_double_add(benchmark_stats().hat_nni_verified_gain, move.gain);
      hat_note_accepted_nni_gain(move.gain);
      hat_nni_prefix_gate_note_accepted(gate, move.gain, split_delta);
      const std::vector<uint64_t> after_affected =
          collect_affected_nni_edge_keys(tree, seeds, 2);
      merge_edge_key_sets(affected, after_affected);
      erase_cached_nni_edges(cache, tree, affected);
      if (unbounded || pass < pass_cap) {
        HatScoredNniData local_scored =
            score_hat_nni_edges(ctx, tree, active_taxa, &affected);
        add_cached_positive_moves(local_scored, cache, min_gain);
      }
      progress_note_hat_nni("v1-local", pass, pass_cap, unbounded,
                            progress_edges(), accepted_moves, accepted_gain,
                            min_gain);
      recorder.record(tree, pass, 1, accepted_moves, move.gain,
                      accepted_gain, split_delta, accepted_split_delta,
                      progress_edges(), min_gain);
    }
    progress_note_hat_nni("v1-local-done", pass, pass_cap, unbounded,
                          progress_edges(), accepted_moves, accepted_gain,
                          min_gain);
    benchmark_timer_add(benchmark_stats().hat_nni_ns,
                        std::chrono::duration_cast<std::chrono::nanoseconds>(
                            std::chrono::steady_clock::now() - started));
    return;
  }

  int pass = 0;
  long long accepted_moves = 0;
  long long accepted_split_delta = 0;
  double accepted_gain = 0.0;
  HatNniPrefixGate gate = make_hat_nni_prefix_gate(ctx);
  HatNniPassRecorder recorder(ctx, "v1-full");
  recorder.record(tree, 0, 0, 0, 0.0, 0.0, 0, 0, progress_edges(), min_gain);
  while (unbounded || pass < pass_cap) {
    progress_note_hat_nni("v1-full-score", pass, pass_cap, unbounded,
                          progress_edges(), accepted_moves, accepted_gain,
                          min_gain);
    HatNniMove move = hat_best_nni_move_hdt(ctx, tree, active_taxa);
    if (move.edge_a < 0 ||
        !hat_nni_gain_meets_threshold(move.gain, min_gain)) {
      break;
    }
    const int split_delta =
        split_support_delta_for_move(split_support, tree, move);
    const HatNniGateDecision decision =
        hat_nni_prefix_gate_decide(ctx, gate, pass, move.gain, split_delta);
    if (decision != HatNniGateDecision::Allow) {
      hat_nni_prefix_gate_note_stop(decision);
      break;
    }
    if (!apply_nni_move(tree, move)) break;
    ++pass;
    ++accepted_moves;
    accepted_gain += move.gain;
    accepted_split_delta += split_delta;
    benchmark_counter_add(benchmark_stats().hat_nni_passes_executed);
    benchmark_counter_add(benchmark_stats().hat_nni_moves_accepted);
    hat_atomic_double_add(benchmark_stats().hat_nni_verified_gain, move.gain);
    hat_note_accepted_nni_gain(move.gain);
    hat_nni_prefix_gate_note_accepted(gate, move.gain, split_delta);
    progress_note_hat_nni("v1-full", pass, pass_cap, unbounded,
                          progress_edges(), accepted_moves, accepted_gain,
                          min_gain);
    recorder.record(tree, pass, 1, accepted_moves, move.gain,
                    accepted_gain, split_delta, accepted_split_delta,
                    progress_edges(), min_gain);
  }
  progress_note_hat_nni("v1-full-done", pass, pass_cap, unbounded,
                        progress_edges(), accepted_moves, accepted_gain,
                        min_gain);
  benchmark_timer_add(benchmark_stats().hat_nni_ns,
                      std::chrono::duration_cast<std::chrono::nanoseconds>(
                          std::chrono::steady_clock::now() - started));
}

}  // namespace wqfm::internal
