#include "internal.hpp"
#include "candidates/split_records.hpp"

#include <array>
#include <climits>
#include <cstdlib>
#include <iterator>
#include <limits>
#include <numeric>
#include <queue>
#include <sstream>

namespace wqfm::internal {
namespace {

std::size_t checked_mre_size_product(std::size_t first,
                                     std::size_t second,
                                     const char* what) {
  if (first != 0 && second > std::numeric_limits<std::size_t>::max() / first) {
    throw std::length_error(what);
  }
  return first * second;
}

std::size_t checked_mre_flat_slots(int first,
                                   int second,
                                   const char* what) {
  if (first < 0 || second < 0) throw std::length_error(what);
  const std::size_t slots = checked_mre_size_product(
      static_cast<std::size_t>(first), static_cast<std::size_t>(second), what);
  if (slots > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
    throw std::length_error(what);
  }
  return slots;
}

void validate_mre_dimensions(const std::vector<GeneIndex>& genes,
                             int taxon_count) {
  if (taxon_count < 0 ||
      taxon_count > (std::numeric_limits<int>::max() - 1) / 2) {
    throw std::length_error("MRE taxon count exceeds its integer domain");
  }
  if (genes.size() >
      static_cast<std::size_t>(std::numeric_limits<int>::max())) {
    throw std::length_error("MRE gene count exceeds its integer domain");
  }
  const int max_nodes = std::max(2, 2 * taxon_count + 1);
  (void)checked_mre_flat_slots(
      max_nodes, static_cast<int>(genes.size()),
      "MRE flat state exceeds its integer index domain");
  constexpr std::size_t kMaxEulerNodes =
      (static_cast<std::size_t>(std::numeric_limits<int>::max()) + 1u) / 2u;
  for (const GeneIndex& gene : genes) {
    if (gene.nodes.size() > kMaxEulerNodes ||
        gene.taxon_to_node.size() >
            static_cast<std::size_t>(std::numeric_limits<int>::max())) {
      throw std::length_error("MRE gene index exceeds its integer domain");
    }
  }
}

using MreEntry = SplitAggregate256;

bool entry_less(const MreEntry& a, const MreEntry& b) {
  return split_priority_less_by_key(a, b);
}

std::vector<MreEntry> count_mre_entries_cpu_range(
    const std::vector<GeneIndex>& genes,
    int taxon_count,
    const SplitTaxonFingerprints256& fingerprints,
    int begin_gene,
    int end_gene) {
  begin_gene = std::max(0, std::min(begin_gene,
                                    static_cast<int>(genes.size())));
  end_gene = std::max(begin_gene,
                      std::min(end_gene, static_cast<int>(genes.size())));
  const int gene_count = end_gene - begin_gene;
  if (gene_count <= 0) return {};
  const int chunks =
      std::max(1, std::min(gene_count, current_inner_parallel_worker_limit()));
  std::vector<std::vector<SplitOccurrence256>> local_records(
      static_cast<std::size_t>(chunks));

  parallel_for_inner(0, chunks, [&](int chunk) {
    const int begin = begin_gene + static_cast<int>(
        static_cast<int64_t>(gene_count) * chunk / chunks);
    const int end = begin_gene + static_cast<int>(
        static_cast<int64_t>(gene_count) * (chunk + 1) / chunks);
    auto& local = local_records[static_cast<std::size_t>(chunk)];
    std::size_t capacity = 0;
    for (int g = begin; g < end; ++g) {
      const std::size_t nodes = genes[static_cast<std::size_t>(g)].nodes.size();
      if (nodes > std::numeric_limits<std::size_t>::max() - capacity) {
        throw std::length_error(
            "MRE split records exceed the addressable size domain");
      }
      capacity += nodes;
    }
    local.reserve(capacity);
    std::vector<uint32_t> subtree_size;
    std::vector<uint8_t> contains_ref;
    std::vector<std::array<uint64_t, 4>> subtree_sum;
    std::vector<std::array<uint64_t, 4>> duplicate_balance;
    std::vector<uint32_t> one_taxon;

    for (int g = begin; g < end; ++g) {
      const GeneIndex& gene = genes[static_cast<std::size_t>(g)];
      if (gene.nodes.empty()) {
        progress_consensus_count_gene();
        continue;
      }
      const GeneTraversal& tr = gene.reference_traversal;
      const int node_count = static_cast<int>(gene.nodes.size());
      const GeneSplitLeafProjection projection =
          make_gene_split_leaf_projection(gene, taxon_count, fingerprints.word);
      subtree_size.assign(static_cast<std::size_t>(node_count), 0);
      contains_ref.assign(static_cast<std::size_t>(node_count), 0);
      subtree_sum.assign(static_cast<std::size_t>(node_count),
                         std::array<uint64_t, 4>{});
      duplicate_balance.assign(static_cast<std::size_t>(node_count),
                               std::array<uint64_t, 4>{});
      one_taxon.assign(static_cast<std::size_t>(node_count), UINT32_MAX);

      for (int node : tr.postorder) {
        subtree_size[static_cast<std::size_t>(node)] =
            projection.unique_size[static_cast<std::size_t>(node)];
        contains_ref[static_cast<std::size_t>(node)] =
            projection.reference_marker[static_cast<std::size_t>(node)];
        subtree_sum[static_cast<std::size_t>(node)] =
            projection.unique_fingerprint[static_cast<std::size_t>(node)];
        duplicate_balance[static_cast<std::size_t>(node)] =
            projection.duplicate_balance[static_cast<std::size_t>(node)];
        const int taxon = gene.nodes[static_cast<std::size_t>(node)].taxon;
        if (projection.unique_size[static_cast<std::size_t>(node)] != 0) {
          one_taxon[static_cast<std::size_t>(node)] =
              static_cast<uint32_t>(taxon);
        }
        for (int nb : gene.nodes[static_cast<std::size_t>(node)].neighbors) {
          if (tr.parent[static_cast<std::size_t>(nb)] != node) continue;
          subtree_size[static_cast<std::size_t>(node)] +=
              subtree_size[static_cast<std::size_t>(nb)];
          contains_ref[static_cast<std::size_t>(node)] =
              static_cast<uint8_t>(
                  contains_ref[static_cast<std::size_t>(node)] |
                  contains_ref[static_cast<std::size_t>(nb)]);
          for (int w = 0; w < 4; ++w) {
            subtree_sum[static_cast<std::size_t>(node)]
                       [static_cast<std::size_t>(w)] +=
                subtree_sum[static_cast<std::size_t>(nb)]
                           [static_cast<std::size_t>(w)];
            duplicate_balance[static_cast<std::size_t>(node)]
                             [static_cast<std::size_t>(w)] +=
                duplicate_balance[static_cast<std::size_t>(nb)]
                                 [static_cast<std::size_t>(w)];
          }
          if (one_taxon[static_cast<std::size_t>(node)] == UINT32_MAX) {
            one_taxon[static_cast<std::size_t>(node)] =
                one_taxon[static_cast<std::size_t>(nb)];
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
        SplitKey256 key;
        key.size = size;
        key.h = subtree_sum[static_cast<std::size_t>(child)];
        local.push_back(SplitOccurrence256{
            key, static_cast<uint32_t>(g), edge_id,
            static_cast<uint32_t>(child),
            one_taxon[static_cast<std::size_t>(child)]});
        ++edge_id;
      }
      progress_consensus_count_gene();
    }
  }, 1);

  std::size_t observed = 0;
  for (const auto& local : local_records) {
    if (local.size() > std::numeric_limits<std::size_t>::max() - observed) {
      throw std::length_error(
          "MRE split records exceed the addressable size domain");
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
  std::vector<MreEntry> entries =
      aggregate_split_occurrences256(occurrences, radix_stats);
  std::vector<SplitOccurrence256>().swap(occurrences);
  radix_sort_split_priority256(entries, radix_stats);
  BenchmarkStats& stats = benchmark_stats();
  stats.mre_radix_passes.fetch_add(radix_stats.total_passes(),
                                   std::memory_order_relaxed);
  stats.mre_radix_items.fetch_add(radix_stats.total_items(),
                                  std::memory_order_relaxed);
  return entries;
}

std::vector<MreEntry> count_mre_entries(
    const std::vector<GeneIndex>& genes,
    int taxon_count,
    const Options& options) {
  const SplitTaxonFingerprints256 fingerprints =
      build_split_taxon_fingerprints256(taxon_count, options);
  const int gene_count = static_cast<int>(genes.size());
  return count_mre_entries_cpu_range(genes, taxon_count, fingerprints, 0,
                                        gene_count);
}

struct StaticLcaIndex {
  std::vector<int> tin;
  std::vector<int> tout;
  std::vector<int> depth;
  std::vector<int> first_euler;
  std::vector<int> euler;
  std::vector<int> taxon_to_node;
  int block_size = 1;
  int block_count = 0;
  std::vector<int> block_signature;
  std::vector<int> block_min_pos;
  std::vector<std::vector<int>> sparse;
  std::vector<uint8_t> micro_tables;

  bool in_subtree(int node, int root) const {
    if (node < 0 || root < 0) return false;
    return tin[static_cast<std::size_t>(root)] <=
               tin[static_cast<std::size_t>(node)] &&
           tin[static_cast<std::size_t>(node)] <
               tout[static_cast<std::size_t>(root)];
  }

  int better(int a, int b) const {
    if (a < 0) return b;
    if (b < 0) return a;
    return depth[static_cast<std::size_t>(a)] <=
                   depth[static_cast<std::size_t>(b)]
               ? a
               : b;
  }

  int lca(int a, int b) const {
    if (a < 0) return b;
    if (b < 0) return a;
    int l = first_euler[static_cast<std::size_t>(a)];
    int r = first_euler[static_cast<std::size_t>(b)];
    if (l > r) std::swap(l, r);
    const int pos = rmq_pos(l, r);
    return euler[static_cast<std::size_t>(pos)];
  }

  int better_pos(int a, int b) const {
    if (a < 0) return b;
    if (b < 0) return a;
    const int na = euler[static_cast<std::size_t>(a)];
    const int nb = euler[static_cast<std::size_t>(b)];
    return depth[static_cast<std::size_t>(na)] <=
                   depth[static_cast<std::size_t>(nb)]
               ? a
               : b;
  }

  int micro_pos(int block, int lo, int hi) const {
    const int start = block * block_size;
    const int len = std::min(block_size,
                             static_cast<int>(euler.size()) - start);
    if (len <= 0) return -1;
    if (len < block_size) {
      int best = start + lo;
      for (int off = lo + 1; off <= hi && off < len; ++off) {
        best = better_pos(best, start + off);
      }
      return best;
    }
    const int sig = block_signature[static_cast<std::size_t>(block)];
    const int off = micro_tables[
        static_cast<std::size_t>((sig * block_size + lo) * block_size + hi)];
    return start + off;
  }

  int sparse_pos(int block_l, int block_r) const {
    if (block_l > block_r) return -1;
    const int len = block_r - block_l + 1;
    const int lg = 31 - __builtin_clz(static_cast<unsigned>(len));
    const int a = sparse[static_cast<std::size_t>(lg)]
                        [static_cast<std::size_t>(block_l)];
    const int b = sparse[static_cast<std::size_t>(lg)]
                        [static_cast<std::size_t>(block_r - (1 << lg) + 1)];
    return better_pos(a, b);
  }

  int rmq_pos(int l, int r) const {
    const int bl = l / block_size;
    const int br = r / block_size;
    const int ol = l - bl * block_size;
    const int orr = r - br * block_size;
    if (bl == br) return micro_pos(bl, ol, orr);
    int best = micro_pos(bl, ol, block_size - 1);
    best = better_pos(best, sparse_pos(bl + 1, br - 1));
    best = better_pos(best, micro_pos(br, 0, orr));
    return best;
  }
};

StaticLcaIndex build_static_lca_index(const GeneIndex& gene, int taxon_count) {
  StaticLcaIndex out;
  if (taxon_count < 0 ||
      gene.nodes.size() >
          (static_cast<std::size_t>(std::numeric_limits<int>::max()) + 1u) /
              2u ||
      gene.taxon_to_node.size() >
          static_cast<std::size_t>(std::numeric_limits<int>::max())) {
    throw std::length_error("MRE static LCA exceeds its integer domain");
  }
  const int node_count = static_cast<int>(gene.nodes.size());
  out.tin.assign(static_cast<std::size_t>(node_count), 0);
  out.tout.assign(static_cast<std::size_t>(node_count), 0);
  out.depth.assign(static_cast<std::size_t>(node_count), 0);
  out.first_euler.assign(static_cast<std::size_t>(node_count), -1);
  out.taxon_to_node = gene.taxon_to_node;
  if (static_cast<int>(out.taxon_to_node.size()) < taxon_count) {
    out.taxon_to_node.resize(static_cast<std::size_t>(taxon_count), -1);
  }
  if (node_count == 0) return out;
  if (gene.reference_root < 0 || gene.reference_root >= node_count) {
    throw std::invalid_argument("MRE gene has an invalid reference root");
  }

  int timer = 0;
  std::vector<std::tuple<int, int, int>> stack;
  stack.push_back({gene.reference_root, gene.reference_root, 0});
  while (!stack.empty()) {
    auto& [node, parent, next] = stack.back();
    if (next == 0) {
      out.tin[static_cast<std::size_t>(node)] = timer++;
      if (out.first_euler[static_cast<std::size_t>(node)] < 0) {
        out.first_euler[static_cast<std::size_t>(node)] =
            static_cast<int>(out.euler.size());
      }
      out.euler.push_back(node);
    }
    if (next < static_cast<int>(gene.nodes[static_cast<std::size_t>(node)]
                                    .neighbors.size())) {
      const int nb = gene.nodes[static_cast<std::size_t>(node)]
                         .neighbors[static_cast<std::size_t>(next++)];
      if (nb == parent) continue;
      out.depth[static_cast<std::size_t>(nb)] =
          out.depth[static_cast<std::size_t>(node)] + 1;
      stack.push_back({nb, node, 0});
      continue;
    }
    out.tout[static_cast<std::size_t>(node)] = timer;
    stack.pop_back();
    if (!stack.empty()) {
      out.euler.push_back(std::get<0>(stack.back()));
    }
  }

  const int euler_size = static_cast<int>(out.euler.size());
  const int lg = std::max(1, 31 - __builtin_clz(static_cast<unsigned>(
                                 std::max(1, euler_size))));
  out.block_size = std::max(2, lg / 2);
  out.block_count =
      euler_size == 0 ? 0 : 1 + (euler_size - 1) / out.block_size;
  const int signature_count = 1 << std::max(0, out.block_size - 1);
  const std::size_t micro_rows = checked_mre_size_product(
      static_cast<std::size_t>(signature_count),
      static_cast<std::size_t>(out.block_size),
      "MRE micro-RMQ table exceeds the addressable size domain");
  const std::size_t micro_slots = checked_mre_size_product(
      micro_rows, static_cast<std::size_t>(out.block_size),
      "MRE micro-RMQ table exceeds the addressable size domain");
  out.micro_tables.assign(micro_slots, 0);
  for (int sig = 0; sig < signature_count; ++sig) {
    std::vector<int> rel(static_cast<std::size_t>(out.block_size), 0);
    for (int i = 1; i < out.block_size; ++i) {
      rel[static_cast<std::size_t>(i)] =
          rel[static_cast<std::size_t>(i - 1)] +
          (((sig >> (i - 1)) & 1) != 0 ? 1 : -1);
    }
    for (int i = 0; i < out.block_size; ++i) {
      int best = i;
      for (int j = i; j < out.block_size; ++j) {
        if (rel[static_cast<std::size_t>(j)] <
            rel[static_cast<std::size_t>(best)]) {
          best = j;
        }
        out.micro_tables[static_cast<std::size_t>(
            (sig * out.block_size + i) * out.block_size + j)] =
            static_cast<uint8_t>(best);
      }
    }
  }

  out.block_signature.assign(static_cast<std::size_t>(out.block_count), 0);
  out.block_min_pos.assign(static_cast<std::size_t>(out.block_count), -1);
  for (int block = 0; block < out.block_count; ++block) {
    const int start = block * out.block_size;
    const int end = std::min(euler_size, start + out.block_size);
    int sig = 0;
    int best = start;
    for (int pos = start + 1; pos < end; ++pos) {
      const int prev_node = out.euler[static_cast<std::size_t>(pos - 1)];
      const int node = out.euler[static_cast<std::size_t>(pos)];
      if (out.depth[static_cast<std::size_t>(node)] >
          out.depth[static_cast<std::size_t>(prev_node)]) {
        sig |= 1 << (pos - start - 1);
      }
      best = out.better_pos(best, pos);
    }
    out.block_signature[static_cast<std::size_t>(block)] = sig;
    out.block_min_pos[static_cast<std::size_t>(block)] = best;
  }

  int levels = 1;
  while ((uint64_t{1} << levels) <=
         static_cast<uint64_t>(std::max(1, out.block_count))) {
    ++levels;
  }
  out.sparse.assign(static_cast<std::size_t>(levels),
                    std::vector<int>(static_cast<std::size_t>(out.block_count),
                                     -1));
  if (out.block_count > 0) {
    out.sparse[0] = out.block_min_pos;
    for (int level = 1; level < levels; ++level) {
      const int span = static_cast<int>(uint32_t{1} << level);
      const int half = span >> 1;
      for (int block = 0; block + span <= out.block_count; ++block) {
        out.sparse[static_cast<std::size_t>(level)]
                  [static_cast<std::size_t>(block)] =
            out.better_pos(out.sparse[static_cast<std::size_t>(level - 1)]
                                     [static_cast<std::size_t>(block)],
                           out.sparse[static_cast<std::size_t>(level - 1)]
                                     [static_cast<std::size_t>(block + half)]);
      }
      for (int block = std::max(0, out.block_count - span + 1);
           block < out.block_count; ++block) {
        if (out.sparse[static_cast<std::size_t>(level)]
                      [static_cast<std::size_t>(block)] < 0) {
          out.sparse[static_cast<std::size_t>(level)]
                    [static_cast<std::size_t>(block)] =
              out.sparse[static_cast<std::size_t>(level - 1)]
                        [static_cast<std::size_t>(block)];
        }
      }
    }
  }
  return out;
}

class LinkCutTree {
 public:
  explicit LinkCutTree(int n = 0) { reset(n); }

  void reset(int n) {
    if (n < 0) throw std::invalid_argument("negative Link-Cut tree size");
    ch_.assign(static_cast<std::size_t>(n), {-1, -1});
    parent_.assign(static_cast<std::size_t>(n), -1);
    rev_.assign(static_cast<std::size_t>(n), 0);
    size_.assign(static_cast<std::size_t>(n), 1);
    push_stack_.clear();
    push_stack_.reserve(static_cast<std::size_t>(n));
  }

  void reserve(int n) {
    if (n < 0) throw std::invalid_argument("negative Link-Cut tree capacity");
    const std::size_t capacity = static_cast<std::size_t>(n);
    ch_.reserve(capacity);
    parent_.reserve(capacity);
    rev_.reserve(capacity);
    size_.reserve(capacity);
    push_stack_.reserve(capacity);
  }

  void add_node() {
    ch_.push_back({-1, -1});
    parent_.push_back(-1);
    rev_.push_back(0);
    size_.push_back(1);
  }

  void link(int a, int b) {
    makeroot(a);
    if (findroot(b) != a) parent_[static_cast<std::size_t>(a)] = b;
  }

  void cut(int a, int b) {
    makeroot(a);
    access(b);
    splay(b);
    if (ch_[static_cast<std::size_t>(b)][0] == a) {
      ch_[static_cast<std::size_t>(b)][0] = -1;
      parent_[static_cast<std::size_t>(a)] = -1;
      pull(b);
    }
  }

  int path_size(int root, int node) {
    makeroot(root);
    access(node);
    splay(node);
    return size_[static_cast<std::size_t>(node)];
  }

  int kth_on_path(int root, int node, int k) {
    makeroot(root);
    access(node);
    splay(node);
    int cur = node;
    while (cur >= 0) {
      push(cur);
      const int left = ch_[static_cast<std::size_t>(cur)][0];
      const int left_size =
          left >= 0 ? size_[static_cast<std::size_t>(left)] : 0;
      if (k < left_size) {
        cur = left;
      } else if (k == left_size) {
        splay(cur);
        return cur;
      } else {
        k -= left_size + 1;
        cur = ch_[static_cast<std::size_t>(cur)][1];
      }
    }
    return node;
  }

 private:
  std::vector<std::array<int, 2>> ch_;
  std::vector<int> parent_;
  std::vector<uint8_t> rev_;
  std::vector<int> size_;
  std::vector<int> push_stack_;

  bool is_root(int x) const {
    const int p = parent_[static_cast<std::size_t>(x)];
    return p < 0 ||
           (ch_[static_cast<std::size_t>(p)][0] != x &&
            ch_[static_cast<std::size_t>(p)][1] != x);
  }

  void pull(int x) {
    int size = 1;
    for (int side = 0; side < 2; ++side) {
      const int c = ch_[static_cast<std::size_t>(x)][side];
      if (c >= 0) size += size_[static_cast<std::size_t>(c)];
    }
    size_[static_cast<std::size_t>(x)] = size;
  }

  void push(int x) {
    if (x < 0 || rev_[static_cast<std::size_t>(x)] == 0) return;
    std::swap(ch_[static_cast<std::size_t>(x)][0],
              ch_[static_cast<std::size_t>(x)][1]);
    for (int side = 0; side < 2; ++side) {
      const int c = ch_[static_cast<std::size_t>(x)][side];
      if (c >= 0) rev_[static_cast<std::size_t>(c)] ^= 1;
    }
    rev_[static_cast<std::size_t>(x)] = 0;
  }

  void push_path(int x) {
    push_stack_.clear();
    push_stack_.push_back(x);
    while (!is_root(x)) {
      x = parent_[static_cast<std::size_t>(x)];
      push_stack_.push_back(x);
    }
    for (auto it = push_stack_.rbegin(); it != push_stack_.rend(); ++it) {
      push(*it);
    }
  }

  void rotate(int x) {
    const int p = parent_[static_cast<std::size_t>(x)];
    const int g = parent_[static_cast<std::size_t>(p)];
    const int side = ch_[static_cast<std::size_t>(p)][1] == x ? 1 : 0;
    const int b = ch_[static_cast<std::size_t>(x)][side ^ 1];
    if (!is_root(p)) {
      ch_[static_cast<std::size_t>(g)]
         [ch_[static_cast<std::size_t>(g)][1] == p ? 1 : 0] = x;
    }
    parent_[static_cast<std::size_t>(x)] = g;
    ch_[static_cast<std::size_t>(x)][side ^ 1] = p;
    parent_[static_cast<std::size_t>(p)] = x;
    ch_[static_cast<std::size_t>(p)][side] = b;
    if (b >= 0) parent_[static_cast<std::size_t>(b)] = p;
    pull(p);
    pull(x);
  }

  void splay(int x) {
    push_path(x);
    while (!is_root(x)) {
      const int p = parent_[static_cast<std::size_t>(x)];
      const int g = parent_[static_cast<std::size_t>(p)];
      if (!is_root(p)) {
        const bool zigzig =
            (ch_[static_cast<std::size_t>(p)][0] == x) ==
            (ch_[static_cast<std::size_t>(g)][0] == p);
        rotate(zigzig ? p : x);
      }
      rotate(x);
    }
  }

  void access(int x) {
    int last = -1;
    for (int y = x; y >= 0; y = parent_[static_cast<std::size_t>(y)]) {
      splay(y);
      ch_[static_cast<std::size_t>(y)][1] = last;
      pull(y);
      last = y;
    }
    splay(x);
  }

  void makeroot(int x) {
    access(x);
    rev_[static_cast<std::size_t>(x)] ^= 1;
    push(x);
  }

  int findroot(int x) {
    access(x);
    while (true) {
      push(x);
      const int left = ch_[static_cast<std::size_t>(x)][0];
      if (left < 0) break;
      x = left;
    }
    splay(x);
    return x;
  }
};

struct ChildAgg {
  int count = 0;
  uint64_t leaf_sum = 0;
  int min_node = -1;
  int max_node = -1;
};

struct ChildKey {
  int tin = 0;
  int child = 0;
};

bool key_less(const ChildKey& a, const ChildKey& b) {
  return a.tin != b.tin ? a.tin < b.tin : a.child < b.child;
}

bool key_equal(const ChildKey& a, const ChildKey& b) {
  return a.tin == b.tin && a.child == b.child;
}

class ChildRangeForest {
 public:
  struct InitialChild {
    int child = -1;
    int key_tin = 0;
    uint32_t leaf_count = 0;
    int lca_node = -1;
  };

  ChildRangeForest(int max_nodes,
                   int genes,
                   const std::vector<StaticLcaIndex>& lca)
      : max_nodes_(max_nodes), genes_(genes), lca_(lca) {
    const std::size_t flat_slots = checked_mre_flat_slots(
        max_nodes_, genes_,
        "MRE child index exceeds its integer index domain");
    roots_.assign(flat_slots, -1);
    // Flat node pool partitioned into fixed per-gene segments
    // [gene*max_nodes_, (gene+1)*max_nodes_). Pre-sized and never reallocated,
    // so distinct genes own disjoint index ranges and can be updated in
    // parallel without data races. Each gene's live AVL-node count is bounded
    // by node_count <= max_nodes_, so a segment never overflows.
    nodes_.assign(flat_slots, Node{});
    free_.assign(static_cast<std::size_t>(genes_), {});
    next_.assign(static_cast<std::size_t>(genes_), 0);
    operations_.assign(static_cast<std::size_t>(genes_), 0);
    rotations_.assign(static_cast<std::size_t>(genes_), 0);
    max_height_.assign(static_cast<std::size_t>(genes_), 0);
    height_violations_.assign(static_cast<std::size_t>(genes_), 0);
    enumeration_visits_.assign(static_cast<std::size_t>(genes_), 0);
    for (int gene = 0; gene < genes_; ++gene) {
      next_[static_cast<std::size_t>(gene)] = gene * max_nodes_;
    }
  }

  void insert_child(int parent,
                    int gene,
                    int child,
                    int key_tin,
                    uint32_t leaf_count,
                    int lca_node) {
    const int node = allocate(gene);
    Node& n = nodes_[static_cast<std::size_t>(node)];
    n.left = n.right = -1;
    n.height = 1;
    n.key = ChildKey{key_tin, child};
    n.child = child;
    n.leaf_count = leaf_count;
    n.lca_node = lca_node;
    update(node, gene);
    int& root = roots_[root_index(parent, gene)];
    root = insert(root, node, gene);
    ++operations_[static_cast<std::size_t>(gene)];
    note_root(root, gene);
  }

  void build_children(int parent,
                      int gene,
                      std::vector<InitialChild>& children) {
    std::sort(children.begin(), children.end(),
              [](const InitialChild& a, const InitialChild& b) {
                if (a.key_tin != b.key_tin) return a.key_tin < b.key_tin;
                return a.child < b.child;
              });
    int& root = roots_[root_index(parent, gene)];
    if (root >= 0) {
      throw std::logic_error("MRE child index rebuilt over a live root");
    }
    std::vector<int> ordered;
    ordered.reserve(children.size());
    for (const InitialChild& child : children) {
      const int node = allocate(gene);
      Node& n = nodes_[static_cast<std::size_t>(node)];
      n.left = n.right = -1;
      n.height = 1;
      n.key = ChildKey{child.key_tin, child.child};
      n.child = child.child;
      n.leaf_count = child.leaf_count;
      n.lca_node = child.lca_node;
      n.agg = ChildAgg{1, child.leaf_count, child.lca_node, child.lca_node};
      ordered.push_back(node);
    }
    root = build_balanced(ordered, 0, static_cast<int>(ordered.size()), gene);
    ++operations_[static_cast<std::size_t>(gene)];
    note_root(root, gene);
  }

  void erase_child(int parent, int gene, int child, int key_tin) {
    int& root = roots_[root_index(parent, gene)];
    bool erased = false;
    root = erase(root, ChildKey{key_tin, child}, gene, erased);
    if (!erased) {
      throw std::logic_error("MRE child index erased a missing key");
    }
    ++operations_[static_cast<std::size_t>(gene)];
    note_root(root, gene);
  }

  ChildAgg range_query(int parent, int gene, int tin_lo, int tin_hi) {
    ++operations_[static_cast<std::size_t>(gene)];
    return range_aggregate(root(parent, gene), ChildKey{tin_lo, INT_MIN},
                           ChildKey{tin_hi, INT_MIN}, gene);
  }

  std::vector<int> enumerate_range(int parent,
                                   int gene,
                                   int tin_lo,
                                   int tin_hi) {
    ++operations_[static_cast<std::size_t>(gene)];
    std::vector<int> out;
    enumerate_range_impl(root(parent, gene), ChildKey{tin_lo, INT_MIN},
                         ChildKey{tin_hi, INT_MIN}, gene, out);
    return out;
  }

  std::vector<int> enumerate_complement(int parent,
                                        int gene,
                                        int tin_lo,
                                        int tin_hi) {
    ++operations_[static_cast<std::size_t>(gene)];
    std::vector<int> out;
    const int parent_root = root(parent, gene);
    enumerate_prefix_impl(parent_root, ChildKey{tin_lo, INT_MIN}, gene, out);
    enumerate_suffix_impl(parent_root, ChildKey{tin_hi, INT_MIN}, gene, out);
    return out;
  }

  int aggregate_lca(int parent, int gene) const {
    const int r = root(parent, gene);
    if (r < 0) return -1;
    const ChildAgg& agg = nodes_[static_cast<std::size_t>(r)].agg;
    return lca_[static_cast<std::size_t>(gene)].lca(agg.min_node,
                                                    agg.max_node);
  }

  long long operation_count() const {
    return std::accumulate(operations_.begin(), operations_.end(), 0LL);
  }

  long long rotation_count() const {
    return std::accumulate(rotations_.begin(), rotations_.end(), 0LL);
  }

  int maximum_height() const {
    return max_height_.empty()
               ? 0
               : *std::max_element(max_height_.begin(), max_height_.end());
  }

  long long height_violation_count() const {
    return std::accumulate(height_violations_.begin(),
                           height_violations_.end(), 0LL);
  }

  long long enumeration_visit_count() const {
    return std::accumulate(enumeration_visits_.begin(),
                           enumeration_visits_.end(), 0LL);
  }

  static int height_limit_for_count(int count) {
    if (count <= 0) return 0;
    uint64_t capacity = 1;
    int ceil_log = 0;
    const uint64_t required = static_cast<uint64_t>(count) + 1;
    while (capacity < required) {
      capacity <<= 1;
      ++ceil_log;
    }
    return 2 * ceil_log + 1;
  }

  bool validate_root_for_tests(int parent, int gene) const {
    return validate_tree(root(parent, gene), gene).valid;
  }

 private:
  struct Node {
    int left = -1;
    int right = -1;
    int height = 1;
    ChildKey key;
    int child = -1;
    uint32_t leaf_count = 0;
    int lca_node = -1;
    ChildAgg agg;
  };

  int max_nodes_ = 0;
  int genes_ = 0;
  const std::vector<StaticLcaIndex>& lca_;
  std::vector<int> roots_;
  std::vector<Node> nodes_;
  std::vector<std::vector<int>> free_;
  std::vector<int> next_;
  std::vector<long long> operations_;
  std::vector<long long> rotations_;
  std::vector<int> max_height_;
  std::vector<long long> height_violations_;
  std::vector<long long> enumeration_visits_;

  std::size_t root_index(int parent, int gene) const {
    // Parent-major: accept_move iterates genes for a fixed parent, so keeping a
    // parent's per-gene roots contiguous preserves serial cache locality (the
    // dominant cost). Per-gene AVL pools (nodes_) are what get
    // updated in parallel and are already in disjoint segments.
    return static_cast<std::size_t>(parent) * static_cast<std::size_t>(genes_) +
           static_cast<std::size_t>(gene);
  }

  int root(int parent, int gene) const {
    return roots_[root_index(parent, gene)];
  }

  int allocate(int gene) {
    std::vector<int>& fr = free_[static_cast<std::size_t>(gene)];
    if (!fr.empty()) {
      const int node = fr.back();
      fr.pop_back();
      return node;
    }
    const int node = next_[static_cast<std::size_t>(gene)]++;
    const int segment_end = (gene + 1) * max_nodes_;
    if (node >= segment_end) {
      throw std::logic_error("MRE child index exceeded its flat segment");
    }
    return node;
  }

  void release(int gene, int node) {
    if (node >= 0) free_[static_cast<std::size_t>(gene)].push_back(node);
  }

  int combine_min_tin_node(int gene, int a, int b) const {
    if (a < 0) return b;
    if (b < 0) return a;
    const StaticLcaIndex& lca = lca_[static_cast<std::size_t>(gene)];
    return lca.tin[static_cast<std::size_t>(a)] <=
                   lca.tin[static_cast<std::size_t>(b)]
               ? a
               : b;
  }

  int combine_max_tin_node(int gene, int a, int b) const {
    if (a < 0) return b;
    if (b < 0) return a;
    const StaticLcaIndex& lca = lca_[static_cast<std::size_t>(gene)];
    return lca.tin[static_cast<std::size_t>(a)] >=
                   lca.tin[static_cast<std::size_t>(b)]
               ? a
               : b;
  }

  void update(int node, int gene) {
    if (node < 0) return;
    Node& n = nodes_[static_cast<std::size_t>(node)];
    n.height = 1 + std::max(node_height(n.left), node_height(n.right));
    n.agg.count = 1;
    n.agg.leaf_sum = n.leaf_count;
    n.agg.min_node = n.lca_node;
    n.agg.max_node = n.lca_node;
    if (n.left >= 0) {
      const ChildAgg& left = nodes_[static_cast<std::size_t>(n.left)].agg;
      n.agg.count += left.count;
      n.agg.leaf_sum += left.leaf_sum;
      n.agg.min_node = combine_min_tin_node(gene, n.agg.min_node,
                                            left.min_node);
      n.agg.max_node = combine_max_tin_node(gene, n.agg.max_node,
                                            left.max_node);
    }
    if (n.right >= 0) {
      const ChildAgg& right = nodes_[static_cast<std::size_t>(n.right)].agg;
      n.agg.count += right.count;
      n.agg.leaf_sum += right.leaf_sum;
      n.agg.min_node = combine_min_tin_node(gene, n.agg.min_node,
                                            right.min_node);
      n.agg.max_node = combine_max_tin_node(gene, n.agg.max_node,
                                            right.max_node);
    }
  }

  ChildAgg aggregate(int node) const {
    return node >= 0 ? nodes_[static_cast<std::size_t>(node)].agg : ChildAgg{};
  }

  int node_height(int node) const {
    return node >= 0 ? nodes_[static_cast<std::size_t>(node)].height : 0;
  }

  int balance_factor(int node) const {
    if (node < 0) return 0;
    const Node& n = nodes_[static_cast<std::size_t>(node)];
    return node_height(n.left) - node_height(n.right);
  }

  ChildAgg combine_aggregates(const ChildAgg& a,
                              const ChildAgg& b,
                              int gene) const {
    ChildAgg out;
    out.count = a.count + b.count;
    out.leaf_sum = a.leaf_sum + b.leaf_sum;
    out.min_node = combine_min_tin_node(gene, a.min_node, b.min_node);
    out.max_node = combine_max_tin_node(gene, a.max_node, b.max_node);
    return out;
  }

  ChildAgg self_aggregate(int node) const {
    if (node < 0) return ChildAgg{};
    const Node& n = nodes_[static_cast<std::size_t>(node)];
    return ChildAgg{1, n.leaf_count, n.lca_node, n.lca_node};
  }

  int rotate_left(int root, int gene) {
    Node& r = nodes_[static_cast<std::size_t>(root)];
    const int next_root = r.right;
    if (next_root < 0) {
      throw std::logic_error("MRE invalid left AVL rotation");
    }
    Node& next = nodes_[static_cast<std::size_t>(next_root)];
    r.right = next.left;
    next.left = root;
    update(root, gene);
    update(next_root, gene);
    ++rotations_[static_cast<std::size_t>(gene)];
    return next_root;
  }

  int rotate_right(int root, int gene) {
    Node& r = nodes_[static_cast<std::size_t>(root)];
    const int next_root = r.left;
    if (next_root < 0) {
      throw std::logic_error("MRE invalid right AVL rotation");
    }
    Node& next = nodes_[static_cast<std::size_t>(next_root)];
    r.left = next.right;
    next.right = root;
    update(root, gene);
    update(next_root, gene);
    ++rotations_[static_cast<std::size_t>(gene)];
    return next_root;
  }

  int rebalance(int root, int gene) {
    if (root < 0) return -1;
    update(root, gene);
    const int balance = balance_factor(root);
    if (balance > 1) {
      Node& r = nodes_[static_cast<std::size_t>(root)];
      if (balance_factor(r.left) < 0) {
        r.left = rotate_left(r.left, gene);
      }
      root = rotate_right(root, gene);
    } else if (balance < -1) {
      Node& r = nodes_[static_cast<std::size_t>(root)];
      if (balance_factor(r.right) > 0) {
        r.right = rotate_right(r.right, gene);
      }
      root = rotate_left(root, gene);
    }
    if (std::abs(balance_factor(root)) > 1) {
      ++height_violations_[static_cast<std::size_t>(gene)];
      throw std::logic_error("MRE AVL rebalance invariant failed");
    }
    return root;
  }

  int insert(int root, int node, int gene) {
    if (root < 0) return node;
    Node& r = nodes_[static_cast<std::size_t>(root)];
    const ChildKey key = nodes_[static_cast<std::size_t>(node)].key;
    if (key_equal(key, r.key)) {
      throw std::logic_error("MRE child index inserted a duplicate key");
    }
    if (key_less(key, r.key)) {
      r.left = insert(r.left, node, gene);
    } else {
      r.right = insert(r.right, node, gene);
    }
    return rebalance(root, gene);
  }

  int build_balanced(const std::vector<int>& ordered,
                     int begin,
                     int end,
                     int gene) {
    if (begin >= end) return -1;
    const int middle = begin + (end - begin) / 2;
    const int root = ordered[static_cast<std::size_t>(middle)];
    Node& r = nodes_[static_cast<std::size_t>(root)];
    r.left = build_balanced(ordered, begin, middle, gene);
    r.right = build_balanced(ordered, middle + 1, end, gene);
    update(root, gene);
    return root;
  }

  int detach_min(int root, int& minimum, int gene) {
    Node& r = nodes_[static_cast<std::size_t>(root)];
    if (r.left < 0) {
      minimum = root;
      return r.right;
    }
    r.left = detach_min(r.left, minimum, gene);
    return rebalance(root, gene);
  }

  int erase(int root,
            const ChildKey& key,
            int gene,
            bool& erased) {
    if (root < 0) return -1;
    Node& r = nodes_[static_cast<std::size_t>(root)];
    if (key_equal(r.key, key)) {
      erased = true;
      const int left = r.left;
      const int right = r.right;
      if (left < 0 || right < 0) {
        const int replacement = left >= 0 ? left : right;
        release(gene, root);
        return replacement;
      }
      int successor = -1;
      const int remaining_right = detach_min(right, successor, gene);
      Node& replacement = nodes_[static_cast<std::size_t>(successor)];
      replacement.left = left;
      replacement.right = remaining_right;
      release(gene, root);
      return rebalance(successor, gene);
    }
    if (key_less(key, r.key)) {
      r.left = erase(r.left, key, gene, erased);
    } else {
      r.right = erase(r.right, key, gene, erased);
    }
    return erased ? rebalance(root, gene) : root;
  }

  ChildAgg suffix_aggregate(int node,
                            const ChildKey& lo,
                            int gene) const {
    if (node < 0) return ChildAgg{};
    const Node& n = nodes_[static_cast<std::size_t>(node)];
    if (key_less(n.key, lo)) {
      return suffix_aggregate(n.right, lo, gene);
    }
    ChildAgg out = aggregate(n.right);
    out = combine_aggregates(out, self_aggregate(node), gene);
    return combine_aggregates(out, suffix_aggregate(n.left, lo, gene),
                              gene);
  }

  ChildAgg prefix_aggregate(int node,
                            const ChildKey& hi,
                            int gene) const {
    if (node < 0) return ChildAgg{};
    const Node& n = nodes_[static_cast<std::size_t>(node)];
    if (!key_less(n.key, hi)) {
      return prefix_aggregate(n.left, hi, gene);
    }
    ChildAgg out = aggregate(n.left);
    out = combine_aggregates(out, self_aggregate(node), gene);
    return combine_aggregates(out, prefix_aggregate(n.right, hi, gene),
                              gene);
  }

  ChildAgg range_aggregate(int root,
                           const ChildKey& lo,
                           const ChildKey& hi,
                           int gene) const {
    int split_node = root;
    while (split_node >= 0) {
      const Node& n = nodes_[static_cast<std::size_t>(split_node)];
      if (key_less(n.key, lo)) {
        split_node = n.right;
      } else if (!key_less(n.key, hi)) {
        split_node = n.left;
      } else {
        ChildAgg out = suffix_aggregate(n.left, lo, gene);
        out = combine_aggregates(out, self_aggregate(split_node), gene);
        return combine_aggregates(out,
                                  prefix_aggregate(n.right, hi, gene), gene);
      }
    }
    return ChildAgg{};
  }

  struct ValidationResult {
    bool valid = true;
    bool empty = true;
    ChildKey min_key;
    ChildKey max_key;
    int height = 0;
    ChildAgg agg;
  };

  static bool aggregate_equal(const ChildAgg& a, const ChildAgg& b) {
    return a.count == b.count && a.leaf_sum == b.leaf_sum &&
           a.min_node == b.min_node && a.max_node == b.max_node;
  }

  ValidationResult validate_tree(int node, int gene) const {
    if (node < 0) return ValidationResult{};
    const Node& n = nodes_[static_cast<std::size_t>(node)];
    const ValidationResult left = validate_tree(n.left, gene);
    const ValidationResult right = validate_tree(n.right, gene);
    ValidationResult out;
    out.empty = false;
    out.height = 1 + std::max(left.height, right.height);
    out.min_key = left.empty ? n.key : left.min_key;
    out.max_key = right.empty ? n.key : right.max_key;
    ChildAgg expected = combine_aggregates(left.agg, self_aggregate(node),
                                           gene);
    expected = combine_aggregates(expected, right.agg, gene);
    out.agg = expected;
    const bool ordered_left = left.empty || key_less(left.max_key, n.key);
    const bool ordered_right = right.empty || key_less(n.key, right.min_key);
    const bool balanced = std::abs(left.height - right.height) <= 1;
    out.valid = left.valid && right.valid && ordered_left && ordered_right &&
                balanced && n.height == out.height &&
                aggregate_equal(n.agg, expected) &&
                out.height <= height_limit_for_count(expected.count);
    return out;
  }

  void note_root(int root, int gene) {
    const int height = node_height(root);
    int& maximum = max_height_[static_cast<std::size_t>(gene)];
    maximum = std::max(maximum, height);
    const int count = root < 0
                          ? 0
                          : nodes_[static_cast<std::size_t>(root)].agg.count;
    if (height > height_limit_for_count(count)) {
      ++height_violations_[static_cast<std::size_t>(gene)];
      throw std::logic_error("MRE AVL height certificate failed");
    }
  }

  void enumerate_range_impl(int node,
                            const ChildKey& lo,
                            const ChildKey& hi,
                            int gene,
                            std::vector<int>& out) {
    if (node < 0) return;
    ++enumeration_visits_[static_cast<std::size_t>(gene)];
    const Node& n = nodes_[static_cast<std::size_t>(node)];
    if (key_less(n.key, lo)) {
      enumerate_range_impl(n.right, lo, hi, gene, out);
      return;
    }
    if (!key_less(n.key, hi)) {
      enumerate_range_impl(n.left, lo, hi, gene, out);
      return;
    }
    enumerate_range_impl(n.left, lo, hi, gene, out);
    out.push_back(n.child);
    enumerate_range_impl(n.right, lo, hi, gene, out);
  }

  void enumerate_all_impl(int node,
                          int gene,
                          std::vector<int>& out) {
    if (node < 0) return;
    ++enumeration_visits_[static_cast<std::size_t>(gene)];
    const Node& n = nodes_[static_cast<std::size_t>(node)];
    enumerate_all_impl(n.left, gene, out);
    out.push_back(n.child);
    enumerate_all_impl(n.right, gene, out);
  }

  void enumerate_prefix_impl(int node,
                             const ChildKey& bound,
                             int gene,
                             std::vector<int>& out) {
    if (node < 0) return;
    ++enumeration_visits_[static_cast<std::size_t>(gene)];
    const Node& n = nodes_[static_cast<std::size_t>(node)];
    if (!key_less(n.key, bound)) {
      enumerate_prefix_impl(n.left, bound, gene, out);
      return;
    }
    enumerate_all_impl(n.left, gene, out);
    out.push_back(n.child);
    enumerate_prefix_impl(n.right, bound, gene, out);
  }

  void enumerate_suffix_impl(int node,
                             const ChildKey& bound,
                             int gene,
                             std::vector<int>& out) {
    if (node < 0) return;
    ++enumeration_visits_[static_cast<std::size_t>(gene)];
    const Node& n = nodes_[static_cast<std::size_t>(node)];
    if (key_less(n.key, bound)) {
      enumerate_suffix_impl(n.right, bound, gene, out);
      return;
    }
    enumerate_suffix_impl(n.left, bound, gene, out);
    out.push_back(n.child);
    enumerate_all_impl(n.right, gene, out);
  }
};

// Mutable consensus topology with one recyclable slot per live parent-child
// edge. Detach deliberately implements the old vector swap-with-last rule, so
// iteration order is unchanged while historical per-parent capacities cannot
// accumulate across nested accepted clusters.
class FlatChildTopology {
 public:
  explicit FlatChildTopology(int max_nodes = 0) { reset(max_nodes); }

  void reset(int max_nodes) {
    if (max_nodes < 0) {
      throw std::invalid_argument("negative MRE topology size");
    }
    max_nodes_ = max_nodes;
    head_.assign(static_cast<std::size_t>(max_nodes_), -1);
    tail_.assign(static_cast<std::size_t>(max_nodes_), -1);
    size_.assign(static_cast<std::size_t>(max_nodes_), 0);
    child_slot_.assign(static_cast<std::size_t>(max_nodes_), -1);
    slots_.assign(static_cast<std::size_t>(max_nodes_), Slot{});
    free_.clear();
    free_.reserve(static_cast<std::size_t>(max_nodes_));
    next_slot_ = 0;
    live_slots_ = 0;
    peak_live_slots_ = 0;
  }

  void reset_node(int node) {
    require_node(node);
    if (head_[static_cast<std::size_t>(node)] >= 0 ||
        tail_[static_cast<std::size_t>(node)] >= 0 ||
        size_[static_cast<std::size_t>(node)] != 0 ||
        child_slot_[static_cast<std::size_t>(node)] >= 0) {
      throw std::logic_error("MRE topology reused a live node");
    }
  }

  void attach(int parent, int child) {
    require_node(parent);
    require_node(child);
    if (parent == child || child_slot_[static_cast<std::size_t>(child)] >= 0) {
      throw std::logic_error("MRE topology attached a non-free child");
    }
    const int slot = allocate_slot();
    const int old_tail = tail_[static_cast<std::size_t>(parent)];
    Slot& record = slots_[static_cast<std::size_t>(slot)];
    record.parent = parent;
    record.child = child;
    record.prev = old_tail;
    record.next = -1;
    if (old_tail >= 0) {
      slots_[static_cast<std::size_t>(old_tail)].next = slot;
    } else {
      head_[static_cast<std::size_t>(parent)] = slot;
    }
    tail_[static_cast<std::size_t>(parent)] = slot;
    child_slot_[static_cast<std::size_t>(child)] = slot;
    ++size_[static_cast<std::size_t>(parent)];
    ++live_slots_;
    peak_live_slots_ = std::max(peak_live_slots_, live_slots_);
  }

  void detach(int parent, int child) {
    require_node(parent);
    require_node(child);
    const int slot = child_slot_[static_cast<std::size_t>(child)];
    const int last = tail_[static_cast<std::size_t>(parent)];
    if (slot < 0 || last < 0 ||
        slots_[static_cast<std::size_t>(slot)].parent != parent) {
      throw std::logic_error("MRE topology detached a missing child");
    }

    if (slot != last) {
      const int replacement_child =
          slots_[static_cast<std::size_t>(last)].child;
      slots_[static_cast<std::size_t>(slot)].child = replacement_child;
      child_slot_[static_cast<std::size_t>(replacement_child)] = slot;
    }
    const int new_tail = slots_[static_cast<std::size_t>(last)].prev;
    if (new_tail >= 0) {
      slots_[static_cast<std::size_t>(new_tail)].next = -1;
    } else {
      head_[static_cast<std::size_t>(parent)] = -1;
    }
    tail_[static_cast<std::size_t>(parent)] = new_tail;
    child_slot_[static_cast<std::size_t>(child)] = -1;
    --size_[static_cast<std::size_t>(parent)];
    --live_slots_;
    release_slot(last);
  }

  int size(int parent) const {
    require_node(parent);
    return size_[static_cast<std::size_t>(parent)];
  }

  template <class F>
  void for_each_child(int parent, F&& body) const {
    require_node(parent);
    for (int slot = head_[static_cast<std::size_t>(parent)]; slot >= 0;
         slot = slots_[static_cast<std::size_t>(slot)].next) {
      body(slots_[static_cast<std::size_t>(slot)].child);
    }
  }

  std::vector<int> children(int parent) const {
    std::vector<int> out;
    out.reserve(static_cast<std::size_t>(size(parent)));
    for_each_child(parent, [&](int child) { out.push_back(child); });
    return out;
  }

  int allocated_slots() const { return next_slot_; }
  int live_slots() const { return live_slots_; }
  int peak_live_slots() const { return peak_live_slots_; }
  int slot_limit() const { return max_nodes_; }

  bool validate_for_tests() const {
    std::vector<uint8_t> seen_slot(static_cast<std::size_t>(next_slot_), 0);
    std::vector<uint8_t> seen_child(static_cast<std::size_t>(max_nodes_), 0);
    int listed_slots = 0;
    for (int parent = 0; parent < max_nodes_; ++parent) {
      int previous = -1;
      int listed_for_parent = 0;
      for (int slot = head_[static_cast<std::size_t>(parent)]; slot >= 0;
           slot = slots_[static_cast<std::size_t>(slot)].next) {
        if (slot >= next_slot_ ||
            seen_slot[static_cast<std::size_t>(slot)] != 0) {
          return false;
        }
        seen_slot[static_cast<std::size_t>(slot)] = 1;
        const Slot& record = slots_[static_cast<std::size_t>(slot)];
        if (record.parent != parent || record.prev != previous ||
            record.child < 0 || record.child >= max_nodes_ ||
            seen_child[static_cast<std::size_t>(record.child)] != 0 ||
            child_slot_[static_cast<std::size_t>(record.child)] != slot) {
          return false;
        }
        seen_child[static_cast<std::size_t>(record.child)] = 1;
        previous = slot;
        ++listed_for_parent;
        ++listed_slots;
        if (listed_slots > live_slots_) return false;
      }
      if (previous != tail_[static_cast<std::size_t>(parent)] ||
          listed_for_parent != size_[static_cast<std::size_t>(parent)] ||
          ((listed_for_parent == 0) !=
           (head_[static_cast<std::size_t>(parent)] < 0))) {
        return false;
      }
    }
    if (listed_slots != live_slots_) return false;
    std::vector<uint8_t> seen_free(static_cast<std::size_t>(next_slot_), 0);
    for (int slot : free_) {
      if (slot < 0 || slot >= next_slot_ ||
          seen_slot[static_cast<std::size_t>(slot)] != 0 ||
          seen_free[static_cast<std::size_t>(slot)] != 0) {
        return false;
      }
      seen_free[static_cast<std::size_t>(slot)] = 1;
    }
    if (static_cast<int>(free_.size()) + live_slots_ != next_slot_) {
      return false;
    }
    for (int child = 0; child < max_nodes_; ++child) {
      const int slot = child_slot_[static_cast<std::size_t>(child)];
      if ((slot >= 0) !=
          (seen_child[static_cast<std::size_t>(child)] != 0)) {
        return false;
      }
    }
    return peak_live_slots_ <= max_nodes_ && next_slot_ <= max_nodes_;
  }

 private:
  struct Slot {
    int parent = -1;
    int child = -1;
    int prev = -1;
    int next = -1;
  };

  int max_nodes_ = 0;
  std::vector<int> head_;
  std::vector<int> tail_;
  std::vector<int> size_;
  std::vector<int> child_slot_;
  std::vector<Slot> slots_;
  std::vector<int> free_;
  int next_slot_ = 0;
  int live_slots_ = 0;
  int peak_live_slots_ = 0;

  void require_node(int node) const {
    if (node < 0 || node >= max_nodes_) {
      throw std::out_of_range("MRE topology node outside flat domain");
    }
  }

  int allocate_slot() {
    if (!free_.empty()) {
      const int slot = free_.back();
      free_.pop_back();
      return slot;
    }
    if (next_slot_ >= max_nodes_) {
      throw std::logic_error("MRE topology exceeded its linear slot bound");
    }
    return next_slot_++;
  }

  void release_slot(int slot) {
    slots_[static_cast<std::size_t>(slot)] = Slot{};
    free_.push_back(slot);
  }
};

struct MreState {
  int n = 0;
  int k = 0;
  int root = 0;
  int node_count = 0;
  int max_nodes = 0;
  std::vector<int> parent;
  std::vector<int> taxon;
  std::vector<uint32_t> leaf_count;
  std::vector<uint32_t> degree;
  FlatChildTopology topology;
  std::vector<int> taxon_to_leaf;
  std::vector<int> lca_flat;
  std::vector<StaticLcaIndex> static_lca;
  LinkCutTree dyn;
  std::unique_ptr<ChildRangeForest> child_index;

  int& lca_in_gene(int node, int gene) {
    return lca_flat[static_cast<std::size_t>(node) *
                        static_cast<std::size_t>(k) +
                    static_cast<std::size_t>(gene)];
  }

  int lca_in_gene(int node, int gene) const {
    return lca_flat[static_cast<std::size_t>(node) *
                        static_cast<std::size_t>(k) +
                    static_cast<std::size_t>(gene)];
  }

  int new_node() {
    if (node_count >= max_nodes) {
      throw std::runtime_error("MRE consensus exceeded laminar node bound");
    }
    const int node = node_count++;
    if (node >= static_cast<int>(dyn_size_initialized)) {
      dyn.add_node();
      ++dyn_size_initialized;
    }
    parent[static_cast<std::size_t>(node)] = -1;
    taxon[static_cast<std::size_t>(node)] = -1;
    leaf_count[static_cast<std::size_t>(node)] = 0;
    degree[static_cast<std::size_t>(node)] = 0;
    topology.reset_node(node);
    return node;
  }

 private:
  int dyn_size_initialized = 0;
};

MreState initialize_state(const std::vector<GeneIndex>& genes,
                            int taxon_count) {
  validate_mre_dimensions(genes, taxon_count);
  MreState st;
  st.n = taxon_count;
  st.k = static_cast<int>(genes.size());
  st.max_nodes = std::max(2, 2 * taxon_count + 1);
  st.parent.assign(static_cast<std::size_t>(st.max_nodes), -1);
  st.taxon.assign(static_cast<std::size_t>(st.max_nodes), -1);
  st.leaf_count.assign(static_cast<std::size_t>(st.max_nodes), 0);
  st.degree.assign(static_cast<std::size_t>(st.max_nodes), 0);
  st.topology.reset(st.max_nodes);
  st.taxon_to_leaf.assign(static_cast<std::size_t>(taxon_count), -1);
  st.lca_flat.assign(
      checked_mre_size_product(
          static_cast<std::size_t>(st.max_nodes),
          static_cast<std::size_t>(std::max(1, st.k)),
          "MRE LCA state exceeds the addressable size domain"),
      -1);
  st.static_lca.resize(static_cast<std::size_t>(st.k));
  const auto static_lca_started = std::chrono::steady_clock::now();
  parallel_for_inner(0, st.k, [&](int g) {
    st.static_lca[static_cast<std::size_t>(g)] =
        build_static_lca_index(genes[static_cast<std::size_t>(g)], taxon_count);
  }, 1);
  debug_timer_add("mre.static_lca",
                  std::chrono::steady_clock::now() - static_lca_started);
  const auto child_index_started = std::chrono::steady_clock::now();
  st.dyn.reset(0);
  st.dyn.reserve(st.max_nodes);
  st.child_index =
      std::make_unique<ChildRangeForest>(st.max_nodes, st.k, st.static_lca);

  st.root = st.new_node();
  st.leaf_count[static_cast<std::size_t>(st.root)] =
      static_cast<uint32_t>(taxon_count);
  std::vector<std::vector<ChildRangeForest::InitialChild>> initial_children(
      static_cast<std::size_t>(st.k));
  for (int g = 0; g < st.k; ++g) {
    initial_children[static_cast<std::size_t>(g)].reserve(
        static_cast<std::size_t>(taxon_count));
  }
  for (int g = 0; g < st.k; ++g) {
    st.lca_in_gene(st.root, g) =
        genes[static_cast<std::size_t>(g)].reference_root;
  }

  for (int taxon = 0; taxon < taxon_count; ++taxon) {
    const int leaf = st.new_node();
    st.taxon[static_cast<std::size_t>(leaf)] = taxon;
    st.leaf_count[static_cast<std::size_t>(leaf)] = 1;
    st.taxon_to_leaf[static_cast<std::size_t>(taxon)] = leaf;
    st.parent[static_cast<std::size_t>(leaf)] = st.root;
    st.topology.attach(st.root, leaf);
    ++st.degree[static_cast<std::size_t>(st.root)];
    st.dyn.link(st.root, leaf);
    for (int g = 0; g < st.k; ++g) {
      const int lca_node =
          st.static_lca[static_cast<std::size_t>(g)]
              .taxon_to_node[static_cast<std::size_t>(taxon)];
      st.lca_in_gene(leaf, g) = lca_node;
      const int tin = lca_node >= 0
                          ? st.static_lca[static_cast<std::size_t>(g)]
                                .tin[static_cast<std::size_t>(lca_node)]
                          : -1;
      initial_children[static_cast<std::size_t>(g)].push_back(
          ChildRangeForest::InitialChild{leaf, tin, 1, lca_node});
    }
  }
  // Build the initial star's child AVL per gene in parallel: each gene writes
  // only its own pool segment and sorts its own initial_children list, so this
  // is race-free and bit-identical to the serial build.
  parallel_for_inner(0, st.k, [&](int g) {
    st.child_index->build_children(
        st.root, g, initial_children[static_cast<std::size_t>(g)]);
  }, 1);
  debug_timer_add("mre.initial_child_index",
                  std::chrono::steady_clock::now() - child_index_started);
  return st;
}

bool subset_node_in_candidate(const MreState& st,
                              int node,
                              int gene,
                              int rep_node) {
  const int lca = st.lca_in_gene(node, gene);
  return st.static_lca[static_cast<std::size_t>(gene)].in_subtree(lca,
                                                                  rep_node);
}

int find_u0(MreState& st, int leaf, int gene, int rep_node) {
  const int depth = st.dyn.path_size(st.root, leaf) - 1;
  int lo = 0;
  int hi = depth;
  while (lo < hi) {
    const int mid = (lo + hi + 1) >> 1;
    const int u = st.dyn.kth_on_path(st.root, leaf, mid);
    if (subset_node_in_candidate(st, u, gene, rep_node)) {
      hi = mid - 1;
    } else {
      lo = mid;
    }
  }
  return st.dyn.kth_on_path(st.root, leaf, lo);
}

void detach_topology(MreState& st, int parent, int child) {
  st.topology.detach(parent, child);
  st.parent[static_cast<std::size_t>(child)] = -1;
  --st.degree[static_cast<std::size_t>(parent)];
  st.dyn.cut(parent, child);
}

void attach_topology(MreState& st, int parent, int child) {
  st.parent[static_cast<std::size_t>(child)] = parent;
  st.topology.attach(parent, child);
  ++st.degree[static_cast<std::size_t>(parent)];
  st.dyn.link(parent, child);
}

int child_key_tin(const MreState& st, int child, int gene) {
  const int lca = st.lca_in_gene(child, gene);
  if (lca < 0) return -1;
  return st.static_lca[static_cast<std::size_t>(gene)]
      .tin[static_cast<std::size_t>(lca)];
}

void erase_child_index(MreState& st, int parent, int gene, int child) {
  st.child_index->erase_child(parent, gene, child,
                              child_key_tin(st, child, gene));
}

void insert_child_index(MreState& st, int parent, int gene, int child) {
  st.child_index->insert_child(
      parent, gene, child, child_key_tin(st, child, gene),
      st.leaf_count[static_cast<std::size_t>(child)],
      st.lca_in_gene(child, gene));
}

// Runs body(g) for every gene. When the estimated per-gene index work is large
// enough to amortise fork/join overhead, distributes the genes across the
// shared executor; otherwise runs serially. Each gene owns a disjoint AVL pool
// (see ChildRangeForest's per-gene segments) and only reads shared state, so
// the bodies never race. AVL updates and bulk-build pivots are deterministic,
// so the parallel result is bit-identical to the serial one.
//
// The threshold is the estimated number of per-gene index operations
// (moved_children * genes) below which fork/join overhead outweighs the gain.
// Tuned so small-k runs stay serial (no regression) while large-k runs (the
// regime where consensus dominates on high core counts) parallelise the
// expensive large-move accepts. Override with WQFM_CONSENSUS_PAR_THRESHOLD.
long long consensus_gene_parallel_threshold() {
  static const long long value = [] {
    if (const char* env = std::getenv("WQFM_CONSENSUS_PAR_THRESHOLD")) {
      char* end = nullptr;
      const long long parsed = std::strtoll(env, &end, 10);
      if (end != env && parsed > 0) return parsed;
    }
    return 1024LL;
  }();
  return value;
}

template <class F>
void for_each_gene(int k, long long work_estimate, F&& body) {
  if (k > 1 && work_estimate >= consensus_gene_parallel_threshold()) {
    parallel_for_inner(0, k, std::forward<F>(body), 1);
  } else {
    for (int g = 0; g < k; ++g) body(g);
  }
}

// Rebuilds a parent's child AVL for one gene from the current topology.
// Uses a thread-local scratch buffer so concurrent gene tasks don't contend.
void build_parent_gene_index(MreState& st, int parent, int g) {
  std::vector<ChildRangeForest::InitialChild> initial;
  initial.reserve(static_cast<std::size_t>(st.topology.size(parent)));
  st.topology.for_each_child(parent, [&](int child) {
    initial.push_back(ChildRangeForest::InitialChild{
        child,
        child_key_tin(st, child, g),
        st.leaf_count[static_cast<std::size_t>(child)],
        st.lca_in_gene(child, g)});
  });
  st.child_index->build_children(parent, g, initial);
}

void accept_move_contained(MreState& st,
                           const MreEntry& entry,
                           int u0,
                           const std::vector<int>& contained) {
  const int z = st.new_node();
  st.leaf_count[static_cast<std::size_t>(z)] = entry.key.size;
  const bool batch_new_index = contained.size() >= 8;

  // Pass 1 (serial): topology only -- the link-cut tree and children/parent/
  // degree arrays are shared, non-per-gene state.
  for (int child : contained) {
    detach_topology(st, u0, child);
    attach_topology(st, z, child);
  }
  attach_topology(st, u0, z);

  // Pass 2 (parallel over genes): each gene's AVL pool is disjoint; st is read
  // only here (children/leaf_count/lca of moved children are intrinsic and
  // unchanged by the topology pass), so this is race-free.
  const long long work = static_cast<long long>(contained.size()) * st.k;
  for_each_gene(st.k, work, [&](int g) {
    for (int child : contained) erase_child_index(st, u0, g, child);
    if (batch_new_index) {
      build_parent_gene_index(st, z, g);
    } else {
      for (int child : contained) insert_child_index(st, z, g, child);
    }
    st.lca_in_gene(z, g) = st.child_index->aggregate_lca(z, g);
    insert_child_index(st, u0, g, z);
  });
}

void accept_move_complement(MreState& st,
                            const MreEntry& entry,
                            int u0,
                            const std::vector<int>& complement) {
  const int old_parent = st.parent[static_cast<std::size_t>(u0)];
  const uint32_t old_leaf_count = st.leaf_count[static_cast<std::size_t>(u0)];

  const int z = st.new_node();
  st.leaf_count[static_cast<std::size_t>(z)] = old_leaf_count;

  // Pass 1 (serial): topology only.
  if (old_parent >= 0) {
    detach_topology(st, old_parent, u0);
    attach_topology(st, old_parent, z);
  } else {
    st.root = z;
  }
  const bool batch_new_index = complement.size() >= 8;
  for (int child : complement) {
    detach_topology(st, u0, child);
    attach_topology(st, z, child);
  }
  st.leaf_count[static_cast<std::size_t>(u0)] = entry.key.size;
  attach_topology(st, z, u0);

  // Pass 2 (parallel over genes): z inherits u0's old per-gene LCA (same leaf
  // set); u0 is recomputed from its remaining (kept) children. Per-gene pools
  // are disjoint and the lca writes for z/u0 are at gene-disjoint slots.
  const long long work = static_cast<long long>(complement.size() + 1) * st.k;
  for_each_gene(st.k, work, [&](int g) {
    const int z_lca = st.lca_in_gene(u0, g);  // u0's old LCA == z's LCA
    if (old_parent >= 0) {
      erase_child_index(st, old_parent, g, u0);
      st.lca_in_gene(z, g) = z_lca;
      insert_child_index(st, old_parent, g, z);
    } else {
      st.lca_in_gene(z, g) = z_lca;
    }
    for (int child : complement) erase_child_index(st, u0, g, child);
    st.lca_in_gene(u0, g) = st.child_index->aggregate_lca(u0, g);
    if (batch_new_index) {
      build_parent_gene_index(st, z, g);
    } else {
      for (int child : complement) insert_child_index(st, z, g, child);
      insert_child_index(st, z, g, u0);
    }
  });
}

void fill_traversal(ConsensusTree& consensus) {
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
                   consensus.nodes[static_cast<std::size_t>(node)].children.size())) {
      stack.push_back({consensus.nodes[static_cast<std::size_t>(node)]
                           .children[static_cast<std::size_t>(next++)],
                       0});
      continue;
    }
    consensus.tout[static_cast<std::size_t>(node)] =
        static_cast<int>(consensus.preorder.size());
    stack.pop_back();
  }
}

ConsensusTree materialize_consensus(const MreState& st) {
  ConsensusTree out;
  out.root = st.root;
  out.nodes.resize(static_cast<std::size_t>(st.node_count));
  out.taxon_to_leaf.assign(static_cast<std::size_t>(st.n), -1);
  for (int node = 0; node < st.node_count; ++node) {
    out.nodes[static_cast<std::size_t>(node)].parent =
        st.parent[static_cast<std::size_t>(node)];
    out.nodes[static_cast<std::size_t>(node)].taxon =
        st.taxon[static_cast<std::size_t>(node)];
    out.nodes[static_cast<std::size_t>(node)].children =
        st.topology.children(node);
    if (st.taxon[static_cast<std::size_t>(node)] >= 0) {
      out.taxon_to_leaf[static_cast<std::size_t>(
          st.taxon[static_cast<std::size_t>(node)])] = node;
    }
  }
  fill_traversal(out);
  return out;
}

}  // namespace

MreComplexityProbeResult mre_complexity_probe_for_tests(int item_count) {
  if (item_count < 2) {
    throw std::invalid_argument(
        "MRE complexity probe requires at least two items");
  }

  std::vector<SplitOccurrence256> occurrences;
  occurrences.reserve(static_cast<std::size_t>(item_count) * 2u);
  for (int i = 0; i < item_count; ++i) {
    const uint64_t index = static_cast<uint64_t>(i);
    SplitOccurrence256 occurrence;
    occurrence.key.size = static_cast<uint32_t>(
        stable_hash_u64(index % 23, index / 23, 0x082efa98ec4e6c89ULL));
    occurrence.key.h[0] =
        stable_hash_u64(index % 29, index / 29, 0x452821e638d01377ULL);
    occurrence.key.h[1] =
        stable_hash_u64(index % 31, index / 31, 0xbe5466cf34e90c6cULL);
    occurrence.key.h[2] =
        stable_hash_u64(index % 41, index / 41, 0xc0ac29b7c97c50ddULL);
    occurrence.key.h[3] = index;
    occurrence.gene = static_cast<uint32_t>((37LL * i + 11) % item_count);
    occurrence.edge = static_cast<uint32_t>(
        stable_hash_u64(index % 97, index / 97, 0xa4093822299f31d0ULL));
    occurrence.node = static_cast<uint32_t>(i);
    occurrence.one_taxon = static_cast<uint32_t>(item_count - i);
    occurrences.push_back(occurrence);
    if (i % 3 == 0) {
      SplitOccurrence256 duplicate = occurrence;
      duplicate.edge += 1;
      duplicate.node ^= 0x5a5a5a5aU;
      occurrences.push_back(duplicate);
    }
    if (i % 5 == 0) {
      SplitOccurrence256 second_gene = occurrence;
      second_gene.gene = static_cast<uint32_t>(
          (static_cast<uint64_t>(occurrence.gene) + 1u) %
          static_cast<uint64_t>(item_count));
      second_gene.edge ^= 0x80000000U;
      occurrences.push_back(second_gene);
    }
  }
  std::reverse(occurrences.begin(), occurrences.end());
  std::rotate(
      occurrences.begin(),
      occurrences.begin() +
          static_cast<std::ptrdiff_t>(occurrences.size() / 3u),
      occurrences.end());

  const auto occurrence_less_oracle = [](const SplitOccurrence256& a,
                                         const SplitOccurrence256& b) {
    if (a.key.size != b.key.size) return a.key.size < b.key.size;
    if (a.key.h != b.key.h) return a.key.h < b.key.h;
    if (a.gene != b.gene) return a.gene < b.gene;
    return a.edge < b.edge;
  };
  std::vector<SplitOccurrence256> ordered_occurrence_oracle = occurrences;
  std::sort(ordered_occurrence_oracle.begin(), ordered_occurrence_oracle.end(),
            occurrence_less_oracle);

  std::vector<MreEntry> aggregate_oracle;
  aggregate_oracle.reserve(static_cast<std::size_t>(item_count));
  std::size_t oracle_begin = 0;
  while (oracle_begin < ordered_occurrence_oracle.size()) {
    const SplitOccurrence256& first = ordered_occurrence_oracle[oracle_begin];
    MreEntry aggregate;
    aggregate.key = first.key;
    aggregate.split_id = static_cast<uint32_t>(aggregate_oracle.size());
    aggregate.first_tree = first.gene;
    aggregate.first_edge = first.edge;
    aggregate.rep_tree = first.gene;
    aggregate.rep_node = first.node;
    aggregate.one_taxon = first.one_taxon;
    uint32_t previous_gene = UINT32_MAX;
    std::size_t oracle_end = oracle_begin;
    while (oracle_end < ordered_occurrence_oracle.size() &&
           ordered_occurrence_oracle[oracle_end].key == first.key) {
      if (ordered_occurrence_oracle[oracle_end].gene != previous_gene) {
        previous_gene = ordered_occurrence_oracle[oracle_end].gene;
        ++aggregate.count;
      }
      ++oracle_end;
    }
    aggregate_oracle.push_back(aggregate);
    oracle_begin = oracle_end;
  }

  SplitRadixStats radix_stats;
  std::vector<MreEntry> entries =
      aggregate_split_occurrences256(occurrences, radix_stats);
  const auto occurrence_equal = [](const SplitOccurrence256& a,
                                   const SplitOccurrence256& b) {
    return a.key == b.key && a.gene == b.gene && a.edge == b.edge &&
           a.node == b.node && a.one_taxon == b.one_taxon;
  };
  const auto entry_equal = [](const MreEntry& a, const MreEntry& b) {
    return a.key == b.key && a.split_id == b.split_id &&
           a.count == b.count && a.first_tree == b.first_tree &&
           a.first_edge == b.first_edge && a.rep_tree == b.rep_tree &&
           a.rep_node == b.rep_node && a.one_taxon == b.one_taxon;
  };
  const bool occurrence_order_matches = std::equal(
      occurrences.begin(), occurrences.end(),
      ordered_occurrence_oracle.begin(), occurrence_equal);
  const bool aggregation_matches =
      entries.size() == aggregate_oracle.size() &&
      std::equal(entries.begin(), entries.end(), aggregate_oracle.begin(),
                 entry_equal);
  bool dense_ids_match = true;
  for (std::size_t i = 0; i < entries.size(); ++i) {
    dense_ids_match = dense_ids_match && entries[i].split_id == i;
  }

  std::vector<MreEntry> priority_oracle = aggregate_oracle;
  std::sort(priority_oracle.begin(), priority_oracle.end(), entry_less);
  radix_sort_split_priority256(entries, radix_stats);
  const bool priority_matches =
      entries.size() == priority_oracle.size() &&
      std::equal(entries.begin(), entries.end(), priority_oracle.begin(),
                 entry_equal);
  const bool radix_matches = occurrence_order_matches && aggregation_matches &&
                             dense_ids_match && priority_matches;

  StaticLcaIndex lca;
  lca.tin.resize(static_cast<std::size_t>(item_count + 1));
  for (int child = 0; child <= item_count; ++child) {
    lca.tin[static_cast<std::size_t>(child)] = child;
  }
  std::vector<StaticLcaIndex> lca_by_gene;
  lca_by_gene.push_back(std::move(lca));
  ChildRangeForest forest(item_count + 1, 1, lca_by_gene);
  std::vector<uint8_t> live(static_cast<std::size_t>(item_count + 1), 0);
  bool avl_valid = true;
  bool aggregates_valid = true;
  long long empty_complement_visits = 0;

  const auto leaf_mass = [](int child) {
    return static_cast<uint32_t>(1 + (child % 11));
  };
  const auto validate = [&]() {
    avl_valid = avl_valid && forest.validate_root_for_tests(0, 0);
  };
  const auto check_range = [&](int lo, int hi) {
    ChildAgg expected;
    std::vector<int> expected_inside;
    std::vector<int> expected_outside;
    for (int child = 1; child <= item_count; ++child) {
      if (live[static_cast<std::size_t>(child)] == 0) continue;
      if (child >= lo && child < hi) {
        ++expected.count;
        expected.leaf_sum += leaf_mass(child);
        if (expected.min_node < 0) expected.min_node = child;
        expected.max_node = child;
        expected_inside.push_back(child);
      } else {
        expected_outside.push_back(child);
      }
    }
    const ChildAgg observed = forest.range_query(0, 0, lo, hi);
    const std::vector<int> inside = forest.enumerate_range(0, 0, lo, hi);
    const long long before_complement = forest.enumeration_visit_count();
    const std::vector<int> outside =
        forest.enumerate_complement(0, 0, lo, hi);
    const long long complement_visits =
        forest.enumeration_visit_count() - before_complement;
    if (expected.count > 0 && expected_outside.empty()) {
      empty_complement_visits =
          std::max(empty_complement_visits, complement_visits);
    }
    aggregates_valid =
        aggregates_valid && observed.count == expected.count &&
        observed.leaf_sum == expected.leaf_sum &&
        observed.min_node == expected.min_node &&
        observed.max_node == expected.max_node && inside == expected_inside &&
        outside == expected_outside;
  };

  // Monotone insertion is adversarial for an unbalanced BST and for a treap
  // whose fixed priorities happen to be monotone in key order.
  for (int child = 1; child <= item_count; ++child) {
    forest.insert_child(0, 0, child, child, leaf_mass(child), child);
    live[static_cast<std::size_t>(child)] = 1;
  }
  validate();
  check_range(1, item_count + 1);
  check_range(item_count / 4, 3 * item_count / 4);
  check_range(item_count + 1, item_count + 2);

  for (int child = 1; child <= item_count; child += 2) {
    forest.erase_child(0, 0, child, child);
    live[static_cast<std::size_t>(child)] = 0;
  }
  validate();
  check_range(item_count / 5, 4 * item_count / 5);
  for (int child = item_count - (item_count & 1); child >= 2; child -= 2) {
    forest.erase_child(0, 0, child, child);
    live[static_cast<std::size_t>(child)] = 0;
  }
  validate();
  check_range(1, item_count + 1);

  for (int child = item_count; child >= 1; --child) {
    forest.insert_child(0, 0, child, child, leaf_mass(child), child);
    live[static_cast<std::size_t>(child)] = 1;
  }
  validate();
  for (int child = 3; child <= item_count; child += 3) {
    forest.erase_child(0, 0, child, child);
    live[static_cast<std::size_t>(child)] = 0;
  }
  validate();
  check_range(item_count / 7, 6 * item_count / 7);
  for (int child = 3; child <= item_count; child += 3) {
    forest.insert_child(0, 0, child, child, leaf_mass(child), child);
    live[static_cast<std::size_t>(child)] = 1;
  }
  validate();
  check_range(1, item_count + 1);

  FlatChildTopology topology(item_count + 2);
  std::vector<int> root_children;
  std::vector<int> nested_children;
  root_children.reserve(static_cast<std::size_t>(item_count + 1));
  nested_children.reserve(static_cast<std::size_t>(item_count));
  for (int child = 1; child <= item_count + 1; ++child) {
    topology.attach(0, child);
    root_children.push_back(child);
  }
  const auto oracle_detach = [](std::vector<int>& children, int child) {
    const auto found = std::find(children.begin(), children.end(), child);
    if (found == children.end()) return false;
    *found = children.back();
    children.pop_back();
    return true;
  };
  for (int child = 2; child <= item_count; child += 3) {
    topology.detach(0, child);
    if (!oracle_detach(root_children, child)) {
      throw std::logic_error("MRE topology probe oracle lost a child");
    }
    topology.attach(1, child);
    nested_children.push_back(child);
  }
  const std::vector<int> nested_snapshot = nested_children;
  for (std::size_t i = 1; i < nested_snapshot.size(); i += 5) {
    const int child = nested_snapshot[i];
    topology.detach(1, child);
    if (!oracle_detach(nested_children, child)) {
      throw std::logic_error("MRE topology probe oracle lost a move");
    }
    topology.attach(0, child);
    root_children.push_back(child);
  }
  const bool topology_matches = topology.validate_for_tests() &&
                                topology.children(0) == root_children &&
                                topology.children(1) == nested_children;

  MreComplexityProbeResult result;
  result.radix_matches_oracle = radix_matches;
  result.avl_invariants_hold =
      avl_valid && forest.height_violation_count() == 0;
  result.range_aggregates_match = aggregates_valid;
  result.radix_passes = static_cast<int>(radix_stats.total_passes());
  result.radix_items = radix_stats.total_items();
  result.occurrence_radix_passes =
      static_cast<int>(radix_stats.occurrence_passes);
  result.priority_radix_passes =
      static_cast<int>(radix_stats.priority_passes);
  result.occurrence_radix_items = radix_stats.occurrence_items;
  result.priority_radix_items = radix_stats.priority_items;
  result.emitted_occurrences = radix_stats.emitted_occurrences;
  result.deduplicated_gene_splits = radix_stats.deduplicated_gene_splits;
  result.distinct_splits = radix_stats.distinct_splits;
  result.dense_split_ids = dense_ids_match;
  result.max_avl_height = forest.maximum_height();
  result.avl_height_limit =
      ChildRangeForest::height_limit_for_count(item_count);
  result.empty_complement_visits = empty_complement_visits;
  result.empty_complement_visit_limit =
      2LL * result.avl_height_limit + 2;
  result.flat_topology_matches = topology_matches;
  result.topology_peak_slots = topology.peak_live_slots();
  result.topology_slot_limit = topology.slot_limit();
  return result;
}

ConsensusTree build_consensus(const std::vector<GeneIndex>& genes,
                                     int taxon_count,
                                     const Options& options) {
  DebugScope scope("solver.build_consensus");
  validate_mre_dimensions(genes, taxon_count);
  progress_consensus_start(options, static_cast<int>(genes.size()),
                           taxon_count);
  const auto count_started = std::chrono::steady_clock::now();
  const std::vector<MreEntry> entries =
      count_mre_entries(genes, taxon_count, options);
  BenchmarkStats& stats = benchmark_stats();
  benchmark_timer_add(stats.consensus_count_ns,
                      std::chrono::duration_cast<std::chrono::nanoseconds>(
                          std::chrono::steady_clock::now() - count_started));
  stats.consensus_distinct_splits.store(static_cast<long long>(entries.size()),
                                        std::memory_order_relaxed);

  const auto assembly_started = std::chrono::steady_clock::now();
  MreState st = initialize_state(genes, taxon_count);
  int accepted = 0;
  long long candidates_checked = 0;
  long long candidates_rejected = 0;
  long long case_a_accepts = 0;
  long long case_b_accepts = 0;
  long long moved_children_total = 0;
  long long max_child_degree = 0;
  uint64_t consensus_hash = 1469598103934665603ULL;
  const int max_accepted_splits = std::max(0, taxon_count - 3);
  const bool debug_mre_timing = debug_timing_enabled();
  progress_consensus_assembly(0, static_cast<long long>(entries.size()), 0,
                              max_accepted_splits);
  for (const MreEntry& entry : entries) {
    if (accepted >= max_accepted_splits) break;
    // Priority order is support-descending, so the first failed strict cutoff
    // also excludes every later entry. The total-gene denominator is the MSC
    // split-probability denominator only for complete single-copy gene trees.
    if (options.hat_one_third_singleton &&
        static_cast<uint64_t>(entry.count) * 3u <=
            static_cast<uint64_t>(genes.size())) {
      break;
    }
    ++candidates_checked;
    const auto note_assembly_progress = [&]() {
      if ((candidates_checked & 1023ll) == 0 ||
          accepted >= max_accepted_splits) {
        progress_consensus_assembly(
            candidates_checked, static_cast<long long>(entries.size()),
            accepted, max_accepted_splits);
      }
    };
    const int gene = static_cast<int>(entry.rep_tree);
    const int rep_node = static_cast<int>(entry.rep_node);
    const int taxon = static_cast<int>(entry.one_taxon);
    if (taxon < 0 || taxon >= taxon_count) {
      ++candidates_rejected;
      note_assembly_progress();
      continue;
    }
    const int leaf = st.taxon_to_leaf[static_cast<std::size_t>(taxon)];
    if (leaf < 0) {
      ++candidates_rejected;
      note_assembly_progress();
      continue;
    }
    std::chrono::steady_clock::time_point find_started;
    if (debug_mre_timing) find_started = std::chrono::steady_clock::now();
    const int u0 = find_u0(st, leaf, gene, rep_node);
    if (debug_mre_timing) {
      debug_timer_add("mre.find_u0",
                      std::chrono::steady_clock::now() - find_started);
    }
    const StaticLcaIndex& lca = st.static_lca[static_cast<std::size_t>(gene)];
    const int lo = lca.tin[static_cast<std::size_t>(rep_node)];
    const int hi = lca.tout[static_cast<std::size_t>(rep_node)];
    std::chrono::steady_clock::time_point range_started;
    if (debug_mre_timing) range_started = std::chrono::steady_clock::now();
    const ChildAgg contained = st.child_index->range_query(u0, gene, lo, hi);
    if (debug_mre_timing) {
      debug_timer_add("mre.range_query",
                      std::chrono::steady_clock::now() - range_started);
    }
    if (contained.count == 0 ||
        contained.leaf_sum != static_cast<uint64_t>(entry.key.size) ||
        contained.leaf_sum == st.leaf_count[static_cast<std::size_t>(u0)]) {
      ++candidates_rejected;
      note_assembly_progress();
      continue;
    }
    max_child_degree =
        std::max<long long>(max_child_degree,
                            st.degree[static_cast<std::size_t>(u0)]);
    if (static_cast<uint64_t>(contained.count) * 2u <=
        st.degree[static_cast<std::size_t>(u0)]) {
      std::chrono::steady_clock::time_point enum_started;
      if (debug_mre_timing) enum_started = std::chrono::steady_clock::now();
      std::vector<int> moved =
          st.child_index->enumerate_range(u0, gene, lo, hi);
      if (debug_mre_timing) {
        debug_timer_add("mre.enumerate",
                        std::chrono::steady_clock::now() - enum_started);
      }
      moved_children_total += static_cast<long long>(moved.size());
      std::chrono::steady_clock::time_point accept_started;
      if (debug_mre_timing) accept_started = std::chrono::steady_clock::now();
      accept_move_contained(st, entry, u0, moved);
      if (debug_mre_timing) {
        debug_timer_add("mre.accept",
                        std::chrono::steady_clock::now() - accept_started);
      }
      ++case_a_accepts;
    } else {
      std::chrono::steady_clock::time_point enum_started;
      if (debug_mre_timing) enum_started = std::chrono::steady_clock::now();
      std::vector<int> moved =
          st.child_index->enumerate_complement(u0, gene, lo, hi);
      if (debug_mre_timing) {
        debug_timer_add("mre.enumerate",
                        std::chrono::steady_clock::now() - enum_started);
      }
      moved_children_total += static_cast<long long>(moved.size());
      std::chrono::steady_clock::time_point accept_started;
      if (debug_mre_timing) accept_started = std::chrono::steady_clock::now();
      accept_move_complement(st, entry, u0, moved);
      if (debug_mre_timing) {
        debug_timer_add("mre.accept",
                        std::chrono::steady_clock::now() - accept_started);
      }
      ++case_b_accepts;
    }
    ++accepted;
    consensus_hash = stable_hash_u64(consensus_hash, entry.count,
                                     entry.key.h[0] ^ entry.key.h[1] ^
                                         entry.key.h[2] ^ entry.key.h[3]);
    note_assembly_progress();
  }
  progress_consensus_finish(accepted, max_accepted_splits);
  benchmark_timer_add(stats.consensus_assembly_ns,
                      std::chrono::duration_cast<std::chrono::nanoseconds>(
                          std::chrono::steady_clock::now() - assembly_started));
  stats.consensus_accepted_splits.store(accepted, std::memory_order_relaxed);
  stats.consensus_hash_low.store(
      static_cast<long long>(consensus_hash & 0x7fffffffffffffffULL),
      std::memory_order_relaxed);
  stats.mre_candidates_checked.fetch_add(candidates_checked,
                                            std::memory_order_relaxed);
  stats.mre_candidates_rejected.fetch_add(candidates_rejected,
                                             std::memory_order_relaxed);
  stats.mre_case_a_accepts.fetch_add(case_a_accepts,
                                        std::memory_order_relaxed);
  stats.mre_case_b_accepts.fetch_add(case_b_accepts,
                                        std::memory_order_relaxed);
  stats.mre_moved_children.fetch_add(moved_children_total,
                                        std::memory_order_relaxed);
  benchmark_counter_max(stats.mre_max_child_degree, max_child_degree);
  stats.mre_child_index_operations.fetch_add(
      st.child_index->operation_count(), std::memory_order_relaxed);
  stats.mre_child_index_rotations.fetch_add(
      st.child_index->rotation_count(), std::memory_order_relaxed);
  benchmark_counter_max(stats.mre_child_index_max_height,
                        st.child_index->maximum_height());
  const long long height_violations =
      st.child_index->height_violation_count();
  stats.mre_child_index_height_violations.fetch_add(
      height_violations, std::memory_order_relaxed);
  stats.mre_child_index_enumeration_visits.fetch_add(
      st.child_index->enumeration_visit_count(), std::memory_order_relaxed);
  if (st.topology.live_slots() != st.node_count - 1 ||
      st.topology.allocated_slots() > st.topology.slot_limit() ||
      st.topology.peak_live_slots() > st.topology.slot_limit()) {
    stats.complexity_contract_satisfied.store(0, std::memory_order_relaxed);
    throw std::logic_error("MRE flat topology certificate failed");
  }
  stats.mre_topology_slots_allocated.fetch_add(
      st.topology.allocated_slots(), std::memory_order_relaxed);
  benchmark_counter_max(stats.mre_topology_peak_live_slots,
                        st.topology.peak_live_slots());
  if (height_violations != 0) {
    stats.complexity_contract_satisfied.store(0, std::memory_order_relaxed);
  }
  debug_counter_add("consensus.mre_distinct_splits",
                    static_cast<long long>(entries.size()));
  debug_counter_add("consensus.mre_accepted_splits", accepted);
  if (!options.debug_trace.empty()) {
    std::ostringstream trace;
    trace << "{\"event\":\"mre_consensus\",\"taxa\":" << taxon_count
          << ",\"distinct_splits\":" << entries.size()
          << ",\"accepted_splits\":" << accepted
          << ",\"strict_one_third\":"
          << (options.hat_one_third_singleton ? "true" : "false")
          << ",\"hash_low\":" << (consensus_hash & 0x7fffffffffffffffULL)
          << "}";
    trace_json_line(options, 0, trace.str());
  }
  return materialize_consensus(st);
}

}  // namespace wqfm::internal
