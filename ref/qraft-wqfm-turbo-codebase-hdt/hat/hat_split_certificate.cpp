#include "hat/hat_split_certificate.hpp"
#include "hat/hdt/hdt_checked.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <numbers>
#include <stdexcept>

namespace wqfm::internal {

int HatGeneDistanceIndex::distance(const GeneIndex& gene, int a, int b) const {
  if (a < 0 || b < 0 || static_cast<std::size_t>(a) >= depth.size() ||
      static_cast<std::size_t>(b) >= depth.size()) {
    throw std::logic_error("split certificate has an invalid distance endpoint");
  }
  const long long endpoints = static_cast<long long>(depth[a]) + depth[b];
  // Each jump crosses a light edge, so there are O(log |G|) jumps. Unlike
  // binary lifting, the two arrays use O(|G|), not O(|G| log |G|), words.
  while (head[a] != head[b]) {
    if (depth[head[a]] < depth[head[b]]) std::swap(a, b);
    a = gene.root_traversal.parent[head[a]];
  }
  const long long result = endpoints - 2LL * std::min(depth[a], depth[b]);
  if (result < 0 || result > std::numeric_limits<int>::max()) {
    throw std::logic_error("split certificate distance exceeds its integer domain");
  }
  return static_cast<int>(result);
}

std::shared_ptr<const HatSplitCertificateIndex> make_hat_split_certificate_index(
    const std::vector<GeneIndex>& genes, int taxon_count) {
  const auto started = std::chrono::steady_clock::now();
  auto index = std::make_shared<HatSplitCertificateIndex>();
  index->taxon_count = taxon_count;
  if (taxon_count < 4 || genes.empty()) return index;
  const int gene_count = checked_hdt_int_size(genes.size(), "too many genes");
  std::vector<unsigned char> valid(genes.size(), 1);
  parallel_for_inner(0, gene_count, [&](int g) {
    const GeneIndex& gene = genes[g];
    std::vector<unsigned char> seen(static_cast<std::size_t>(taxon_count), 0);
    int count = 0;
    for (const GeneIndex::Node& node : gene.nodes) {
      if (node.taxon < 0) continue;
      // Do not silently select one physical occurrence using taxon_to_node.
      // The distance theorem used here requires one copy of every taxon.
      if (node.taxon >= taxon_count || seen[node.taxon] ||
          node.neighbors.size() != 1) {
        valid[g] = 0;
        break;
      }
      seen[node.taxon] = 1;
      ++count;
    }
    if (count != taxon_count) valid[g] = 0;
  }, 1);
  index->complete_single_copy =
      std::all_of(valid.begin(), valid.end(), [](auto x) { return x != 0; });
  if (index->complete_single_copy) {
    index->genes.resize(genes.size());
    parallel_for_inner(0, gene_count, [&](int g) {
      const GeneIndex& gene = genes[g];
      HatGeneDistanceIndex& out = index->genes[g];
      const int size = checked_hdt_int_size(gene.nodes.size(), "gene too large");
      out.depth.assign(gene.nodes.size(), 0);
      out.head.assign(gene.nodes.size(), gene.root);
      std::vector<int> subtree(gene.nodes.size(), 1);
      std::vector<int> heavy(gene.nodes.size(), -1);
      for (int u : gene.root_traversal.postorder) {
        const int parent = gene.root_traversal.parent[u];
        if (u == gene.root) continue;
        if (parent < 0 || parent >= size) {
          throw std::logic_error("split certificate requires a complete traversal");
        }
        subtree[parent] += subtree[u];
        if (heavy[parent] < 0 || subtree[u] > subtree[heavy[parent]]) {
          heavy[parent] = u;
        }
      }
      for (int u : gene.root_traversal.preorder) {
        if (u == gene.root) continue;
        const int parent = gene.root_traversal.parent[u];
        out.depth[u] = out.depth[parent] + 1;
        out.head[u] = heavy[parent] == u ? out.head[parent] : u;
      }
    }, 1);
  }
  benchmark_timer_add(benchmark_stats().hat_certificate_index_ns,
                      std::chrono::steady_clock::now() - started);
  return index;
}

double hat_split_certificate_threshold(int gene_count, int taxon_count,
                                       double coefficient, double crossover) {
  if (gene_count <= 0 || taxon_count < 2 || !std::isfinite(coefficient) ||
      coefficient <= 0.0 || !std::isfinite(crossover) || crossover < 0.0) {
    throw std::invalid_argument(
        "split certificate needs genes, taxa, a finite positive coefficient, "
        "and a finite nonnegative crossover");
  }
  if (gene_count == crossover) return 0.0;
  const long double k = gene_count;
  const long double k0 = crossover;
  const long double e = std::numbers::e_v<long double>;
  const long double log_k = std::log(k + e);
  const long double log_k0 = std::log(k0 + e);
  // Rationalize the square-root difference. Near the crossover, log1p
  // avoids cancellation and preserves the sign on either side of K0.
  const long double log_ratio = std::abs(k - k0) <= (k0 + e) / 2.0L
      ? std::log1p((k - k0) / (k0 + e)) : log_k - log_k0;
  const long double range = 2.0L * (taxon_count - 1);
  const double threshold = static_cast<double>(
      static_cast<long double>(coefficient) * range / std::sqrt(k) *
      log_ratio / (std::sqrt(log_k) + std::sqrt(log_k0)));
  // Zero is legitimate only at the crossover. A rounded-away negative or
  // positive threshold would silently change the configured decision rule.
  if (!std::isfinite(threshold) || threshold == 0.0 ||
      (threshold > 0.0) != (gene_count > crossover)) {
    throw std::invalid_argument("split certificate threshold is not representable");
  }
  return threshold;
}

HatSplitCertificateResult hat_check_split_certificate(
    const SolverContext& ctx, const std::vector<int>& active_taxa,
    const Bitset& split) {
  HatSplitCertificateResult result;
  if (!ctx.split_certificate_index ||
      !ctx.split_certificate_index->complete_single_copy) return result;
  if (ctx.options.hat_certificate_anchors < 1 ||
      ctx.options.hat_certificate_anchors > 4) {
    throw std::invalid_argument("split certificate requires 1..4 anchor pairs");
  }
  const auto started = std::chrono::steady_clock::now();
  const auto& index = *ctx.split_certificate_index;
  const auto& genes = base_genes(ctx);
  const int k = checked_hdt_int_size(genes.size(), "too many certificate genes");
  const int m = checked_hdt_int_size(active_taxa.size(), "too many active taxa");
  if (index.genes.size() != genes.size() || split.size() != m) {
    throw std::logic_error("split certificate index or split size mismatch");
  }
  int a = -1;
  int b = -1;
  int last_a = -1;
  int last_b = -1;
  for (int i = 0; i < m; ++i) {
    const int taxon = active_taxa[i];
    if (taxon < 0 || taxon >= index.taxon_count) {
      throw std::logic_error("split certificate has an invalid active taxon");
    }
    int& anchor = split.test(i) ? a : b;
    if (anchor < 0) anchor = taxon;
    (split.test(i) ? last_a : last_b) = taxon;
  }
  if (a < 0 || b < 0) {
    throw std::logic_error("split certificate needs two nonempty sides");
  }
  // Integer per-gene rows and a fixed gene-order sum: no rounding before the
  // final division, no atomics, and no dependence on task completion order.
  std::vector<int> rows(checked_hdt_size_mul(
      genes.size(), active_taxa.size(), "split certificate rows overflow"));
  const std::array<std::pair<int, int>, 4> anchor_pairs{{
      {a, b}, {last_a, last_b}, {a, last_b}, {last_a, b}}};
  result.available = true;
  result.gap = -std::numeric_limits<double>::infinity();
  result.threshold = hat_split_certificate_threshold(
      k, index.taxon_count, ctx.options.hat_certificate_coefficient,
      ctx.options.hat_certificate_crossover);
  std::int64_t best_gap_numerator = std::numeric_limits<std::int64_t>::min();
  for (int pair = 0; pair < ctx.options.hat_certificate_anchors; ++pair) {
    const auto anchors = anchor_pairs[pair];
    if (std::find(anchor_pairs.begin(), anchor_pairs.begin() + pair, anchors) !=
        anchor_pairs.begin() + pair) {
      continue;
    }
    parallel_for_inner(0, k, [&](int g) {
      const GeneIndex& gene = genes[g];
      const HatGeneDistanceIndex& distances = index.genes[g];
      const int an = gene.taxon_to_node[anchors.first];
      const int bn = gene.taxon_to_node[anchors.second];
      int* row = rows.data() + static_cast<std::size_t>(g) * m;
      for (int i = 0; i < m; ++i) {
        const int x = gene.taxon_to_node[active_taxa[i]];
        row[i] = distances.distance(gene, an, x) - distances.distance(gene, bn, x);
      }
    }, 1);
    std::vector<std::int64_t> sums(active_taxa.size(), 0);
    for (int g = 0; g < k; ++g) {
      const int* row = rows.data() + static_cast<std::size_t>(g) * m;
      for (int i = 0; i < m; ++i) sums[i] += row[i];
    }
    // |Y| <= |G|-1 <= INT_MAX and k <= INT_MAX, so the signed 64-bit
    // accumulation and the difference of two sums cannot overflow.
    std::int64_t max_a = std::numeric_limits<std::int64_t>::min();
    std::int64_t min_b = std::numeric_limits<std::int64_t>::max();
    for (int i = 0; i < m; ++i) {
      if (split.test(i)) {
        max_a = std::max(max_a, sums[i]);
      } else {
        min_b = std::min(min_b, sums[i]);
      }
    }
    const std::int64_t gap_numerator = min_b - max_a;
    const double gap = static_cast<double>(gap_numerator) / (2.0 * k);
    if (!std::isfinite(gap)) break;
    best_gap_numerator = std::max(best_gap_numerator, gap_numerator);
    result.gap = std::max(result.gap, gap);
    benchmark_counter_add(benchmark_stats().hat_certificate_distance_queries,
                          2LL * k * m);
    // Acceptance by ANY of a fixed number of anchors is safe: the uniform
    // population theorem holds for all anchors, including adaptive choices.
    if (std::isfinite(result.threshold) && result.gap > result.threshold) {
      result.accepted = true;
      break;
    }
  }
  benchmark_counter_add(benchmark_stats().hat_certificate_tests);
  benchmark_counter_add(result.accepted
                            ? benchmark_stats().hat_certificate_accepts
                            : benchmark_stats().hat_certificate_rejects);
  if (best_gap_numerator != std::numeric_limits<std::int64_t>::min() &&
      best_gap_numerator < 0) {
    benchmark_counter_max(
        benchmark_stats().hat_certificate_max_negative_gap_twice_k,
        -best_gap_numerator);
    if (result.accepted) {
      benchmark_counter_add(benchmark_stats().hat_certificate_negative_accepts);
    }
  }
  benchmark_timer_add(benchmark_stats().hat_certificate_ns,
                      std::chrono::steady_clock::now() - started);
  return result;
}

}  // namespace wqfm::internal
