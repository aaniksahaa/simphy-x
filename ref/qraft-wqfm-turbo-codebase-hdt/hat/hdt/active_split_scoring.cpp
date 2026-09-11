// Fixed-palette active-split scoring over immutable source trees.

#include "hat/hat_internal.hpp"
#include "hat/hdt/hdt_checked.hpp"
#include "hat/hdt/hdt_traversal.hpp"
#include "hat/hdt/owner_projection.hpp"

#include <algorithm>
#include <chrono>

namespace wqfm::internal {
namespace {

void score_active_split_cpu_gene_range(
    const std::vector<GeneIndex>& genes,
    const OwnerProjection& owners,
    const HatQueryTree& query,
    const std::vector<unsigned char>& score_node,
    const std::vector<int>& node_to_edge,
    int edge_count,
    int begin_gene,
    int end_gene,
    bool induce_source,
    std::vector<double>& scores);


void score_active_split_cpu_gene_range(
    const std::vector<GeneIndex>& genes,
    const OwnerProjection& owners,
    const HatQueryTree& query,
    const std::vector<unsigned char>& score_node,
    const std::vector<int>& node_to_edge,
    int edge_count,
    int begin_gene,
    int end_gene,
    bool induce_source,
    std::vector<double>& scores) {
  begin_gene = std::max(0, begin_gene);
  end_gene = std::max(begin_gene,
                      std::min(end_gene, static_cast<int>(genes.size())));
  const int cpu_gene_count = end_gene - begin_gene;
  if (cpu_gene_count <= 0 || edge_count <= 0) return;

  const int chunks = hat_score_gene_partition_count(cpu_gene_count);
  const std::size_t partial_count = checked_hdt_size_mul(
      static_cast<std::size_t>(chunks), static_cast<std::size_t>(edge_count),
      "active-split partial score table exceeds the addressable size domain");
  std::vector<double> partial(partial_count, 0.0);
  benchmark_counter_max(
      benchmark_stats().hat_score_partial_peak_bytes,
      checked_hdt_long_long_size(
          checked_hdt_size_mul(
              partial_count, sizeof(double),
              "active-split partial score bytes exceed the addressable size "
              "domain"),
          "active-split partial score bytes exceed the telemetry domain"));
  std::vector<long long> local_nodes(static_cast<std::size_t>(chunks), 0);
  std::vector<long long> local_recolors(static_cast<std::size_t>(chunks), 0);
  std::vector<long long> local_recomputed(static_cast<std::size_t>(chunks), 0);
  std::vector<HdtTraversalStats> local_hdt(
      static_cast<std::size_t>(chunks));
  const auto build_started = std::chrono::steady_clock::now();
  parallel_for_inner(0, chunks, [&](int chunk) {
    const int begin = begin_gene + chunk;
    const int end = begin + 1;
    double* row = partial.data() +
                  static_cast<std::size_t>(chunk) *
                      static_cast<std::size_t>(edge_count);
    long long visits = 0;
    HdtTraversalStats hdt_stats;
    HdtExtractWorkspace hdt_workspace;
    // Active split passes can visit every gene against a large root query.
    // Keep this workspace scoped to the chunk so Taskflow worker threads do not
    // retain large HDT buffers across later recursive calls.
    for (int g = begin; g < end; ++g) {
      const GeneIndex& gene = genes[static_cast<std::size_t>(g)];
      HdtColorPolicy colors;
      colors.default_current = 0;
      colors.default_background = 0;
      colors.current_by_owner.assign(
          static_cast<std::size_t>(owners.owner_count), 1);
      colors.background_by_owner.assign(
          static_cast<std::size_t>(owners.owner_count), 0);
      HdtInitialStateResult initial = hdt_build_initial_state(
          induce_source
              ? build_induced_base_component_graph(gene, owners, 2, colors)
              : build_base_component_graph(gene, owners, 2, colors));
      hdt_account_initial_state(initial, hdt_stats);
      hdt_traverse_active_splits(
          std::move(initial.state), query, query.leaf_active_order, score_node,
          hdt_workspace, hdt_stats,
          [&](int query_node, HdtState& current) {
            const int edge =
                node_to_edge[static_cast<std::size_t>(query_node)];
            if (edge < 0) return;
            const std::array<double, 3> counts =
                hdt_root_counts(current);
            row[static_cast<std::size_t>(edge)] += counts[0] - counts[1];
            ++visits;
          },
          induce_source ? 16 : 4);
    }
    local_nodes[static_cast<std::size_t>(chunk)] = visits;
    local_recolors[static_cast<std::size_t>(chunk)] =
        hdt_stats.recolor_changed_leaves;
    local_recomputed[static_cast<std::size_t>(chunk)] =
        hdt_stats.recolor_marked_nodes;
    local_hdt[static_cast<std::size_t>(chunk)] = hdt_stats;
  }, 1);

  for (int chunk = 0; chunk < chunks; ++chunk) {
    const double* row = partial.data() +
                        static_cast<std::size_t>(chunk) *
                            static_cast<std::size_t>(edge_count);
    for (int e = 0; e < edge_count; ++e) {
      scores[static_cast<std::size_t>(e)] += row[static_cast<std::size_t>(e)];
    }
    benchmark_counter_add(benchmark_stats().hat_hdt_queries,
                          local_nodes[static_cast<std::size_t>(chunk)]);
    benchmark_counter_add(benchmark_stats().hat_hdt_recolor_ops,
                          local_recolors[static_cast<std::size_t>(chunk)]);
    benchmark_counter_add(benchmark_stats().hat_hdt_nodes_recomputed,
                          local_recomputed[static_cast<std::size_t>(chunk)]);
    const HdtTraversalStats& traversal =
        local_hdt[static_cast<std::size_t>(chunk)];
    benchmark_merge_hdt_traversal_stats(traversal);
  }
  benchmark_counter_add(benchmark_stats().hat_hdt_banks_built,
                        static_cast<long long>(cpu_gene_count));
  benchmark_timer_add(
      benchmark_stats().hat_hdt_bank_build_ns,
      std::chrono::duration_cast<std::chrono::nanoseconds>(
          std::chrono::steady_clock::now() - build_started));
}

}  // namespace

std::vector<double> hat_score_active_split_edges_hdt(
    const SolverContext& ctx,
    const HatActiveSet& active,
    const HatQueryTree& query,
    const std::vector<int>& edge_child_nodes) {
  const int edge_count = static_cast<int>(edge_child_nodes.size());
  std::vector<double> scores(static_cast<std::size_t>(edge_count), 0.0);
  if (ctx.options.hat_split_mode == HatSplitMode::ConsensusBalanced ||
      ctx.options.hat_split_mode == HatSplitMode::ConsensusPolytomyBalanced) {
    return scores;
  }
  if (edge_count == 0 || query.root < 0 ||
      static_cast<int>(active.taxa.size()) < 4) {
    return scores;
  }

  const std::vector<GeneIndex>& genes = scoring_genes(ctx);
  const int gene_count = static_cast<int>(genes.size());
  if (gene_count == 0) return scores;
  note_hat_source_bank_scoring_use(genes);
  std::vector<unsigned char> score_node(query.nodes.size(), 0);
  std::vector<int> node_to_edge(query.nodes.size(), -1);
  int valid_edge_count = 0;
  for (int edge = 0; edge < edge_count; ++edge) {
    const int node = edge_child_nodes[static_cast<std::size_t>(edge)];
    if (node < 0 || node >= static_cast<int>(query.nodes.size())) continue;
    score_node[static_cast<std::size_t>(node)] = 1;
    node_to_edge[static_cast<std::size_t>(node)] = edge;
    ++valid_edge_count;
  }
  if (valid_edge_count == 0) return scores;

  const OwnerProjection owners = make_owner_projection(
      ctx, active.active_by_taxon, ctx.active_gene_identity_taxa.get());

  score_active_split_cpu_gene_range(
      genes, owners, query, score_node, node_to_edge, edge_count, 0,
      gene_count,
      ctx.options.hat_one_third_singleton || ctx.split_certificate_index != nullptr,
      scores);
  return scores;
}

}  // namespace wqfm::internal
