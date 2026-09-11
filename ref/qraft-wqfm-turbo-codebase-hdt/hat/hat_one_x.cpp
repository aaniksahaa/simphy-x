#include "hat/hat_internal.hpp"
#include "hat/hat_one_x_internal.hpp"
#include "hat/hdt/hdt_checked.hpp"
#include "hat/hdt/hdt_traversal.hpp"
#include "hat/hdt/owner_projection.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <exception>
#include <limits>
#include <mutex>
#include <string>

namespace wqfm::internal {
namespace {

constexpr double kOneXScoreTieEps = 1e-9;
struct HatAttachView {
  std::vector<int> parent;
  std::vector<int> preorder;
  std::vector<int> postorder;
  std::vector<int> tin;
  std::vector<int> tout;
  std::vector<int> subtree_taxa_count;
  std::vector<int> edge_child;
  std::vector<int> internal_nodes;
  int root = -1;
};

int active_degree(const SolverTree& tree, int node) {
  int degree = 0;
  for (int nb : tree.nodes[static_cast<std::size_t>(node)].neighbors) {
    if (nb >= 0 && nb < static_cast<int>(tree.nodes.size()) &&
        tree.nodes[static_cast<std::size_t>(nb)].active) {
      ++degree;
    }
  }
  return degree;
}

HatAttachView build_hat_attach_view(const SolverTree& tree,
                                    const std::vector<int>& active_taxa,
                                    int taxon_count) {
  HatAttachView view;
  const int n = static_cast<int>(tree.nodes.size());
  view.parent.assign(static_cast<std::size_t>(n), -1);
  view.tin.assign(static_cast<std::size_t>(n), 0);
  view.tout.assign(static_cast<std::size_t>(n), 0);
  view.subtree_taxa_count.assign(static_cast<std::size_t>(n), 0);
  if (n == 0) return view;

  const HatTaxonIndex active_by_taxon =
      make_hat_taxon_index(active_taxa, taxon_count);
  std::vector<int> active_leaf_node(active_by_taxon.taxa.size(), -1);

  for (int i = 0; i < n; ++i) {
    if (!tree.nodes[static_cast<std::size_t>(i)].active) continue;
    const int64_t id = tree.nodes[static_cast<std::size_t>(i)].identity_id;
    if (id >= 0 && id <= std::numeric_limits<int>::max()) {
      const int active = active_by_taxon.find(static_cast<int>(id));
      if (active >= 0) {
        active_leaf_node[static_cast<std::size_t>(active)] = i;
      }
    }
    if (view.root < 0 && id < 0 && active_degree(tree, i) >= 2) {
      view.root = i;
    }
  }
  if (view.root < 0) {
    for (int node : active_leaf_node) {
      if (node >= 0) {
        view.root = node;
        break;
      }
    }
  }
  if (view.root < 0) return view;

  std::vector<std::pair<int, int>> stack{{view.root, 0}};
  view.parent[static_cast<std::size_t>(view.root)] = view.root;
  while (!stack.empty()) {
    const int node = stack.back().first;
    int& next = stack.back().second;
    if (next == 0) {
      view.tin[static_cast<std::size_t>(node)] =
          static_cast<int>(view.preorder.size());
      view.preorder.push_back(node);
    }
    if (next < static_cast<int>(
                   tree.nodes[static_cast<std::size_t>(node)].neighbors.size())) {
      const int nb =
          tree.nodes[static_cast<std::size_t>(node)].neighbors
              [static_cast<std::size_t>(next++)];
      if (!tree.nodes[static_cast<std::size_t>(nb)].active ||
          nb == view.parent[static_cast<std::size_t>(node)]) {
        continue;
      }
      view.parent[static_cast<std::size_t>(nb)] = node;
      view.edge_child.push_back(nb);
      stack.push_back({nb, 0});
      continue;
    }
    view.tout[static_cast<std::size_t>(node)] =
        static_cast<int>(view.preorder.size());
    view.postorder.push_back(node);
    stack.pop_back();
  }

  for (int node : view.postorder) {
    int count = 0;
    const int64_t id = tree.nodes[static_cast<std::size_t>(node)].identity_id;
    if (id >= 0 && id <= std::numeric_limits<int>::max() &&
        active_by_taxon.find(static_cast<int>(id)) >= 0) {
      count = 1;
    }
    for (int nb : tree.nodes[static_cast<std::size_t>(node)].neighbors) {
      if (nb >= 0 && nb < n && tree.nodes[static_cast<std::size_t>(nb)].active &&
          view.parent[static_cast<std::size_t>(nb)] == node) {
        count += view.subtree_taxa_count[static_cast<std::size_t>(nb)];
      }
    }
    view.subtree_taxa_count[static_cast<std::size_t>(node)] = count;
  }

  for (int node : view.preorder) {
    if (tree.nodes[static_cast<std::size_t>(node)].identity_id < 0 &&
        active_degree(tree, node) == 3) {
      view.internal_nodes.push_back(node);
    }
  }
  return view;
}

std::vector<MedianColorPlan> build_median_color_plans(
    const SolverTree& tree,
    const HatAttachView& view,
    const std::vector<int>& active_taxa,
    int taxon_count) {
  (void)active_taxa;
  (void)taxon_count;
  std::vector<MedianColorPlan> plans;
  plans.reserve(view.internal_nodes.size());
  const int total_active =
      view.root >= 0 &&
              view.root < static_cast<int>(view.subtree_taxa_count.size())
          ? view.subtree_taxa_count[static_cast<std::size_t>(view.root)]
          : 0;
  for (int median : view.internal_nodes) {
    MedianColorPlan plan;
    plan.median_node = median;
    const int median_subtree =
        median >= 0 && median < static_cast<int>(view.subtree_taxa_count.size())
            ? view.subtree_taxa_count[static_cast<std::size_t>(median)]
            : 0;
    int direction = 0;
    for (int nb : tree.nodes[static_cast<std::size_t>(median)].neighbors) {
      if (!tree.nodes[static_cast<std::size_t>(nb)].active) continue;
      if (direction >= 3) break;
      int branch_size = 0;
      if (view.parent[static_cast<std::size_t>(nb)] == median) {
        branch_size = view.subtree_taxa_count[static_cast<std::size_t>(nb)];
      } else if (view.parent[static_cast<std::size_t>(median)] == nb) {
        branch_size = total_active - median_subtree;
      }
      plan.branch_size[static_cast<std::size_t>(direction)] = branch_size;
      ++direction;
    }
    plans.push_back(std::move(plan));
  }
  return plans;
}

int direction_for_neighbor(const SolverTree& tree, int median, int neighbor);

OneXScoreWork prepare_one_x_score_work(
    const SolverTree& child_tree,
    const std::vector<int>& active_taxa,
    const std::vector<int>& context_taxa_ids,
    int taxon_count,
    std::vector<MedianColorPlan> plans) {
  OneXScoreWork work;
  work.child_tree = &child_tree;
  work.active_taxa_ids = &active_taxa;
  work.context_taxa_ids = &context_taxa_ids;
  work.query =
      make_hat_query_tree_from_solver_tree(child_tree, active_taxa, taxon_count);
  work.plans = std::move(plans);
  work.totals.assign(work.plans.size(), MedianOneXScores{});
  work.source_to_plan.assign(child_tree.nodes.size(), -1);
  for (int index = 0; index < static_cast<int>(work.plans.size()); ++index) {
    const int source = work.plans[static_cast<std::size_t>(index)].median_node;
    if (source >= 0 && source < static_cast<int>(work.source_to_plan.size())) {
      work.source_to_plan[static_cast<std::size_t>(source)] = index;
    }
  }

  work.direction_by_query_median.assign(
      work.query.nodes.size(), std::array<int, 3>{-1, -1, -1});
  for (int query_median = 0;
       query_median < static_cast<int>(work.query.nodes.size());
       ++query_median) {
    const int source =
        work.query.nodes[static_cast<std::size_t>(query_median)].source_node;
    if (source < 0 || source >= static_cast<int>(work.source_to_plan.size()) ||
        work.source_to_plan[static_cast<std::size_t>(source)] < 0) {
      continue;
    }
    std::array<int, 3> branch_nodes{-1, -1, -1};
    const int child_begin =
        work.query.child_first[static_cast<std::size_t>(query_median)];
    const int child_count =
        work.query.child_count[static_cast<std::size_t>(query_median)];
    if (query_median == work.query.root) {
      if (child_count != 3) continue;
      branch_nodes = std::array<int, 3>{
          work.query.child_nodes[static_cast<std::size_t>(child_begin)],
          work.query.child_nodes[static_cast<std::size_t>(child_begin + 1)],
          work.query.child_nodes[static_cast<std::size_t>(child_begin + 2)]};
    } else {
      if (child_count != 2) continue;
      branch_nodes =
          hat_binary_median_branch_coloring(work.query, query_median)
              .branch_nodes;
    }
    for (int i = 0; i < 3; ++i) {
      const int branch_node = branch_nodes[static_cast<std::size_t>(i)];
      const int branch_source =
          branch_node >= 0 &&
                  branch_node < static_cast<int>(work.query.nodes.size())
              ? work.query.nodes[static_cast<std::size_t>(branch_node)]
                    .source_node
              : -1;
      work.direction_by_query_median[static_cast<std::size_t>(query_median)]
                                    [static_cast<std::size_t>(i)] =
          direction_for_neighbor(child_tree, source, branch_source);
    }
  }

  work.query_active_to_owner.assign(active_taxa.size(), -1);
  for (int active = 0; active < static_cast<int>(active_taxa.size());
       ++active) {
    work.query_active_to_owner[static_cast<std::size_t>(active)] =
        active_taxa[static_cast<std::size_t>(active)];
  }
  work.query_leaf_owner_order =
      query_leaf_owner_order(work.query, work.query_active_to_owner);
  work.context_owner_ids = context_taxa_ids;
  return work;
}

void range_add(std::vector<double>& diff,
               const HatAttachView& view,
               int node,
               double value) {
  if (node < 0 || node >= static_cast<int>(view.tin.size())) return;
  diff[static_cast<std::size_t>(view.tin[static_cast<std::size_t>(node)])] +=
      value;
  diff[static_cast<std::size_t>(view.tout[static_cast<std::size_t>(node)])] -=
      value;
}

int direction_for_neighbor(const SolverTree& tree, int median, int neighbor) {
  int direction = 0;
  for (int nb : tree.nodes[static_cast<std::size_t>(median)].neighbors) {
    if (!tree.nodes[static_cast<std::size_t>(nb)].active) continue;
    if (nb == neighbor) return direction;
    ++direction;
  }
  return -1;
}

void remap_one_x_jobs_to_owners(std::vector<OneXScoreWork>& jobs,
                                const HatTaxonIndex& local_by_taxon) {
  for (OneXScoreWork& job : jobs) {
    if (job.active_taxa_ids != nullptr) {
      const std::vector<int>& active_taxa = *job.active_taxa_ids;
      job.query_active_to_owner.assign(active_taxa.size(), -1);
      for (int active = 0; active < static_cast<int>(active_taxa.size());
           ++active) {
        job.query_active_to_owner[static_cast<std::size_t>(active)] =
            local_by_taxon.find(
                active_taxa[static_cast<std::size_t>(active)]);
      }
      job.query_leaf_owner_order =
          query_leaf_owner_order(job.query,
                                            job.query_active_to_owner);
    }
    job.context_owner_ids.clear();
    if (job.context_taxa_ids != nullptr) {
      job.context_owner_ids.reserve(job.context_taxa_ids->size());
      for (int taxon : *job.context_taxa_ids) {
        const int owner = local_by_taxon.find(taxon);
        if (owner >= 0) job.context_owner_ids.push_back(owner);
      }
    }
  }
}

void add_one_x_scores(MedianOneXScores& total,
                      const MedianOneXScores& part) {
  for (int direction = 0; direction < 3; ++direction) {
    total.raw[static_cast<std::size_t>(direction)] +=
        part.raw[static_cast<std::size_t>(direction)];
  }
}

void score_one_x_jobs_cpu_gene_range(
    const SolverContext& ctx,
    std::vector<OneXScoreWork>& jobs,
    const OwnerProjection& owners,
    int begin_gene,
    int end_gene,
    std::vector<std::vector<MedianOneXScores>>& totals_by_job) {
  if (jobs.empty()) return;
  const std::vector<GeneIndex>& genes = scoring_genes(ctx);
  begin_gene = std::max(0, begin_gene);
  end_gene = std::max(begin_gene,
                      std::min(end_gene, static_cast<int>(genes.size())));
  const int gene_count = end_gene - begin_gene;
  if (gene_count <= 0) return;

  // One row per gene keeps score grouping independent of the worker limit.
  // Work can still be scheduled dynamically, while the final reduction follows
  // stable input-gene order and no fixed chunk count limits useful workers.
  const int chunks = hat_score_gene_partition_count(gene_count);

  std::vector<std::vector<MedianOneXScores>> partials;
  partials.reserve(jobs.size());
  std::size_t partial_count = 0;
  for (const OneXScoreWork& job : jobs) {
    const std::size_t job_partial_count = checked_hdt_size_mul(
        static_cast<std::size_t>(chunks), job.plans.size(),
        "one-X partial score table exceeds the addressable size domain");
    partial_count = checked_hdt_size_add(
        partial_count, job_partial_count,
        "one-X partial score tables exceed the addressable size domain");
    partials.emplace_back(job_partial_count);
  }
  benchmark_counter_max(
      benchmark_stats().hat_score_partial_peak_bytes,
      checked_hdt_long_long_size(
          checked_hdt_size_mul(
              partial_count, sizeof(MedianOneXScores),
              "one-X partial score bytes exceed the addressable size domain"),
          "one-X partial score bytes exceed the telemetry domain"));
  std::vector<long long> local_queries(static_cast<std::size_t>(chunks), 0);
  std::vector<long long> local_recolors(static_cast<std::size_t>(chunks), 0);
  std::vector<long long> local_recomputed(static_cast<std::size_t>(chunks), 0);
  std::vector<HdtTraversalStats> local_hdt(
      static_cast<std::size_t>(chunks));
  parallel_for_dynamic(0, chunks, [&](int chunk) {
    const int begin = begin_gene + chunk;
    const int end = begin + 1;
    long long queries = 0;
    HdtTraversalStats hdt_stats;
    HdtExtractWorkspace hdt_workspace;
    // Keep fixed-vector workspaces bounded to this chunk instead of pinning
    // root-scale buffers in a worker thread-local.
    for (int g = begin; g < end; ++g) {
      const GeneIndex& gene = genes[static_cast<std::size_t>(g)];
      GeneIndex induced_gene;
      const GeneIndex* graph_gene = &gene;
      if (ctx.options.hat_one_third_singleton ||
          ctx.split_certificate_index != nullptr) {
        bool changed = false;
        induced_gene =
            induce_gene_for_owner_projection(gene, owners, changed);
        if (changed) graph_gene = &induced_gene;
      }
      for (int job_index = 0; job_index < static_cast<int>(jobs.size());
           ++job_index) {
        OneXScoreWork& job = jobs[static_cast<std::size_t>(job_index)];
        MedianOneXScores* row =
            partials[static_cast<std::size_t>(job_index)].data() +
            static_cast<std::size_t>(chunk) * job.plans.size();
        auto record_median = [&](int query_median,
                                 const std::array<std::uint8_t, 3>&
                                     branch_colors,
                                 const HatFourColorCounts& counts) {
              const int source =
                  job.query.nodes[static_cast<std::size_t>(query_median)]
                      .source_node;
              if (source < 0 ||
                  source >= static_cast<int>(job.source_to_plan.size())) {
                return;
              }
              const int plan =
                  job.source_to_plan[static_cast<std::size_t>(source)];
              if (plan < 0) return;
              const double q12_03 = counts.q03_12;
              const std::array<double, 3> raw_by_color{
                  q12_03 - 2.0 * (counts.q02_13 + counts.q01_23),
                  counts.q02_13 - 2.0 * (q12_03 + counts.q01_23),
                  counts.q01_23 - 2.0 * (q12_03 + counts.q02_13)};
              for (int i = 0; i < 3; ++i) {
                const int direction =
                    job.direction_by_query_median
                        [static_cast<std::size_t>(query_median)]
                        [static_cast<std::size_t>(i)];
                const int color = branch_colors[static_cast<std::size_t>(i)];
                if (direction >= 0 && direction < 3 && color >= 0 &&
                    color < 3) {
                  MedianOneXScores& score =
                      row[static_cast<std::size_t>(plan)];
                  score.raw[static_cast<std::size_t>(direction)] +=
                      raw_by_color[static_cast<std::size_t>(color)];
                }
              }
              ++queries;
            };
        HdtColorPolicy colors;
        colors.default_current = 0;
        colors.default_background = 0;
        colors.current_by_owner.assign(
            static_cast<std::size_t>(owners.owner_count), 0);
        colors.background_by_owner.assign(
            static_cast<std::size_t>(owners.owner_count), 0);
        for (int owner : job.query_leaf_owner_order) {
          if (owner < 0 || owner >= owners.owner_count) continue;
          colors.current_by_owner[static_cast<std::size_t>(owner)] = 1;
        }
        for (int owner : job.context_owner_ids) {
          if (owner < 0 || owner >= owners.owner_count) continue;
          colors.current_by_owner[static_cast<std::size_t>(owner)] = 3;
          colors.background_by_owner[static_cast<std::size_t>(owner)] = 3;
        }
        HdtInitialStateResult initial = hdt_build_initial_state(
            build_base_component_graph(*graph_gene, owners, 4, colors));
        hdt_account_initial_state(initial, hdt_stats);
        hdt_traverse_one_x_medians(
            std::move(initial.state), job.query, job.query_leaf_owner_order,
            hdt_workspace, hdt_stats,
            [&](int query_median, const std::array<int, 3>&,
                const std::array<std::uint8_t, 3>& branch_colors,
                HdtState& current) {
              const std::array<double, 3> values =
                  hdt_root_counts(current);
              record_median(
                  query_median, branch_colors,
                  HatFourColorCounts{values[0], values[1], values[2]});
            },
            (ctx.options.hat_one_third_singleton ||
             ctx.split_certificate_index != nullptr)
                ? 8 : 4);
      }
    }
    local_queries[static_cast<std::size_t>(chunk)] = queries;
    local_recolors[static_cast<std::size_t>(chunk)] =
        hdt_stats.recolor_changed_leaves;
    local_recomputed[static_cast<std::size_t>(chunk)] =
        hdt_stats.recolor_marked_nodes;
    local_hdt[static_cast<std::size_t>(chunk)] = hdt_stats;
  }, 1);

  long long four_color_passes = 0;
  for (OneXScoreWork& job : jobs) {
    four_color_passes += static_cast<long long>(job.plans.size());
  }
  for (int chunk = 0; chunk < chunks; ++chunk) {
    for (int job_index = 0; job_index < static_cast<int>(jobs.size());
         ++job_index) {
      OneXScoreWork& job = jobs[static_cast<std::size_t>(job_index)];
      const MedianOneXScores* row =
          partials[static_cast<std::size_t>(job_index)].data() +
          static_cast<std::size_t>(chunk) * job.plans.size();
      for (int index = 0; index < static_cast<int>(job.plans.size());
           ++index) {
        add_one_x_scores(
            totals_by_job[static_cast<std::size_t>(job_index)]
                         [static_cast<std::size_t>(index)],
            row[static_cast<std::size_t>(index)]);
      }
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
}

void score_one_x_work_items_hdt(
    const SolverContext& ctx,
    std::vector<OneXScoreWork>& jobs,
    const std::vector<int>& live_taxa_ids,
    int taxon_count,
    const std::vector<int>* identity_taxa) {
  if (jobs.empty()) return;

  const std::vector<GeneIndex>& genes = scoring_genes(ctx);
  note_hat_source_bank_scoring_use(genes);
  const HatTaxonIndex local_by_taxon =
      make_hat_taxon_index(live_taxa_ids, taxon_count);
  const OwnerProjection owners =
      make_owner_projection(ctx, local_by_taxon, identity_taxa);
  remap_one_x_jobs_to_owners(jobs, local_by_taxon);
  std::vector<std::vector<MedianOneXScores>> totals_by_job;
  totals_by_job.reserve(jobs.size());
  for (const OneXScoreWork& job : jobs) {
    totals_by_job.emplace_back(job.plans.size());
  }
  score_one_x_jobs_cpu_gene_range(
      ctx, jobs, owners, 0, static_cast<int>(genes.size()), totals_by_job);

  long long four_color_passes = 0;
  for (int job_index = 0; job_index < static_cast<int>(jobs.size());
       ++job_index) {
    OneXScoreWork& job = jobs[static_cast<std::size_t>(job_index)];
    four_color_passes += static_cast<long long>(job.plans.size());
    for (int index = 0; index < static_cast<int>(job.plans.size()); ++index) {
      add_one_x_scores(job.totals[static_cast<std::size_t>(index)],
                       totals_by_job[static_cast<std::size_t>(job_index)]
                                    [static_cast<std::size_t>(index)]);
    }
  }
  benchmark_counter_add(benchmark_stats().attachment_four_color_passes,
                        four_color_passes);
}

AttachmentEdge singleton_attachment_edge(const HatAttachView& view) {
  AttachmentEdge edge;
  edge.a = view.root;
  edge.b = -1;
  edge.ties = 1;
  edge.second_rank_score = edge.rank_score;
  edge.rank_margin = std::numeric_limits<double>::infinity();
  return edge;
}

struct RankedOneXEdge {
  AttachmentEdge edge;
  int input_order = -1;
};

// Maximization policy for non-finite results: +Inf, finite, -Inf, then NaN.
int one_x_score_class(double score) {
  if (std::isnan(score)) return 0;
  if (std::isinf(score)) return std::signbit(score) ? 1 : 3;
  return 2;
}

bool one_x_stable_edge_less(const RankedOneXEdge& a,
                            const RankedOneXEdge& b) {
  const int a_first = std::min(a.edge.a, a.edge.b);
  const int a_second = std::max(a.edge.a, a.edge.b);
  const int b_first = std::min(b.edge.a, b.edge.b);
  const int b_second = std::max(b.edge.a, b.edge.b);
  if (a_first != b_first) return a_first < b_first;
  if (a_second != b_second) return a_second < b_second;
  if (a.edge.a != b.edge.a) return a.edge.a < b.edge.a;
  if (a.edge.b != b.edge.b) return a.edge.b < b.edge.b;
  return a.input_order < b.input_order;
}

bool one_x_exact_rank_less(const RankedOneXEdge& a,
                           const RankedOneXEdge& b) {
  const int a_class = one_x_score_class(a.edge.rank_score);
  const int b_class = one_x_score_class(b.edge.rank_score);
  if (a_class != b_class) return a_class > b_class;
  if (a_class == 2 && a.edge.rank_score != b.edge.rank_score) {
    return a.edge.rank_score > b.edge.rank_score;
  }
  return one_x_stable_edge_less(a, b);
}

bool one_x_score_in_best_band(double score, double exact_best) {
  const int score_class = one_x_score_class(score);
  const int best_class = one_x_score_class(exact_best);
  if (score_class != best_class) return false;
  if (best_class != 2) return true;
  return std::fabs(score - exact_best) <= kOneXScoreTieEps;
}

double one_x_score_margin(double best, double next) {
  const int best_class = one_x_score_class(best);
  const int next_class = one_x_score_class(next);
  if (best_class != next_class) return std::numeric_limits<double>::infinity();
  if (best_class != 2) return 0.0;
  const double margin = best - next;
  return std::isfinite(margin) ? margin
                               : std::numeric_limits<double>::infinity();
}

AttachmentEdge finish_one_x_edge_from_scores(
    const SolverContext& ctx,
    const SolverTree& child_tree,
    const HatAttachView& view,
    const std::vector<MedianColorPlan>& plans,
    const std::vector<MedianOneXScores>& scores,
    std::size_t context_taxa_count) {
  std::vector<std::array<double, 3>> normalized_coefficients(plans.size());
  for (int index = 0; index < static_cast<int>(plans.size()); ++index) {
    std::array<double, 3> normalized =
        scores[static_cast<std::size_t>(index)].raw;
    const MedianColorPlan& plan = plans[static_cast<std::size_t>(index)];
    const double z =
        static_cast<double>(plan.branch_size[0]) *
        static_cast<double>(plan.branch_size[1]) *
        static_cast<double>(plan.branch_size[2]) *
        static_cast<double>(context_taxa_count);
    if (z > 0.0) {
      const double scale = ctx.options.hat_norm_alpha == 1.0
                               ? 1.0 / z
                               : std::pow(z, -ctx.options.hat_norm_alpha);
      for (double& value : normalized) value *= scale;
      if (ctx.options.hat_one_x_mode == HatOneXMode::Normalized) {
        benchmark_counter_add(
            benchmark_stats().attachment_normalized_one_x_calls);
      }
    }
    normalized_coefficients[static_cast<std::size_t>(index)] = normalized;
  }

  const auto range_started = std::chrono::steady_clock::now();
  std::vector<double> raw_diff(view.preorder.size() + 1, 0.0);
  std::vector<double> normalized_diff(view.preorder.size() + 1, 0.0);
  double raw_global = 0.0;
  double normalized_global = 0.0;
  for (int index = 0; index < static_cast<int>(plans.size()); ++index) {
    const int median = plans[static_cast<std::size_t>(index)].median_node;
    int direction = 0;
    for (int nb : child_tree.nodes[static_cast<std::size_t>(median)].neighbors) {
      if (!child_tree.nodes[static_cast<std::size_t>(nb)].active) continue;
      if (direction >= 3) break;
      const double raw_value =
          scores[static_cast<std::size_t>(index)]
              .raw[static_cast<std::size_t>(direction)];
      const double normalized_value =
          normalized_coefficients[static_cast<std::size_t>(index)]
                                 [static_cast<std::size_t>(direction)];
      ++direction;
      if (view.parent[static_cast<std::size_t>(nb)] == median) {
        range_add(raw_diff, view, nb, raw_value);
        range_add(normalized_diff, view, nb, normalized_value);
      } else if (view.parent[static_cast<std::size_t>(median)] == nb) {
        raw_global += raw_value;
        normalized_global += normalized_value;
        for (int child :
             child_tree.nodes[static_cast<std::size_t>(median)].neighbors) {
          if (child_tree.nodes[static_cast<std::size_t>(child)].active &&
              view.parent[static_cast<std::size_t>(child)] == median) {
            range_add(raw_diff, view, child, -raw_value);
            range_add(normalized_diff, view, child, -normalized_value);
          }
        }
      }
    }
  }

  std::vector<double> raw_prefix(view.preorder.size(), 0.0);
  std::vector<double> normalized_prefix(view.preorder.size(), 0.0);
  double raw_running = 0.0;
  double normalized_running = 0.0;
  for (int i = 0; i < static_cast<int>(view.preorder.size()); ++i) {
    raw_running += raw_diff[static_cast<std::size_t>(i)];
    normalized_running += normalized_diff[static_cast<std::size_t>(i)];
    raw_prefix[static_cast<std::size_t>(i)] = raw_running + raw_global;
    normalized_prefix[static_cast<std::size_t>(i)] =
        normalized_running + normalized_global;
  }

  std::vector<RankedOneXEdge> ranked;
  ranked.reserve(view.edge_child.size());
  for (int child : view.edge_child) {
    AttachmentEdge edge;
    edge.a = view.parent[static_cast<std::size_t>(child)];
    edge.b = child;
    edge.raw_one_x_score =
        raw_prefix[static_cast<std::size_t>(
            view.tin[static_cast<std::size_t>(child)])];
    edge.normalized_one_x_score =
        normalized_prefix[static_cast<std::size_t>(
            view.tin[static_cast<std::size_t>(child)])];
    switch (ctx.options.hat_one_x_mode) {
      case HatOneXMode::Raw:
        edge.rank_score = edge.raw_one_x_score;
        break;
      case HatOneXMode::Normalized:
        edge.rank_score = edge.normalized_one_x_score;
        break;
    }
    ranked.push_back(
        RankedOneXEdge{edge, static_cast<int>(ranked.size())});
  }
  std::sort(ranked.begin(), ranked.end(), one_x_exact_rank_less);
  if (ranked.empty()) return AttachmentEdge{};
  const double exact_best = ranked.front().edge.rank_score;
  int best_ties = 0;
  while (best_ties < static_cast<int>(ranked.size()) &&
         one_x_score_in_best_band(
             ranked[static_cast<std::size_t>(best_ties)].edge.rank_score,
             exact_best)) {
    ++best_ties;
  }
  if (best_ties > 1) {
    benchmark_counter_add(benchmark_stats().attachment_edge_ties,
                          best_ties - 1);
  }
  const auto selected = std::min_element(
      ranked.begin(), ranked.begin() + best_ties, one_x_stable_edge_less);
  AttachmentEdge best = selected->edge;
  best.rank = 0;
  best.ties = best_ties;
  if (best_ties < static_cast<int>(ranked.size())) {
    best.second_rank_score =
        ranked[static_cast<std::size_t>(best_ties)].edge.rank_score;
    best.rank_margin =
        one_x_score_margin(exact_best, best.second_rank_score);
  } else {
    best.second_rank_score = best.rank_score;
    best.rank_margin = std::numeric_limits<double>::infinity();
  }

  benchmark_timer_add(benchmark_stats().attachment_range_add_ns,
                      std::chrono::steady_clock::now() - range_started);
  return best;
}

}  // namespace

AttachmentEdge hat_best_one_x_edge_hdt(
    const SolverContext& ctx,
    const SolverTree& child_tree,
    const std::vector<int>& active_taxa,
    const std::vector<int>& context_taxa_ids,
    const std::vector<int>* identity_taxa) {
  DebugScope scope("hat.one_x");
  const auto total_started = std::chrono::steady_clock::now();
  benchmark_counter_add(benchmark_stats().hat_one_x_calls);

  const int taxon_count = static_cast<int>(context_taxa(ctx).size());
  HatAttachView view =
      build_hat_attach_view(child_tree, active_taxa, taxon_count);
  if (view.root < 0) return AttachmentEdge{};
  if (view.edge_child.empty()) {
    benchmark_timer_add(benchmark_stats().hat_one_x_ns,
                        std::chrono::duration_cast<std::chrono::nanoseconds>(
                            std::chrono::steady_clock::now() - total_started));
    return singleton_attachment_edge(view);
  }
  if (view.internal_nodes.empty()) {
    AttachmentEdge edge;
    const int child = view.edge_child.front();
    edge.a = view.parent[static_cast<std::size_t>(child)];
    edge.b = child;
    edge.ties = 1;
    edge.second_rank_score = edge.rank_score;
    edge.rank_margin = std::numeric_limits<double>::infinity();
    benchmark_timer_add(benchmark_stats().hat_one_x_ns,
                        std::chrono::duration_cast<std::chrono::nanoseconds>(
                            std::chrono::steady_clock::now() - total_started));
    return edge;
  }

  std::vector<MedianColorPlan> plans =
      build_median_color_plans(child_tree, view, active_taxa, taxon_count);
  benchmark_counter_add(benchmark_stats().hat_one_x_medians,
                        static_cast<long long>(plans.size()));
  benchmark_counter_add(benchmark_stats().attachment_internal_nodes,
                        static_cast<long long>(plans.size()));

  const auto score_started = std::chrono::steady_clock::now();
  std::vector<OneXScoreWork> jobs;
  jobs.push_back(prepare_one_x_score_work(child_tree, active_taxa,
                                          context_taxa_ids, taxon_count,
                                          std::move(plans)));
  std::vector<int> live_taxa;
  live_taxa.reserve(active_taxa.size() + context_taxa_ids.size());
  live_taxa.insert(live_taxa.end(), active_taxa.begin(), active_taxa.end());
  live_taxa.insert(live_taxa.end(), context_taxa_ids.begin(),
                   context_taxa_ids.end());
  score_one_x_work_items_hdt(ctx, jobs, live_taxa, taxon_count, identity_taxa);
  benchmark_timer_add(benchmark_stats().attachment_score_ns,
                      std::chrono::steady_clock::now() - score_started);

  const OneXScoreWork& job = jobs.front();
  AttachmentEdge best = finish_one_x_edge_from_scores(
      ctx, child_tree, view, job.plans, job.totals, context_taxa_ids.size());
  benchmark_timer_add(benchmark_stats().hat_one_x_ns,
                      std::chrono::duration_cast<std::chrono::nanoseconds>(
                          std::chrono::steady_clock::now() - total_started));
  return best;
}

std::pair<AttachmentEdge, AttachmentEdge> hat_best_one_x_edges_hdt(
    const SolverContext& ctx,
    const SolverTree& left_tree,
    const SolverTree& right_tree,
    const std::vector<int>& left_taxa,
    const std::vector<int>& right_taxa,
    const std::vector<int>* identity_taxa) {
  DebugScope scope("hat.one_x_pair");
  const auto total_started = std::chrono::steady_clock::now();
  benchmark_counter_add(benchmark_stats().hat_one_x_calls, 2);

  const int taxon_count = static_cast<int>(context_taxa(ctx).size());
  HatAttachView left_view =
      build_hat_attach_view(left_tree, left_taxa, taxon_count);
  HatAttachView right_view =
      build_hat_attach_view(right_tree, right_taxa, taxon_count);

  AttachmentEdge left_edge;
  AttachmentEdge right_edge;
  bool left_needs_score = false;
  bool right_needs_score = false;

  auto initialize_edge = [](const SolverTree& tree,
                            const HatAttachView& view,
                            AttachmentEdge& edge) -> bool {
    (void)tree;
    if (view.root < 0) {
      edge = AttachmentEdge{};
      return false;
    }
    if (view.edge_child.empty()) {
      edge = singleton_attachment_edge(view);
      return false;
    }
    if (view.internal_nodes.empty()) {
      const int child = view.edge_child.front();
      edge.a = view.parent[static_cast<std::size_t>(child)];
      edge.b = child;
      edge.ties = 1;
      edge.second_rank_score = edge.rank_score;
      edge.rank_margin = std::numeric_limits<double>::infinity();
      return false;
    }
    return true;
  };

  left_needs_score = initialize_edge(left_tree, left_view, left_edge);
  right_needs_score = initialize_edge(right_tree, right_view, right_edge);

  std::vector<OneXScoreWork> jobs;
  jobs.reserve(2);
  int left_job = -1;
  int right_job = -1;
  if (left_needs_score) {
    std::vector<MedianColorPlan> plans =
        build_median_color_plans(left_tree, left_view, left_taxa, taxon_count);
    benchmark_counter_add(benchmark_stats().hat_one_x_medians,
                          static_cast<long long>(plans.size()));
    benchmark_counter_add(benchmark_stats().attachment_internal_nodes,
                          static_cast<long long>(plans.size()));
    left_job = static_cast<int>(jobs.size());
    jobs.push_back(prepare_one_x_score_work(left_tree, left_taxa, right_taxa,
                                            taxon_count, std::move(plans)));
  }
  if (right_needs_score) {
    std::vector<MedianColorPlan> plans =
        build_median_color_plans(right_tree, right_view, right_taxa,
                                 taxon_count);
    benchmark_counter_add(benchmark_stats().hat_one_x_medians,
                          static_cast<long long>(plans.size()));
    benchmark_counter_add(benchmark_stats().attachment_internal_nodes,
                          static_cast<long long>(plans.size()));
    right_job = static_cast<int>(jobs.size());
    jobs.push_back(prepare_one_x_score_work(right_tree, right_taxa, left_taxa,
                                            taxon_count, std::move(plans)));
  }

  if (!jobs.empty()) {
    long long median_count = 0;
    for (const OneXScoreWork& job : jobs) {
      median_count += static_cast<long long>(job.plans.size());
    }
    progress_note_stage("hat-one-x L=" + std::to_string(left_taxa.size()) +
                        " R=" + std::to_string(right_taxa.size()) +
                        " med=" + std::to_string(median_count));
    std::vector<int> live_taxa;
    live_taxa.reserve(left_taxa.size() + right_taxa.size());
    live_taxa.insert(live_taxa.end(), left_taxa.begin(), left_taxa.end());
    live_taxa.insert(live_taxa.end(), right_taxa.begin(), right_taxa.end());
    const auto score_started = std::chrono::steady_clock::now();
    score_one_x_work_items_hdt(ctx, jobs, live_taxa, taxon_count,
                               identity_taxa);
    benchmark_timer_add(benchmark_stats().attachment_score_ns,
                        std::chrono::steady_clock::now() - score_started);
  }

  if (left_job >= 0) {
    const OneXScoreWork& job = jobs[static_cast<std::size_t>(left_job)];
    left_edge = finish_one_x_edge_from_scores(
        ctx, left_tree, left_view, job.plans, job.totals, right_taxa.size());
  }
  if (right_job >= 0) {
    const OneXScoreWork& job = jobs[static_cast<std::size_t>(right_job)];
    right_edge = finish_one_x_edge_from_scores(
        ctx, right_tree, right_view, job.plans, job.totals, left_taxa.size());
  }

  benchmark_timer_add(benchmark_stats().hat_one_x_ns,
                      std::chrono::duration_cast<std::chrono::nanoseconds>(
                          std::chrono::steady_clock::now() - total_started));
  return {left_edge, right_edge};
}

}  // namespace wqfm::internal
