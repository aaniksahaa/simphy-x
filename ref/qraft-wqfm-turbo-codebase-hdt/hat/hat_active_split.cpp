#include "hat/hat_internal.hpp"
#include "hat/hat_score_selection.hpp"
#include "hat/hat_split_certificate.hpp"

#include <algorithm>
#include <cmath>
#include <iterator>
#include <limits>
#include <tuple>

namespace wqfm::internal {
namespace {

constexpr int kHatMinSide = 1;
constexpr double kHatScoreTieEps = 1e-9;

struct HatConsensusEdgeCandidate {
  struct LeafSet {
    int first_begin = 0;
    int first_end = 0;
    int second_begin = 0;
    int second_end = 0;
    int count = 0;

    bool operator==(const LeafSet& other) const {
      return first_begin == other.first_begin && first_end == other.first_end &&
             second_begin == other.second_begin &&
             second_end == other.second_end && count == other.count;
    }
  };

  LeafSet split;
  Candidate materialized_candidate;
  int edge_child_node = -1;
  int representative_order = std::numeric_limits<int>::max();
  bool materialized = false;
  bool rejected = false;
};

int hat_effective_min_side(const SolverContext& ctx, int active_count) {
  (void)ctx;
  return std::min(kHatMinSide, std::max(1, active_count / 2));
}

int split_balance_from_count(int left, int active_count) {
  return std::min(left, active_count - left);
}

Candidate trivial_completion_candidate(int m) {
  Candidate candidate;
  candidate.split = Bitset(m);
  if (m > 0) candidate.split.set(0);
  candidate.representative_order = std::numeric_limits<int>::max();
  return candidate;
}

void add_leaf_interval(HatConsensusEdgeCandidate::LeafSet& set,
                       int begin,
                       int end) {
  if (begin >= end) return;
  if (set.first_begin == set.first_end) {
    set.first_begin = begin;
    set.first_end = end;
  } else {
    set.second_begin = begin;
    set.second_end = end;
  }
  set.count += end - begin;
}

bool leaf_set_contains(const HatConsensusEdgeCandidate::LeafSet& set,
                       int leaf_pos) {
  return (set.first_begin <= leaf_pos && leaf_pos < set.first_end) ||
         (set.second_begin <= leaf_pos && leaf_pos < set.second_end);
}

bool leaf_set_lex_less(const HatConsensusEdgeCandidate::LeafSet& a,
                       const HatConsensusEdgeCandidate::LeafSet& b,
                       const std::vector<int>& active_to_leaf_pos) {
  for (int active = 0; active < static_cast<int>(active_to_leaf_pos.size());
       ++active) {
    const int leaf_pos = active_to_leaf_pos[static_cast<std::size_t>(active)];
    if (leaf_pos < 0) continue;
    const bool av = leaf_set_contains(a, leaf_pos);
    const bool bv = leaf_set_contains(b, leaf_pos);
    if (av != bv) return !av && bv;
  }
  return false;
}

bool leaf_set_storage_less(const HatConsensusEdgeCandidate::LeafSet& a,
                           const HatConsensusEdgeCandidate::LeafSet& b) {
  return std::tie(a.first_begin, a.first_end, a.second_begin, a.second_end,
                  a.count) <
         std::tie(b.first_begin, b.first_end, b.second_begin, b.second_end,
                  b.count);
}

Bitset materialize_leaf_set(const HatQueryTree& query,
                            const HatConsensusEdgeCandidate::LeafSet& set,
                            int m) {
  Bitset split(m);
  auto add_interval = [&](int begin, int end) {
    begin = std::clamp(begin, 0, static_cast<int>(query.leaf_active_order.size()));
    end = std::clamp(end, begin, static_cast<int>(query.leaf_active_order.size()));
    for (int pos = begin; pos < end; ++pos) {
      const int active = query.leaf_active_order[static_cast<std::size_t>(pos)];
      if (active >= 0 && active < m) split.set(active);
    }
  };
  add_interval(set.first_begin, set.first_end);
  add_interval(set.second_begin, set.second_end);
  return split;
}

void add_query_subtree_to_split(const HatQueryTree& query,
                                int node,
                                int m,
                                Bitset& split) {
  if (node < 0 ||
      node >= static_cast<int>(query.subtree_leaf_first.size()) ||
      node >= static_cast<int>(query.subtree_leaf_count.size())) {
    return;
  }
  const int begin = query.subtree_leaf_first[static_cast<std::size_t>(node)];
  const int end =
      begin + query.subtree_leaf_count[static_cast<std::size_t>(node)];
  for (int pos = begin;
       pos < end && pos < static_cast<int>(query.leaf_active_order.size());
       ++pos) {
    const int active = query.leaf_active_order[static_cast<std::size_t>(pos)];
    if (active >= 0 && active < m) split.set(active);
  }
}

std::vector<HatConsensusEdgeCandidate> polytomy_group_candidates_from_hat_query_tree(
    const HatQueryTree& query,
    int active_count,
    int min_side) {
  std::vector<HatConsensusEdgeCandidate> candidates;
  if (active_count < 4 || query.root < 0) return candidates;
  const int target = active_count / 2;
  if (target < min_side) return candidates;

  for (int node : query.preorder) {
    if (node < 0 || node >= static_cast<int>(query.child_count.size()) ||
        node >= static_cast<int>(query.child_first.size()) ||
        node >= static_cast<int>(query.tin.size())) {
      continue;
    }
    const int children = query.child_count[static_cast<std::size_t>(node)];
    const int min_polytomy_children = node == query.root ? 4 : 3;
    if (children < min_polytomy_children) continue;

    const int child_begin = query.child_first[static_cast<std::size_t>(node)];
    std::vector<int> weights(static_cast<std::size_t>(children), 0);
    for (int i = 0; i < children; ++i) {
      const int child =
          query.child_nodes[static_cast<std::size_t>(child_begin + i)];
      if (child >= 0 &&
          child < static_cast<int>(query.subtree_leaf_count.size())) {
        weights[static_cast<std::size_t>(i)] =
            query.subtree_leaf_count[static_cast<std::size_t>(child)];
      }
    }

    std::vector<unsigned char> reachable(static_cast<std::size_t>(target + 1),
                                         0);
    std::vector<int> previous_sum(static_cast<std::size_t>(target + 1), -1);
    std::vector<int> previous_child(static_cast<std::size_t>(target + 1), -1);
    reachable[0] = 1;
    for (int child_index = 0; child_index < children; ++child_index) {
      const int weight = weights[static_cast<std::size_t>(child_index)];
      if (weight <= 0 || weight > target) continue;
      for (int sum = target - weight; sum >= 0; --sum) {
        if (!reachable[static_cast<std::size_t>(sum)]) continue;
        const int next = sum + weight;
        if (reachable[static_cast<std::size_t>(next)]) continue;
        reachable[static_cast<std::size_t>(next)] = 1;
        previous_sum[static_cast<std::size_t>(next)] = sum;
        previous_child[static_cast<std::size_t>(next)] = child_index;
      }
    }

    int best_sum = -1;
    for (int sum = target; sum >= min_side; --sum) {
      if (reachable[static_cast<std::size_t>(sum)]) {
        best_sum = sum;
        break;
      }
    }
    if (best_sum < 0 || active_count - best_sum < min_side) continue;

    std::vector<unsigned char> selected(static_cast<std::size_t>(children), 0);
    for (int sum = best_sum; sum > 0;) {
      const int child_index = previous_child[static_cast<std::size_t>(sum)];
      const int prev = previous_sum[static_cast<std::size_t>(sum)];
      if (child_index < 0 || prev < 0) break;
      selected[static_cast<std::size_t>(child_index)] = 1;
      sum = prev;
    }

    Candidate grouped;
    grouped.split = Bitset(active_count);
    for (int i = 0; i < children; ++i) {
      if (!selected[static_cast<std::size_t>(i)]) continue;
      const int child =
          query.child_nodes[static_cast<std::size_t>(child_begin + i)];
      add_query_subtree_to_split(query, child, active_count, grouped.split);
    }
    const int count = grouped.split.count();
    if (count < min_side || active_count - count < min_side) continue;
    grouped.representative_order = query.tin[static_cast<std::size_t>(node)];

    HatConsensusEdgeCandidate candidate;
    candidate.materialized = true;
    candidate.materialized_candidate = std::move(grouped);
    candidate.split.count = count;
    candidate.representative_order =
        candidate.materialized_candidate.representative_order;
    candidate.edge_child_node = node;
    candidates.push_back(std::move(candidate));
  }
  return candidates;
}

std::vector<HatConsensusEdgeCandidate> candidates_from_hat_query_tree(
    const HatQueryTree& query) {
  const int m = query.leaf_count;
  if (m < 2 || query.root < 0) return {};

  std::vector<int> active_to_leaf_pos(static_cast<std::size_t>(m), -1);
  for (int pos = 0; pos < static_cast<int>(query.leaf_active_order.size());
       ++pos) {
    const int active = query.leaf_active_order[static_cast<std::size_t>(pos)];
    if (active >= 0 && active < m) {
      active_to_leaf_pos[static_cast<std::size_t>(active)] = pos;
    }
  }
  const int zero_leaf_pos = !active_to_leaf_pos.empty()
                                ? active_to_leaf_pos.front()
                                : -1;

  std::vector<HatConsensusEdgeCandidate> candidates;
  for (int node : query.preorder) {
    if (node == query.root) continue;
    if (node < 0 || node >= static_cast<int>(query.subtree_leaf_first.size())) {
      continue;
    }
    const int begin = query.subtree_leaf_first[static_cast<std::size_t>(node)];
    const int count = query.subtree_leaf_count[static_cast<std::size_t>(node)];
    const int end = begin + count;
    if (count < 2 || count > m - 2) continue;
    HatConsensusEdgeCandidate::LeafSet split;
    if (zero_leaf_pos >= begin && zero_leaf_pos < end) {
      add_leaf_interval(split, 0, begin);
      add_leaf_interval(split, end, m);
    } else {
      add_leaf_interval(split, begin, end);
    }

    const int order = query.tin[static_cast<std::size_t>(node)];
    HatConsensusEdgeCandidate candidate;
    candidate.split = split;
    candidate.representative_order = order;
    candidate.edge_child_node = node;
    candidates.push_back(candidate);
  }

  std::sort(candidates.begin(), candidates.end(),
            [](const HatConsensusEdgeCandidate& a,
               const HatConsensusEdgeCandidate& b) {
              if (leaf_set_storage_less(a.split, b.split)) return true;
              if (leaf_set_storage_less(b.split, a.split)) return false;
              return a.representative_order < b.representative_order;
            });
  candidates.erase(
      std::unique(candidates.begin(), candidates.end(),
                  [](const HatConsensusEdgeCandidate& a,
                     const HatConsensusEdgeCandidate& b) {
                    return a.split == b.split;
                  }),
      candidates.end());
  std::sort(candidates.begin(), candidates.end(),
            [&](const HatConsensusEdgeCandidate& a,
                const HatConsensusEdgeCandidate& b) {
              if (a.representative_order != b.representative_order) {
                return a.representative_order < b.representative_order;
              }
              return leaf_set_lex_less(a.split, b.split, active_to_leaf_pos);
            });
  return candidates;
}

Candidate materialize_candidate(const HatQueryTree& query,
                                const HatConsensusEdgeCandidate& item,
                                int m) {
  if (item.materialized) return item.materialized_candidate;
  Candidate candidate;
  candidate.split = materialize_leaf_set(query, item.split, m);
  candidate.representative_order = item.representative_order;
  return candidate;
}

HatSplitChoice choose_balanced_consensus_split(
    const SolverContext& ctx,
    const HatQueryTree& query,
    const std::vector<HatConsensusEdgeCandidate>& edge_candidates,
    int active_count) {
  HatSplitChoice best;
  best.score = 0.0;
  best.representative_order = std::numeric_limits<int>::max();
  best.candidate_count = static_cast<int>(edge_candidates.size());
  int best_balance = -1;
  int best_index = -1;
  const int min_side = hat_effective_min_side(ctx, active_count);
  for (int i = 0; i < static_cast<int>(edge_candidates.size()); ++i) {
    const HatConsensusEdgeCandidate& candidate =
        edge_candidates[static_cast<std::size_t>(i)];
    const int left_count = candidate.split.count;
    const int right_count = active_count - left_count;
    if (candidate.rejected || left_count < min_side || right_count < min_side) {
      continue;
    }
    const int balance = split_balance_from_count(left_count, active_count);
    const int rep = candidate.representative_order;
    if (balance > best_balance ||
        (balance == best_balance && rep < best.representative_order)) {
      best_balance = balance;
      best.representative_order = rep;
      best_index = i;
    }
  }

  if (best_index >= 0) {
    best.candidate_index = best_index;
    best.candidate = materialize_candidate(
        query, edge_candidates[static_cast<std::size_t>(best_index)],
        active_count);
  } else {
    best.candidate = trivial_completion_candidate(active_count);
    best.representative_order = best.candidate.representative_order;
  }
  const std::vector<Side> assignment =
      assignment_from_split(best.candidate.split, active_count);
  best.left_count =
      static_cast<int>(std::count(assignment.begin(), assignment.end(),
                                  Side::Left));
  best.right_count = static_cast<int>(assignment.size()) - best.left_count;
  return best;
}

int select_scored_active_split_candidate(
    const std::vector<HatConsensusEdgeCandidate>& edge_candidates,
    const std::vector<double>& scores,
    int active_count,
    int min_side) {
  if (scores.size() != edge_candidates.size()) {
    throw std::logic_error("active-split candidate/score size mismatch");
  }

  std::vector<int> eligible;
  eligible.reserve(edge_candidates.size());
  for (int i = 0; i < static_cast<int>(edge_candidates.size()); ++i) {
    const HatConsensusEdgeCandidate& candidate =
        edge_candidates[static_cast<std::size_t>(i)];
    const int left_count = candidate.split.count;
    const int right_count = active_count - left_count;
    if (candidate.rejected || left_count < min_side || right_count < min_side) {
      continue;
    }
    eligible.push_back(i);
  }

  const auto selection = hat_select_score_maximum(
      eligible.begin(), eligible.end(), kHatScoreTieEps,
      [&](int index) { return scores[static_cast<std::size_t>(index)]; },
      [&](int index) {
        const HatConsensusEdgeCandidate& candidate =
            edge_candidates[static_cast<std::size_t>(index)];
        const int balance = split_balance_from_count(candidate.split.count,
                                                     active_count);
        return std::tuple{-balance,
                          candidate.representative_order,
                          candidate.edge_child_node,
                          candidate.split.first_begin,
                          candidate.split.first_end,
                          candidate.split.second_begin,
                          candidate.split.second_end,
                          index};
      });
  return selection.selected == eligible.end() ? -1 : *selection.selected;
}

}  // namespace

int hat_select_active_split_candidate_for_tests(
    const std::vector<HatActiveSplitSelectionProbe>& probes,
    int active_count,
    int min_side) {
  if (active_count < 0 || min_side < 0) {
    throw std::invalid_argument(
        "active-split probe sizes must be non-negative");
  }

  std::vector<HatConsensusEdgeCandidate> candidates;
  std::vector<double> scores;
  candidates.reserve(probes.size());
  scores.reserve(probes.size());
  for (const HatActiveSplitSelectionProbe& probe : probes) {
    HatConsensusEdgeCandidate candidate;
    candidate.split.first_begin = probe.leaf_intervals[0];
    candidate.split.first_end = probe.leaf_intervals[1];
    candidate.split.second_begin = probe.leaf_intervals[2];
    candidate.split.second_end = probe.leaf_intervals[3];
    candidate.split.count = probe.left_count;
    candidate.representative_order = probe.representative_order;
    candidate.edge_child_node = probe.edge_child_node;
    candidates.push_back(candidate);
    scores.push_back(probe.score);
  }
  return select_scored_active_split_candidate(candidates, scores, active_count,
                                              min_side);
}

HatSplitChoice hat_choose_active_split_hdt(
    const SolverContext& ctx,
    const std::vector<Identity>& identities,
    uint64_t path,
    int depth) {
  (void)depth;
  DebugScope scope("hat.choose_active_split");
  const auto started = std::chrono::steady_clock::now();
  progress_note_phase(path, "hat-active-split");

  const HatActiveSet active =
      make_hat_active_set(identities, context_taxa(ctx).size());
  const HatQueryTree query = make_hat_query_tree_from_consensus_induced(
      context_consensus(ctx), active);
  std::vector<HatConsensusEdgeCandidate> edge_candidates =
      candidates_from_hat_query_tree(query);
  if (ctx.options.hat_split_mode == HatSplitMode::ConsensusPolytomyBalanced) {
    std::vector<HatConsensusEdgeCandidate> grouped =
        polytomy_group_candidates_from_hat_query_tree(
            query, static_cast<int>(identities.size()),
            hat_effective_min_side(ctx, static_cast<int>(identities.size())));
    edge_candidates.insert(edge_candidates.end(),
                           std::make_move_iterator(grouped.begin()),
                           std::make_move_iterator(grouped.end()));
  }
  const bool trivial_completion = edge_candidates.empty();
  if (edge_candidates.empty()) {
    if (ctx.options.hat_one_third_singleton) {
      bool complete_permutation =
          query.leaf_count == static_cast<int>(identities.size()) &&
          query.leaf_active_order.size() == identities.size();
      std::vector<std::uint8_t> seen(identities.size(), 0);
      for (int active : query.leaf_active_order) {
        if (active < 0 || active >= static_cast<int>(seen.size()) ||
            seen[static_cast<std::size_t>(active)] != 0) {
          complete_permutation = false;
          break;
        }
        seen[static_cast<std::size_t>(active)] = 1;
      }
      if (!complete_permutation) {
        throw std::logic_error(
            "unresolved HAT query does not contain every active identity");
      }
    }
    HatConsensusEdgeCandidate completion_item;
    completion_item.materialized = true;
    completion_item.materialized_candidate =
        trivial_completion_candidate(static_cast<int>(identities.size()));
    completion_item.split.count =
        completion_item.materialized_candidate.split.count();
    completion_item.representative_order =
        completion_item.materialized_candidate.representative_order;
    edge_candidates.push_back(std::move(completion_item));
  }
  benchmark_counter_add(benchmark_stats().candidate_sets);
  benchmark_counter_add(benchmark_stats().active_consensus_candidate_sets);
  benchmark_counter_add(benchmark_stats().candidates,
                        static_cast<long long>(edge_candidates.size()));
  benchmark_counter_add(benchmark_stats().active_consensus_candidates,
                        static_cast<long long>(edge_candidates.size()));
  benchmark_counter_max(benchmark_stats().max_candidates,
                        static_cast<long long>(edge_candidates.size()));

  if (trivial_completion) {
    HatSplitChoice completion;
    completion.candidate = materialize_candidate(
        query, edge_candidates.front(), static_cast<int>(identities.size()));
    completion.score = 0.0;
    completion.representative_order = completion.candidate.representative_order;
    completion.candidate_count = 1;
    completion.left_count = 1;
    completion.right_count = static_cast<int>(identities.size()) - 1;
    const std::vector<Side> assignment = assignment_from_split(
        completion.candidate.split, static_cast<int>(identities.size()));
    if (!ctx.options.benchmark_json.empty()) {
      benchmark_subproblem_set_initial(
          path, "hat_trivial_completion", completion.candidate_count,
          completion.score, completion.representative_order, assignment,
          identities);
    }
    benchmark_subproblem_add_phase(
        path, "hat_trivial_completion",
        std::chrono::steady_clock::now() - started);
    return completion;
  }

  const bool balanced = ctx.options.hat_split_mode != HatSplitMode::Raw;
  const char* source = balanced ? "hat_consensus_balanced" : "hat_active_hdt";
  const char* phase = balanced ? "hat_consensus_balanced_split" : "hat_active_split";
  if (ctx.options.hat_split_mode == HatSplitMode::ConsensusPolytomyBalanced) {
    source = "hat_consensus_polytomy_balanced";
    phase = "hat_consensus_polytomy_balanced_split";
  }
  std::vector<double> scores;
  if (!balanced) {
    benchmark_counter_add(benchmark_stats().hat_consensus_edges_scored,
                          static_cast<long long>(edge_candidates.size()));
    std::vector<int> edge_child_nodes;
    edge_child_nodes.reserve(edge_candidates.size());
    for (const HatConsensusEdgeCandidate& item : edge_candidates) {
      edge_child_nodes.push_back(item.edge_child_node);
    }
    const auto score_started = std::chrono::steady_clock::now();
    scores = hat_score_active_split_edges_hdt(ctx, active, query, edge_child_nodes);
    benchmark_timer_add(benchmark_stats().hat_split_score_ns,
                        std::chrono::steady_clock::now() - score_started);
  }

  HatSplitChoice best;
  const int active_count = static_cast<int>(identities.size());
  bool accepted = false;
  bool rejected_any = false;
  // The bounded retry count is independent of n. Repeated exact-maximum scans
  // preserve the epsilon-band rule without an epsilon sorting comparator.
  for (int attempt = 0; attempt < ctx.options.hat_certificate_candidates;
       ++attempt) {
    if (balanced) {
      best = choose_balanced_consensus_split(
          ctx, query, edge_candidates, active_count);
    } else {
      best.candidate_index = select_scored_active_split_candidate(
          edge_candidates, scores, active_count,
          hat_effective_min_side(ctx, active_count));
      if (best.candidate_index >= 0) {
        const int i = best.candidate_index;
        best.score = scores[i];
        best.representative_order = edge_candidates[i].representative_order;
        best.candidate = materialize_candidate(
            query, edge_candidates[i], active_count);
      }
    }
    if (best.candidate_index < 0) break;
    const auto certificate = hat_check_split_certificate(
        ctx, active.taxa, best.candidate.split);
    if (!certificate.available || certificate.accepted) {
      accepted = true;
      break;
    }
    edge_candidates[best.candidate_index].rejected = true;
    rejected_any = true;
  }
  best.candidate_count = static_cast<int>(edge_candidates.size());
  if (!accepted) {
    if (ctx.options.hat_one_third_singleton) {
      throw std::logic_error(
          "HAT one-third/singleton mode found no eligible active split");
    }
    best.candidate = trivial_completion_candidate(active_count);
    best.candidate_index = -1;
    best.score = 0.0;
    best.representative_order = best.candidate.representative_order;
    if (rejected_any) {
      benchmark_counter_add(benchmark_stats().hat_certificate_singletons);
      source = "hat_certificate_singleton";
    }
  }
  const std::vector<Side> assignment =
      assignment_from_split(best.candidate.split,
                            static_cast<int>(identities.size()));
  best.left_count =
      static_cast<int>(std::count(assignment.begin(), assignment.end(),
                                  Side::Left));
  best.right_count = static_cast<int>(assignment.size()) - best.left_count;

  if (!ctx.options.benchmark_json.empty()) {
    benchmark_subproblem_set_initial(path, source,
                                     best.candidate_count, best.score,
                                     best.representative_order, assignment,
                                     identities);
  }
  benchmark_subproblem_add_phase(path, phase,
                                 std::chrono::steady_clock::now() - started);
  return best;
}

}  // namespace wqfm::internal
