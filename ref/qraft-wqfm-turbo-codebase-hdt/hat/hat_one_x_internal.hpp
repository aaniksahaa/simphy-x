#pragma once

#include "hat/hat_internal.hpp"

#include <array>

namespace wqfm::internal {

struct MedianColorPlan {
  int median_node = -1;
  std::array<int, 3> branch_size{0, 0, 0};
};

struct MedianOneXScores {
  std::array<double, 3> raw{0.0, 0.0, 0.0};
};

struct OneXScoreWork {
  const SolverTree* child_tree = nullptr;
  const std::vector<int>* active_taxa_ids = nullptr;
  const std::vector<int>* context_taxa_ids = nullptr;
  HatQueryTree query;
  std::vector<MedianColorPlan> plans;
  std::vector<int> source_to_plan;
  std::vector<std::array<int, 3>> direction_by_query_median;
  std::vector<int> query_active_to_owner;
  std::vector<int> query_leaf_owner_order;
  std::vector<int> context_owner_ids;
  std::vector<MedianOneXScores> totals;
};

}  // namespace wqfm::internal
