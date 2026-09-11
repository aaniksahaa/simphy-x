#pragma once

#include <cmath>
#include <cstddef>
#include <functional>
#include <stdexcept>
#include <vector>

namespace wqfm::internal {

// Maximization order shared with one-X: +Inf, finite, -Inf, then NaN.
enum class HatScoreClass : unsigned char {
  NotANumber = 0,
  NegativeInfinity = 1,
  Finite = 2,
  PositiveInfinity = 3,
};

inline HatScoreClass hat_score_class(double score) noexcept {
  if (std::isnan(score)) return HatScoreClass::NotANumber;
  if (std::isinf(score)) {
    return std::signbit(score) ? HatScoreClass::NegativeInfinity
                               : HatScoreClass::PositiveInfinity;
  }
  return HatScoreClass::Finite;
}

// Strict ordering of score values before any epsilon decision band is applied.
inline bool hat_exact_score_less(double lhs, double rhs) noexcept {
  const HatScoreClass lhs_class = hat_score_class(lhs);
  const HatScoreClass rhs_class = hat_score_class(rhs);
  if (lhs_class != rhs_class) return lhs_class < rhs_class;
  if (lhs_class == HatScoreClass::Finite) return lhs < rhs;
  return false;
}

inline bool hat_score_in_maximum_band(double score,
                                      double exact_maximum,
                                      double epsilon) noexcept {
  const HatScoreClass score_class = hat_score_class(score);
  const HatScoreClass maximum_class = hat_score_class(exact_maximum);
  if (score_class != maximum_class) return false;
  if (maximum_class != HatScoreClass::Finite) return true;
  return std::fabs(score - exact_maximum) <= epsilon;
}

template <typename ForwardIt>
struct HatScoreSelectionResult {
  ForwardIt exact_maximum;
  ForwardIt selected;
  std::size_t band_size = 0;
};

// Selects in two independent passes. The first pass finds the exact numerical
// maximum. The second applies the stable key only within that maximum's epsilon
// band. stable_key_of must include a deterministic final fallback when distinct
// candidates can otherwise have the same key.
template <typename ForwardIt,
          typename ScoreOf,
          typename StableKeyOf,
          typename Eligible>
HatScoreSelectionResult<ForwardIt> hat_select_score_maximum(
    ForwardIt first,
    ForwardIt last,
    double epsilon,
    ScoreOf score_of,
    StableKeyOf stable_key_of,
    Eligible eligible) {
  if (!std::isfinite(epsilon) || epsilon < 0.0) {
    throw std::invalid_argument(
        "HAT score-selection epsilon must be finite and non-negative");
  }

  HatScoreSelectionResult<ForwardIt> result{last, last, 0};
  double exact_maximum = 0.0;
  for (ForwardIt candidate = first; candidate != last; ++candidate) {
    if (!std::invoke(eligible, *candidate)) continue;
    const double score = std::invoke(score_of, *candidate);
    if (result.exact_maximum == last ||
        hat_exact_score_less(exact_maximum, score)) {
      result.exact_maximum = candidate;
      exact_maximum = score;
    }
  }

  if (result.exact_maximum == last) return result;

  for (ForwardIt candidate = first; candidate != last; ++candidate) {
    if (!std::invoke(eligible, *candidate)) continue;
    const double score = std::invoke(score_of, *candidate);
    if (!hat_score_in_maximum_band(score, exact_maximum, epsilon)) continue;
    ++result.band_size;
    if (result.selected == last ||
        std::invoke(stable_key_of, *candidate) <
            std::invoke(stable_key_of, *result.selected)) {
      result.selected = candidate;
    }
  }
  return result;
}

template <typename ForwardIt, typename ScoreOf, typename StableKeyOf>
HatScoreSelectionResult<ForwardIt> hat_select_score_maximum(
    ForwardIt first,
    ForwardIt last,
    double epsilon,
    ScoreOf score_of,
    StableKeyOf stable_key_of) {
  return hat_select_score_maximum(
      first, last, epsilon, score_of, stable_key_of,
      [](const auto&) noexcept { return true; });
}

struct HatNniMove;
struct SolverTree;

void hat_select_nni_paths_for_tests(
    const SolverTree& tree,
    const std::vector<HatNniMove>& candidates,
    double min_gain,
    HatNniMove& cached_result,
    HatNniMove& full_rescore_result);

}  // namespace wqfm::internal
