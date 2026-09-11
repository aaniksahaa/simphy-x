#pragma once

#include <algorithm>
#include <array>
#include <cstdint>
#include <span>
#include <vector>

// Fixed-palette boundary-DP payload: the MINIMAL observable-quotient coordinates
// for weighted-quartet scoring (see Vector_formulation.md and
// scripts/complexity/gen_boundary_dp.py). Each HDT component carries a fixed-size
// coordinate vector `q = P_T . Sigma(component)`; every grammar operation is a
// fixed (bi)linear tensor contraction, so each merge is O(1) for fixed palette K.
//
// This module exposes the SCORING (center-fold) grammar operations
//   LEAF_TO_C, C_TO_G, GG_TO_G, IG_TF (G->C1, non-root), IG_FF (G->Closed, root)
// plus the root readout, driven by the exact-rational tensors emitted into
// src/hat/hdt/generated/boundary_dp_scoring_k{2,4}.hpp. Palette sizes:
//   K = 2  -> active-split objective   (2 channels: q00_11, q01_01)
//   K = 4  -> one-X / NNI objective    (3 channels: q01_23, q02_13, q03_12)
//
// Interfaces and their quotient dimensions (K=2 / K=4):
//   Closed  2 / 3   (root readout coords)
//   C1     12 / 26  (subtree by one up-port)
//   G      15 / 32  (star center)
namespace wqfm {

struct BoundaryTensors {
  int palette = 0;  // K
  int dim_closed = 0, dim_c1 = 0, dim_g = 0, n_channels = 0;
  std::vector<double> c_to_g;        // [dim_g * dim_c1]        C1 -> G
  std::vector<double> gg;            // [dim_g * dim_g * dim_g] (G,G) -> G
  std::vector<double> ig_ff;         // [dim_closed * dim_g]    G -> Closed
  std::vector<double> ig_tf;         // [dim_c1 * dim_g]        G -> C1
  std::vector<double> readout;       // [n_channels * dim_closed]
  std::vector<double> leaf_absent;   // [dim_c1]
  std::vector<double> leaf_present;  // [K * dim_c1]  (row c = a color-c leaf, weight 1)
};

// K must be 2 or 4. Loaded once (thread-safe static init); returns an empty
// tensor set for any other K.
const BoundaryTensors& boundary_tensors(int palette);

using BdpVec = std::vector<double>;

// Typed scoring-grammar operations on quotient-coordinate payloads.
BdpVec bdp_leaf(const BoundaryTensors& t, int color, double weight);     // -> C1
BdpVec bdp_c_to_g(const BoundaryTensors& t, const BdpVec& c1);           // C1 -> G
BdpVec bdp_gg(const BoundaryTensors& t, const BdpVec& a, const BdpVec& b);  // (G,G) -> G
BdpVec bdp_ig_tf(const BoundaryTensors& t, const BdpVec& g);             // G -> C1
BdpVec bdp_ig_ff(const BoundaryTensors& t, const BdpVec& g);             // G -> Closed

// Readout: the (<=3) resolved-quartet channel totals from a Closed payload.
std::array<double, 3> bdp_readout(const BoundaryTensors& t, const BdpVec& closed);

// Complete Contract/rebuild algebra. Unlike the scoring-only API above, this
// includes both C1 orientations, C2 corridors, CLOSE, REVERSE, and CC_TO_C.
// The largest generated quotient has 56 coordinates, so values are stored
// inline and every transition is allocation-free.
enum class FullBdpInterface : std::uint8_t {
  Closed,
  C1,
  C2,
  G,
};

enum class FullBdpOrientation : std::uint8_t {
  None,
  U,
  D,
};

inline constexpr std::size_t kFullBdpMaxDimension = 56;

struct FullBdpVec {
  FullBdpInterface interface = FullBdpInterface::Closed;
  FullBdpOrientation orientation = FullBdpOrientation::None;
  std::uint8_t dimension = 0;
  std::array<double, kFullBdpMaxDimension> values;

  FullBdpVec() = default;

  FullBdpVec(const FullBdpVec& other)
      : interface(other.interface),
        orientation(other.orientation),
        dimension(other.dimension) {
    std::copy_n(other.values.data(), dimension, values.data());
  }

  FullBdpVec& operator=(const FullBdpVec& other) {
    if (this == &other) return *this;
    interface = other.interface;
    orientation = other.orientation;
    dimension = other.dimension;
    std::copy_n(other.values.data(), dimension, values.data());
    return *this;
  }

  FullBdpVec(FullBdpVec&& other) noexcept
      : interface(other.interface),
        orientation(other.orientation),
        dimension(other.dimension) {
    std::copy_n(other.values.data(), dimension, values.data());
  }

  FullBdpVec& operator=(FullBdpVec&& other) noexcept {
    if (this == &other) return *this;
    interface = other.interface;
    orientation = other.orientation;
    dimension = other.dimension;
    std::copy_n(other.values.data(), dimension, values.data());
    return *this;
  }

  std::span<double> active() {
    return std::span<double>(values.data(), dimension);
  }
  std::span<const double> active() const {
    return std::span<const double>(values.data(), dimension);
  }
};

struct FullBdpDual {
  FullBdpVec current;
  FullBdpVec background;
  bool current_equals_background = false;
};

struct FullBdpConstView {
  FullBdpInterface interface = FullBdpInterface::Closed;
  FullBdpOrientation orientation = FullBdpOrientation::None;
  std::uint8_t dimension = 0;
  const double* values = nullptr;
};

struct FullBdpMutableView {
  FullBdpInterface interface = FullBdpInterface::Closed;
  FullBdpOrientation orientation = FullBdpOrientation::None;
  std::uint8_t dimension = 0;
  double* values = nullptr;
};

inline FullBdpConstView full_bdp_view(const FullBdpVec& value) {
  return FullBdpConstView{value.interface, value.orientation, value.dimension,
                          value.values.data()};
}

inline FullBdpMutableView full_bdp_view(FullBdpVec& value) {
  return FullBdpMutableView{value.interface, value.orientation,
                            value.dimension, value.values.data()};
}

int full_bdp_dimension(int palette, FullBdpInterface interface);
int full_bdp_max_dimension(int palette);
int full_bdp_channel_count(int palette);
long long full_bdp_thread_scalar_fmas();

FullBdpVec full_bdp_leaf(int palette, int color, double weight);
FullBdpDual full_bdp_leaf_dual(int palette,
                               int current_color,
                               int background_color,
                               double weight);
FullBdpVec full_bdp_empty_center(int palette);
FullBdpVec full_bdp_identity_c2(int palette);
FullBdpVec full_bdp_empty_closed(int palette);

FullBdpVec full_bdp_c_to_g(int palette, const FullBdpVec& c1_u);
FullBdpVec full_bdp_ig_ff(int palette, const FullBdpVec& g);
FullBdpVec full_bdp_ig_tf(int palette, const FullBdpVec& g);
FullBdpVec full_bdp_ig_tt(int palette, const FullBdpVec& g);
FullBdpVec full_bdp_close(int palette, const FullBdpVec& c1_u);
FullBdpVec full_bdp_reverse(int palette, const FullBdpVec& c);
FullBdpVec full_bdp_gg(int palette,
                       const FullBdpVec& left,
                       const FullBdpVec& right);
FullBdpVec full_bdp_cc(int palette,
                       const FullBdpVec& higher,
                       const FullBdpVec& lower);

std::array<double, 3> full_bdp_readout(int palette,
                                       const FullBdpVec& closed);

void full_bdp_leaf_into(int palette,
                        int color,
                        double weight,
                        FullBdpMutableView output);
void full_bdp_c_to_g_into(int palette,
                          FullBdpConstView input,
                          FullBdpMutableView output);
void full_bdp_ig_ff_into(int palette,
                         FullBdpConstView input,
                         FullBdpMutableView output);
void full_bdp_ig_tf_into(int palette,
                         FullBdpConstView input,
                         FullBdpMutableView output);
void full_bdp_ig_tt_into(int palette,
                         FullBdpConstView input,
                         FullBdpMutableView output);
void full_bdp_close_into(int palette,
                         FullBdpConstView input,
                         FullBdpMutableView output);
void full_bdp_reverse_into(int palette,
                           FullBdpConstView input,
                           FullBdpMutableView output);
void full_bdp_gg_into(int palette,
                      FullBdpConstView left,
                      FullBdpConstView right,
                      FullBdpMutableView output);
void full_bdp_cc_into(int palette,
                      FullBdpConstView higher,
                      FullBdpConstView lower,
                      FullBdpMutableView output);
void full_bdp_c_to_g_background_into(int palette,
                                     FullBdpConstView input,
                                     FullBdpMutableView output);
void full_bdp_ig_ff_background_into(int palette,
                                    FullBdpConstView input,
                                    FullBdpMutableView output);
void full_bdp_ig_tf_background_into(int palette,
                                    FullBdpConstView input,
                                    FullBdpMutableView output);
void full_bdp_ig_tt_background_into(int palette,
                                    FullBdpConstView input,
                                    FullBdpMutableView output);
void full_bdp_close_background_into(int palette,
                                    FullBdpConstView input,
                                    FullBdpMutableView output);
void full_bdp_reverse_background_into(int palette,
                                      FullBdpConstView input,
                                      FullBdpMutableView output);
void full_bdp_gg_background_into(int palette,
                                 FullBdpConstView left,
                                 FullBdpConstView right,
                                 FullBdpMutableView output);
void full_bdp_cc_background_into(int palette,
                                 FullBdpConstView higher,
                                 FullBdpConstView lower,
                                 FullBdpMutableView output);

// Exact color-graded kernels used by the HDT. Active coordinates are
// packed in increasing full-coordinate order. Packing changes storage only;
// every generated scalar reduction retains its original term order.
bool full_bdp_coordinate_active(int palette,
                                FullBdpInterface interface,
                                int coordinate,
                                std::uint8_t color_mask);
std::uint64_t full_bdp_active_coordinate_bits(
    int palette,
    FullBdpInterface interface,
    std::uint8_t color_mask);
int full_bdp_active_coordinate_count(int palette,
                                     FullBdpInterface interface,
                                     std::uint8_t color_mask);
int full_bdp_packed_coordinate_rank(int palette,
                                    FullBdpInterface interface,
                                    std::uint8_t color_mask,
                                    int coordinate);
std::uint8_t full_bdp_value_color_mask(int palette,
                                       const FullBdpVec& value);
void full_bdp_leaf_color_mask_into(int palette,
                                   int color,
                                   double weight,
                                   std::uint8_t color_mask,
                                   double* output);
void full_bdp_c_to_g_color_mask_into(int palette,
                                     std::uint8_t color_mask,
                                     const double* input,
                                     double* output);
void full_bdp_ig_ff_color_mask_into(int palette,
                                    std::uint8_t color_mask,
                                    const double* input,
                                    double* output);
void full_bdp_ig_tf_color_mask_into(int palette,
                                    std::uint8_t color_mask,
                                    const double* input,
                                    double* output);
void full_bdp_ig_tt_color_mask_into(int palette,
                                    std::uint8_t color_mask,
                                    const double* input,
                                    double* output);
void full_bdp_close_color_mask_into(int palette,
                                    std::uint8_t color_mask,
                                    const double* input,
                                    double* output);
void full_bdp_reverse_c1_color_mask_into(int palette,
                                         std::uint8_t color_mask,
                                         const double* input,
                                         double* output);
void full_bdp_reverse_c2_color_mask_into(int palette,
                                         std::uint8_t color_mask,
                                         const double* input,
                                         double* output);
void full_bdp_gg_color_mask_into(int palette,
                                 std::uint8_t left_mask,
                                 std::uint8_t right_mask,
                                 const double* left,
                                 const double* right,
                                 double* output);
void full_bdp_cc_color_mask_into(int palette,
                                 FullBdpInterface higher_interface,
                                 FullBdpInterface lower_interface,
                                 std::uint8_t higher_mask,
                                 std::uint8_t lower_mask,
                                 const double* higher,
                                 const double* lower,
                                 double* output);
std::array<double, 3> full_bdp_color_mask_readout(
    int palette,
    std::uint8_t color_mask,
    const double* closed);
std::array<double, 3> full_bdp_readout(int palette,
                                       FullBdpConstView closed);

FullBdpDual full_bdp_c_to_g(int palette, const FullBdpDual& c1_u);
FullBdpDual full_bdp_ig_ff(int palette, const FullBdpDual& g);
FullBdpDual full_bdp_ig_tf(int palette, const FullBdpDual& g);
FullBdpDual full_bdp_ig_tt(int palette, const FullBdpDual& g);
FullBdpDual full_bdp_close(int palette, const FullBdpDual& c1_u);
FullBdpDual full_bdp_reverse(int palette, const FullBdpDual& c);
FullBdpDual full_bdp_gg(int palette,
                        const FullBdpDual& left,
                        const FullBdpDual& right);
FullBdpDual full_bdp_cc(int palette,
                        const FullBdpDual& higher,
                        const FullBdpDual& lower);

}  // namespace wqfm
