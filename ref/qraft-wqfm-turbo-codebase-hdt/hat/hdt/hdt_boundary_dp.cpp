#include "hat/hdt/hdt_boundary_dp.hpp"

#include "generated/boundary_dp_scoring_k2.hpp"
#include "generated/boundary_dp_scoring_k4.hpp"
#include "generated/boundary_dp_k2.hpp"
#include "generated/boundary_dp_k4.hpp"

#include <bit>
#include <stdexcept>
#include <string>
#include <utility>

namespace wqfm {

namespace {

thread_local long long g_full_bdp_scalar_fmas = 0;

// Load a generated namespace's exact-rational tensors as doubles.
#define WQFM_BUILD_BOUNDARY_TENSORS(NS, K)                                     \
  [] {                                                                         \
    namespace n = NS;                                                          \
    BoundaryTensors t;                                                         \
    t.palette = (K);                                                           \
    t.dim_closed = n::dim_Closed;                                              \
    t.dim_c1 = n::dim_C1;                                                      \
    t.dim_g = n::dim_G;                                                        \
    t.n_channels = n::n_channels;                                             \
    auto cp = [](std::vector<double>& d, const int64_t* num,                   \
                 const int64_t* den, int c) {                                  \
      d.resize(static_cast<std::size_t>(c));                                   \
      for (int i = 0; i < c; ++i)                                              \
        d[static_cast<std::size_t>(i)] =                                       \
            static_cast<double>(num[i]) / static_cast<double>(den[i]);         \
    };                                                                         \
    cp(t.c_to_g, n::C_TO_G_num, n::C_TO_G_den, t.dim_g * t.dim_c1);            \
    cp(t.gg, n::GG_TO_G_num, n::GG_TO_G_den, t.dim_g * t.dim_g * t.dim_g);     \
    cp(t.ig_ff, n::IG_FF_num, n::IG_FF_den, t.dim_closed * t.dim_g);           \
    cp(t.ig_tf, n::IG_TF_num, n::IG_TF_den, t.dim_c1 * t.dim_g);               \
    cp(t.readout, n::readout_num, n::readout_den,                             \
       t.n_channels * t.dim_closed);                                          \
    cp(t.leaf_absent, n::leaf_absent_num, n::leaf_absent_den, t.dim_c1);       \
    cp(t.leaf_present, n::leaf_present_num, n::leaf_present_den,               \
       (K) * t.dim_c1);                                                        \
    return t;                                                                  \
  }()

// M[rout x rin] applied to x[rin] -> y[rout].
BdpVec unary(const std::vector<double>& M, int rout, int rin, const BdpVec& x) {
  BdpVec y(static_cast<std::size_t>(rout), 0.0);
  for (int k = 0; k < rout; ++k) {
    double s = 0.0;
    const int base = k * rin;
    for (int p = 0; p < rin; ++p) s += M[static_cast<std::size_t>(base + p)] * x[static_cast<std::size_t>(p)];
    y[static_cast<std::size_t>(k)] = s;
  }
  return y;
}

FullBdpVec make_full_descriptor(int palette,
                                FullBdpInterface interface,
                                FullBdpOrientation orientation) {
  FullBdpVec out;
  out.interface = interface;
  out.orientation = orientation;
  out.dimension = static_cast<std::uint8_t>(
      full_bdp_dimension(palette, interface));
  if (out.dimension == 0) {
    throw std::invalid_argument("unsupported boundary-DP palette");
  }
  return out;
}

FullBdpVec make_full_value(int palette,
                           FullBdpInterface interface,
                           FullBdpOrientation orientation) {
  FullBdpVec out = make_full_descriptor(palette, interface, orientation);
  std::fill_n(out.values.data(), out.dimension, 0.0);
  return out;
}

void require_full_value(int palette,
                        const FullBdpVec& value,
                        FullBdpInterface interface,
                        FullBdpOrientation orientation,
                        const char* operation) {
  const int expected = full_bdp_dimension(palette, interface);
  if (expected <= 0 || value.interface != interface ||
      value.orientation != orientation || value.dimension != expected) {
    throw std::invalid_argument(std::string(operation) +
                                ": invalid boundary-DP interface");
  }
}

void require_full_view(int palette,
                       FullBdpConstView value,
                       FullBdpInterface interface,
                       FullBdpOrientation orientation,
                       const char* operation) {
  const int expected = full_bdp_dimension(palette, interface);
  if (expected <= 0 || value.interface != interface ||
      value.orientation != orientation || value.dimension != expected ||
      value.values == nullptr) {
    throw std::invalid_argument(std::string(operation) +
                                ": invalid boundary-DP view");
  }
}

[[maybe_unused]] void require_full_view(
    int palette,
    FullBdpMutableView value,
    FullBdpInterface interface,
    FullBdpOrientation orientation,
    const char* operation) {
  require_full_view(
      palette,
      FullBdpConstView{value.interface, value.orientation, value.dimension,
                       value.values},
      interface, orientation, operation);
}

void require_full_hot_view(int palette,
                           FullBdpConstView value,
                           FullBdpInterface interface,
                           FullBdpOrientation orientation,
                           const char* operation) {
#ifndef NDEBUG
  require_full_view(palette, value, interface, orientation, operation);
#else
  (void)palette;
  (void)value;
  (void)interface;
  (void)orientation;
  (void)operation;
#endif
}

void require_full_hot_view(int palette,
                           FullBdpMutableView value,
                           FullBdpInterface interface,
                           FullBdpOrientation orientation,
                           const char* operation) {
#ifndef NDEBUG
  require_full_view(palette, value, interface, orientation, operation);
#else
  (void)palette;
  (void)value;
  (void)interface;
  (void)orientation;
  (void)operation;
#endif
}

template <auto Eval, class Terms>
void generated_full_unary_into(FullBdpConstView input,
                               FullBdpMutableView output,
                               const Terms& terms) {
  g_full_bdp_scalar_fmas += static_cast<long long>(terms.size());
  Eval(input.values, output.values);
}

template <auto Eval, class Terms>
void generated_full_binary_into(FullBdpConstView left,
                                FullBdpConstView right,
                                FullBdpMutableView output,
                                const Terms& terms) {
  g_full_bdp_scalar_fmas += static_cast<long long>(terms.size());
  Eval(left.values, right.values, output.values);
}

template <class Terms>
void sparse_full_constant_into(FullBdpMutableView output,
                               const Terms& terms) {
  std::fill_n(output.values, output.dimension, 0.0);
  g_full_bdp_scalar_fmas += static_cast<long long>(terms.size());
  for (const auto& term : terms) {
    output.values[term.index] = term.sign > 0 ? 1.0 : -1.0;
  }
}

template <auto Eval, class Terms>
FullBdpVec generated_full_unary(int palette,
                                const FullBdpVec& input,
                                FullBdpInterface output_interface,
                                FullBdpOrientation output_orientation,
                                const Terms& terms) {
  FullBdpVec output =
      make_full_descriptor(palette, output_interface, output_orientation);
  g_full_bdp_scalar_fmas += static_cast<long long>(terms.size());
  Eval(input.values.data(), output.values.data());
  return output;
}

template <auto Eval, class Terms>
FullBdpVec generated_full_binary(int palette,
                                 const FullBdpVec& left,
                                 const FullBdpVec& right,
                                 FullBdpInterface output_interface,
                                 FullBdpOrientation output_orientation,
                                 const Terms& terms) {
  FullBdpVec output =
      make_full_descriptor(palette, output_interface, output_orientation);
  g_full_bdp_scalar_fmas += static_cast<long long>(terms.size());
  Eval(left.values.data(), right.values.data(), output.values.data());
  return output;
}

template <auto Eval, class Terms>
FullBdpDual generated_full_unary_dual(
    int palette,
    const FullBdpDual& input,
    FullBdpInterface output_interface,
    FullBdpOrientation output_orientation,
    const Terms& terms) {
  FullBdpDual output;
  output.current =
      make_full_descriptor(palette, output_interface, output_orientation);
  output.background =
      make_full_descriptor(palette, output_interface, output_orientation);
  output.current_equals_background = false;
  g_full_bdp_scalar_fmas +=
      2 * static_cast<long long>(terms.size());
  Eval(input.current.values.data(), input.background.values.data(),
       output.current.values.data(), output.background.values.data());
  return output;
}

template <auto Eval, class Terms>
FullBdpDual generated_full_binary_dual(
    int palette,
    const FullBdpDual& left,
    const FullBdpDual& right,
    FullBdpInterface output_interface,
    FullBdpOrientation output_orientation,
    const Terms& terms) {
  FullBdpDual output;
  output.current =
      make_full_descriptor(palette, output_interface, output_orientation);
  output.background =
      make_full_descriptor(palette, output_interface, output_orientation);
  output.current_equals_background = false;
  g_full_bdp_scalar_fmas +=
      2 * static_cast<long long>(terms.size());
  Eval(left.current.values.data(), left.background.values.data(),
       right.current.values.data(), right.background.values.data(),
       output.current.values.data(), output.background.values.data());
  return output;
}

template <class Terms>
FullBdpVec sparse_full_constant(int palette,
                                FullBdpInterface interface,
                                FullBdpOrientation orientation,
                                const Terms& terms) {
  FullBdpVec output = make_full_value(palette, interface, orientation);
  g_full_bdp_scalar_fmas += static_cast<long long>(terms.size());
  for (const auto& term : terms) {
    output.values[term.index] = term.sign > 0 ? 1.0 : -1.0;
  }
  return output;
}

template <class Terms>
std::array<double, 3> sparse_full_readout(const FullBdpVec& closed,
                                          const Terms& terms) {
  std::array<double, 3> output{0.0, 0.0, 0.0};
  g_full_bdp_scalar_fmas += static_cast<long long>(terms.size());
  if (terms.empty()) return output;
  std::uint8_t active_out = terms.front().out;
  double sum = 0.0;
  for (const auto& term : terms) {
    if (term.out != active_out) {
      output[active_out] += sum;
      active_out = term.out;
      sum = 0.0;
    }
    const double value = closed.values[term.in];
    sum += term.sign > 0 ? value : -value;
  }
  output[active_out] += sum;
  return output;
}

}  // namespace

const BoundaryTensors& boundary_tensors(int palette) {
  static const BoundaryTensors k2 = WQFM_BUILD_BOUNDARY_TENSORS(boundary_dp_k2, 2);
  static const BoundaryTensors k4 = WQFM_BUILD_BOUNDARY_TENSORS(boundary_dp_k4, 4);
  static const BoundaryTensors empty{};
  if (palette == 2) return k2;
  if (palette == 4) return k4;
  return empty;
}

BdpVec bdp_leaf(const BoundaryTensors& t, int color, double weight) {
  BdpVec q(static_cast<std::size_t>(t.dim_c1));
  const int base = color * t.dim_c1;
  for (int i = 0; i < t.dim_c1; ++i)
    q[static_cast<std::size_t>(i)] =
        t.leaf_absent[static_cast<std::size_t>(i)] +
        weight * t.leaf_present[static_cast<std::size_t>(base + i)];
  return q;
}

BdpVec bdp_c_to_g(const BoundaryTensors& t, const BdpVec& c1) {
  return unary(t.c_to_g, t.dim_g, t.dim_c1, c1);
}

BdpVec bdp_gg(const BoundaryTensors& t, const BdpVec& a, const BdpVec& b) {
  const int d = t.dim_g;
  BdpVec y(static_cast<std::size_t>(d), 0.0);
  for (int k = 0; k < d; ++k) {
    double s = 0.0;
    for (int p = 0; p < d; ++p) {
      const double ap = a[static_cast<std::size_t>(p)];
      if (ap == 0.0) continue;
      const int base = (k * d + p) * d;
      for (int q = 0; q < d; ++q)
        s += t.gg[static_cast<std::size_t>(base + q)] * ap * b[static_cast<std::size_t>(q)];
    }
    y[static_cast<std::size_t>(k)] = s;
  }
  return y;
}

BdpVec bdp_ig_tf(const BoundaryTensors& t, const BdpVec& g) {
  return unary(t.ig_tf, t.dim_c1, t.dim_g, g);
}

BdpVec bdp_ig_ff(const BoundaryTensors& t, const BdpVec& g) {
  return unary(t.ig_ff, t.dim_closed, t.dim_g, g);
}

std::array<double, 3> bdp_readout(const BoundaryTensors& t, const BdpVec& closed) {
  std::array<double, 3> out{0.0, 0.0, 0.0};
  for (int ch = 0; ch < t.n_channels && ch < 3; ++ch) {
    double s = 0.0;
    const int base = ch * t.dim_closed;
    for (int i = 0; i < t.dim_closed; ++i)
      s += t.readout[static_cast<std::size_t>(base + i)] * closed[static_cast<std::size_t>(i)];
    out[static_cast<std::size_t>(ch)] = s;
  }
  return out;
}

int full_bdp_dimension(int palette, FullBdpInterface interface) {
  if (palette == 2) {
    namespace n = boundary_dp_full_k2;
    switch (interface) {
      case FullBdpInterface::Closed: return n::dim_Closed;
      case FullBdpInterface::C1: return n::dim_C1;
      case FullBdpInterface::C2: return n::dim_C2;
      case FullBdpInterface::G: return n::dim_G;
    }
  }
  if (palette == 4) {
    namespace n = boundary_dp_full_k4;
    switch (interface) {
      case FullBdpInterface::Closed: return n::dim_Closed;
      case FullBdpInterface::C1: return n::dim_C1;
      case FullBdpInterface::C2: return n::dim_C2;
      case FullBdpInterface::G: return n::dim_G;
    }
  }
  return 0;
}

int full_bdp_max_dimension(int palette) {
  return std::max(
      {full_bdp_dimension(palette, FullBdpInterface::Closed),
       full_bdp_dimension(palette, FullBdpInterface::C1),
       full_bdp_dimension(palette, FullBdpInterface::C2),
       full_bdp_dimension(palette, FullBdpInterface::G)});
}

int full_bdp_channel_count(int palette) {
  if (palette == 2) return boundary_dp_full_k2::n_channels;
  if (palette == 4) return boundary_dp_full_k4::n_channels;
  return 0;
}

long long full_bdp_thread_scalar_fmas() {
  return g_full_bdp_scalar_fmas;
}

FullBdpVec full_bdp_leaf(int palette, int color, double weight) {
  if (color < 0 || color >= palette) {
    throw std::invalid_argument("LEAF_TO_C: color outside palette");
  }
  FullBdpVec output;
  if (palette == 2) {
    namespace n = boundary_dp_full_k2;
    output = sparse_full_constant(palette, FullBdpInterface::C1,
                                  FullBdpOrientation::U,
                                  n::LEAF_ABSENT_C1);
    const auto add = [&](const auto& terms) {
      g_full_bdp_scalar_fmas += static_cast<long long>(terms.size());
      for (const auto& term : terms) {
        output.values[term.index] += term.sign > 0 ? weight : -weight;
      }
    };
    if (color == 0) add(n::LEAF_PRESENT_0_C1);
    else add(n::LEAF_PRESENT_1_C1);
    return output;
  }
  if (palette == 4) {
    namespace n = boundary_dp_full_k4;
    output = sparse_full_constant(palette, FullBdpInterface::C1,
                                  FullBdpOrientation::U,
                                  n::LEAF_ABSENT_C1);
    const auto add = [&](const auto& terms) {
      g_full_bdp_scalar_fmas += static_cast<long long>(terms.size());
      for (const auto& term : terms) {
        output.values[term.index] += term.sign > 0 ? weight : -weight;
      }
    };
    switch (color) {
      case 0: add(n::LEAF_PRESENT_0_C1); break;
      case 1: add(n::LEAF_PRESENT_1_C1); break;
      case 2: add(n::LEAF_PRESENT_2_C1); break;
      case 3: add(n::LEAF_PRESENT_3_C1); break;
      default: break;
    }
    return output;
  }
  throw std::invalid_argument("unsupported boundary-DP palette");
}

FullBdpDual full_bdp_leaf_dual(int palette,
                               int current_color,
                               int background_color,
                               double weight) {
  return FullBdpDual{full_bdp_leaf(palette, current_color, weight),
                     full_bdp_leaf(palette, background_color, weight),
                     current_color == background_color || weight == 0.0};
}

FullBdpVec full_bdp_empty_center(int palette) {
  if (palette == 2) {
    return sparse_full_constant(palette, FullBdpInterface::G,
                                FullBdpOrientation::None,
                                boundary_dp_full_k2::CENTER_EMPTY_G);
  }
  if (palette == 4) {
    return sparse_full_constant(palette, FullBdpInterface::G,
                                FullBdpOrientation::None,
                                boundary_dp_full_k4::CENTER_EMPTY_G);
  }
  throw std::invalid_argument("unsupported boundary-DP palette");
}

FullBdpVec full_bdp_identity_c2(int palette) {
  if (palette == 2) {
    return sparse_full_constant(palette, FullBdpInterface::C2,
                                FullBdpOrientation::None,
                                boundary_dp_full_k2::IDENTITY_C2);
  }
  if (palette == 4) {
    return sparse_full_constant(palette, FullBdpInterface::C2,
                                FullBdpOrientation::None,
                                boundary_dp_full_k4::IDENTITY_C2);
  }
  throw std::invalid_argument("unsupported boundary-DP palette");
}

FullBdpVec full_bdp_empty_closed(int palette) {
  if (palette == 2) {
    return sparse_full_constant(palette, FullBdpInterface::Closed,
                                FullBdpOrientation::None,
                                boundary_dp_full_k2::CLOSED_EMPTY);
  }
  if (palette == 4) {
    return sparse_full_constant(palette, FullBdpInterface::Closed,
                                FullBdpOrientation::None,
                                boundary_dp_full_k4::CLOSED_EMPTY);
  }
  throw std::invalid_argument("unsupported boundary-DP palette");
}

FullBdpVec full_bdp_c_to_g(int palette, const FullBdpVec& c1_u) {
  require_full_value(palette, c1_u, FullBdpInterface::C1,
                     FullBdpOrientation::U, "C_TO_G");
  if (palette == 2) {
    return generated_full_unary<
        boundary_dp_full_k2::eval_C_TO_G_C1_G>(
        palette, c1_u, FullBdpInterface::G, FullBdpOrientation::None,
        boundary_dp_full_k2::C_TO_G_C1_G);
  }
  return generated_full_unary<boundary_dp_full_k4::eval_C_TO_G_C1_G>(
      palette, c1_u, FullBdpInterface::G, FullBdpOrientation::None,
      boundary_dp_full_k4::C_TO_G_C1_G);
}

FullBdpVec full_bdp_ig_ff(int palette, const FullBdpVec& g) {
  require_full_value(palette, g, FullBdpInterface::G,
                     FullBdpOrientation::None, "IG_FF");
  if (palette == 2) {
    return generated_full_unary<
        boundary_dp_full_k2::eval_IG_FF_G_CLOSED>(
        palette, g, FullBdpInterface::Closed, FullBdpOrientation::None,
        boundary_dp_full_k2::IG_FF_G_CLOSED);
  }
  return generated_full_unary<boundary_dp_full_k4::eval_IG_FF_G_CLOSED>(
      palette, g, FullBdpInterface::Closed, FullBdpOrientation::None,
      boundary_dp_full_k4::IG_FF_G_CLOSED);
}

FullBdpVec full_bdp_ig_tf(int palette, const FullBdpVec& g) {
  require_full_value(palette, g, FullBdpInterface::G,
                     FullBdpOrientation::None, "IG_TF");
  if (palette == 2) {
    return generated_full_unary<boundary_dp_full_k2::eval_IG_TF_G_C1>(
        palette, g, FullBdpInterface::C1, FullBdpOrientation::U,
        boundary_dp_full_k2::IG_TF_G_C1);
  }
  return generated_full_unary<boundary_dp_full_k4::eval_IG_TF_G_C1>(
      palette, g, FullBdpInterface::C1, FullBdpOrientation::U,
      boundary_dp_full_k4::IG_TF_G_C1);
}

FullBdpVec full_bdp_ig_tt(int palette, const FullBdpVec& g) {
  require_full_value(palette, g, FullBdpInterface::G,
                     FullBdpOrientation::None, "IG_TT");
  if (palette == 2) {
    return generated_full_unary<boundary_dp_full_k2::eval_IG_TT_G_C2>(
        palette, g, FullBdpInterface::C2, FullBdpOrientation::None,
        boundary_dp_full_k2::IG_TT_G_C2);
  }
  return generated_full_unary<boundary_dp_full_k4::eval_IG_TT_G_C2>(
      palette, g, FullBdpInterface::C2, FullBdpOrientation::None,
      boundary_dp_full_k4::IG_TT_G_C2);
}

FullBdpVec full_bdp_close(int palette, const FullBdpVec& c1_u) {
  require_full_value(palette, c1_u, FullBdpInterface::C1,
                     FullBdpOrientation::U, "CLOSE");
  if (palette == 2) {
    return generated_full_unary<
        boundary_dp_full_k2::eval_CLOSE_C1_CLOSED>(
        palette, c1_u, FullBdpInterface::Closed, FullBdpOrientation::None,
        boundary_dp_full_k2::CLOSE_C1_CLOSED);
  }
  return generated_full_unary<boundary_dp_full_k4::eval_CLOSE_C1_CLOSED>(
      palette, c1_u, FullBdpInterface::Closed, FullBdpOrientation::None,
      boundary_dp_full_k4::CLOSE_C1_CLOSED);
}

FullBdpVec full_bdp_reverse(int palette, const FullBdpVec& c) {
  if (c.interface == FullBdpInterface::C1) {
    if (c.orientation != FullBdpOrientation::U &&
        c.orientation != FullBdpOrientation::D) {
      throw std::invalid_argument("REVERSE C1: missing orientation");
    }
    require_full_value(palette, c, FullBdpInterface::C1, c.orientation,
                       "REVERSE C1");
    const FullBdpOrientation output_orientation =
        c.orientation == FullBdpOrientation::U ? FullBdpOrientation::D
                                                : FullBdpOrientation::U;
    if (palette == 2) {
      return generated_full_unary<boundary_dp_full_k2::eval_REV_C1_C1>(
          palette, c, FullBdpInterface::C1, output_orientation,
          boundary_dp_full_k2::REV_C1_C1);
    }
    return generated_full_unary<boundary_dp_full_k4::eval_REV_C1_C1>(
        palette, c, FullBdpInterface::C1, output_orientation,
        boundary_dp_full_k4::REV_C1_C1);
  }
  require_full_value(palette, c, FullBdpInterface::C2,
                     FullBdpOrientation::None, "REVERSE C2");
  if (palette == 2) {
    return generated_full_unary<boundary_dp_full_k2::eval_REV_C2_C2>(
        palette, c, FullBdpInterface::C2, FullBdpOrientation::None,
        boundary_dp_full_k2::REV_C2_C2);
  }
  return generated_full_unary<boundary_dp_full_k4::eval_REV_C2_C2>(
      palette, c, FullBdpInterface::C2, FullBdpOrientation::None,
      boundary_dp_full_k4::REV_C2_C2);
}

FullBdpVec full_bdp_gg(int palette,
                       const FullBdpVec& left,
                       const FullBdpVec& right) {
  require_full_value(palette, left, FullBdpInterface::G,
                     FullBdpOrientation::None, "GG_TO_G left");
  require_full_value(palette, right, FullBdpInterface::G,
                     FullBdpOrientation::None, "GG_TO_G right");
  if (palette == 2) {
    return generated_full_binary<boundary_dp_full_k2::eval_GG_TO_G_G_G_G>(
        palette, left, right, FullBdpInterface::G, FullBdpOrientation::None,
        boundary_dp_full_k2::GG_TO_G_G_G_G);
  }
  return generated_full_binary<boundary_dp_full_k4::eval_GG_TO_G_G_G_G>(
      palette, left, right, FullBdpInterface::G, FullBdpOrientation::None,
      boundary_dp_full_k4::GG_TO_G_G_G_G);
}

FullBdpVec full_bdp_cc(int palette,
                       const FullBdpVec& higher,
                       const FullBdpVec& lower) {
  const bool higher_c1 = higher.interface == FullBdpInterface::C1;
  const bool lower_c1 = lower.interface == FullBdpInterface::C1;
  if (higher_c1) {
    require_full_value(palette, higher, FullBdpInterface::C1,
                       FullBdpOrientation::D, "CC_TO_C higher");
  } else {
    require_full_value(palette, higher, FullBdpInterface::C2,
                       FullBdpOrientation::None, "CC_TO_C higher");
  }
  if (lower_c1) {
    require_full_value(palette, lower, FullBdpInterface::C1,
                       FullBdpOrientation::U, "CC_TO_C lower");
  } else {
    require_full_value(palette, lower, FullBdpInterface::C2,
                       FullBdpOrientation::None, "CC_TO_C lower");
  }

  const FullBdpInterface output_interface =
      higher_c1 && lower_c1 ? FullBdpInterface::Closed
                            : (higher_c1 || lower_c1
                                   ? FullBdpInterface::C1
                                   : FullBdpInterface::C2);
  const FullBdpOrientation output_orientation =
      output_interface != FullBdpInterface::C1
          ? FullBdpOrientation::None
          : (higher_c1 ? FullBdpOrientation::D : FullBdpOrientation::U);

  if (palette == 2) {
    namespace n = boundary_dp_full_k2;
    if (higher_c1 && lower_c1) {
      return generated_full_binary<n::eval_CC_TO_C_C1_C1_CLOSED>(
          palette, higher, lower, output_interface, output_orientation,
          n::CC_TO_C_C1_C1_CLOSED);
    }
    if (higher_c1) {
      return generated_full_binary<n::eval_CC_TO_C_C1_C2_C1>(
          palette, higher, lower, output_interface, output_orientation,
          n::CC_TO_C_C1_C2_C1);
    }
    if (lower_c1) {
      return generated_full_binary<n::eval_CC_TO_C_C2_C1_C1>(
          palette, higher, lower, output_interface, output_orientation,
          n::CC_TO_C_C2_C1_C1);
    }
    return generated_full_binary<n::eval_CC_TO_C_C2_C2_C2>(
        palette, higher, lower, output_interface, output_orientation,
        n::CC_TO_C_C2_C2_C2);
  }
  namespace n = boundary_dp_full_k4;
  if (higher_c1 && lower_c1) {
    return generated_full_binary<n::eval_CC_TO_C_C1_C1_CLOSED>(
        palette, higher, lower, output_interface, output_orientation,
        n::CC_TO_C_C1_C1_CLOSED);
  }
  if (higher_c1) {
    return generated_full_binary<n::eval_CC_TO_C_C1_C2_C1>(
        palette, higher, lower, output_interface, output_orientation,
        n::CC_TO_C_C1_C2_C1);
  }
  if (lower_c1) {
    return generated_full_binary<n::eval_CC_TO_C_C2_C1_C1>(
        palette, higher, lower, output_interface, output_orientation,
        n::CC_TO_C_C2_C1_C1);
  }
  return generated_full_binary<n::eval_CC_TO_C_C2_C2_C2>(
      palette, higher, lower, output_interface, output_orientation,
      n::CC_TO_C_C2_C2_C2);
}

std::array<double, 3> full_bdp_readout(int palette,
                                       const FullBdpVec& closed) {
  require_full_value(palette, closed, FullBdpInterface::Closed,
                     FullBdpOrientation::None, "READOUT");
  if (palette == 2) {
    return sparse_full_readout(closed,
                               boundary_dp_full_k2::READOUT_CLOSED);
  }
  return sparse_full_readout(closed,
                             boundary_dp_full_k4::READOUT_CLOSED);
}

void full_bdp_leaf_into(int palette,
                        int color,
                        double weight,
                        FullBdpMutableView output) {
  if (color < 0 || color >= palette) {
    throw std::invalid_argument("LEAF_TO_C: color outside palette");
  }
  require_full_hot_view(palette, output, FullBdpInterface::C1,
                        FullBdpOrientation::U, "LEAF_TO_C output");
  const auto add = [&](const auto& terms) {
    g_full_bdp_scalar_fmas += static_cast<long long>(terms.size());
    for (const auto& term : terms) {
      output.values[term.index] += term.sign > 0 ? weight : -weight;
    }
  };
  if (palette == 2) {
    namespace n = boundary_dp_full_k2;
    sparse_full_constant_into(output, n::LEAF_ABSENT_C1);
    if (color == 0) add(n::LEAF_PRESENT_0_C1);
    else add(n::LEAF_PRESENT_1_C1);
    return;
  }
  namespace n = boundary_dp_full_k4;
  sparse_full_constant_into(output, n::LEAF_ABSENT_C1);
  switch (color) {
    case 0: add(n::LEAF_PRESENT_0_C1); break;
    case 1: add(n::LEAF_PRESENT_1_C1); break;
    case 2: add(n::LEAF_PRESENT_2_C1); break;
    case 3: add(n::LEAF_PRESENT_3_C1); break;
    default: break;
  }
}

#define WQFM_FULL_BDP_UNARY_INTO(NAME, INPUT_INTERFACE, INPUT_ORIENTATION,   \
                                 OUTPUT_INTERFACE, OUTPUT_ORIENTATION, ID)  \
  void NAME(int palette, FullBdpConstView input,                            \
            FullBdpMutableView output) {                                    \
    require_full_hot_view(palette, input, INPUT_INTERFACE, INPUT_ORIENTATION,\
                          #NAME " input");                                  \
    require_full_hot_view(palette, output, OUTPUT_INTERFACE,                \
                          OUTPUT_ORIENTATION, #NAME " output");             \
    if (palette == 2) {                                                       \
      namespace n = boundary_dp_full_k2;                                    \
      generated_full_unary_into<n::eval_##ID>(input, output, n::ID);        \
      return;                                                                \
    }                                                                        \
    namespace n = boundary_dp_full_k4;                                      \
    generated_full_unary_into<n::eval_##ID>(input, output, n::ID);          \
  }

WQFM_FULL_BDP_UNARY_INTO(
    full_bdp_c_to_g_into, FullBdpInterface::C1, FullBdpOrientation::U,
    FullBdpInterface::G, FullBdpOrientation::None, C_TO_G_C1_G)
WQFM_FULL_BDP_UNARY_INTO(
    full_bdp_ig_ff_into, FullBdpInterface::G, FullBdpOrientation::None,
    FullBdpInterface::Closed, FullBdpOrientation::None, IG_FF_G_CLOSED)
WQFM_FULL_BDP_UNARY_INTO(
    full_bdp_ig_tf_into, FullBdpInterface::G, FullBdpOrientation::None,
    FullBdpInterface::C1, FullBdpOrientation::U, IG_TF_G_C1)
WQFM_FULL_BDP_UNARY_INTO(
    full_bdp_ig_tt_into, FullBdpInterface::G, FullBdpOrientation::None,
    FullBdpInterface::C2, FullBdpOrientation::None, IG_TT_G_C2)
WQFM_FULL_BDP_UNARY_INTO(
    full_bdp_close_into, FullBdpInterface::C1, FullBdpOrientation::U,
    FullBdpInterface::Closed, FullBdpOrientation::None, CLOSE_C1_CLOSED)

#undef WQFM_FULL_BDP_UNARY_INTO

#define WQFM_FULL_BDP_BACKGROUND_UNARY_INTO(                               \
    NAME, INPUT_INTERFACE, INPUT_ORIENTATION, OUTPUT_INTERFACE,             \
    OUTPUT_ORIENTATION, ID)                                                 \
  void NAME(int palette, FullBdpConstView input,                            \
            FullBdpMutableView output) {                                    \
    require_full_hot_view(palette, input, INPUT_INTERFACE, INPUT_ORIENTATION,\
                          #NAME " input");                                  \
    require_full_hot_view(palette, output, OUTPUT_INTERFACE,                \
                          OUTPUT_ORIENTATION, #NAME " output");             \
    if (palette == 2) {                                                      \
      namespace n = boundary_dp_full_k2;                                    \
      generated_full_unary_into<n::eval_background_##ID>(                  \
          input, output, n::BACKGROUND_##ID);                               \
      return;                                                               \
    }                                                                       \
    namespace n = boundary_dp_full_k4;                                      \
    generated_full_unary_into<n::eval_background_##ID>(                    \
        input, output, n::BACKGROUND_##ID);                                 \
  }

WQFM_FULL_BDP_BACKGROUND_UNARY_INTO(
    full_bdp_c_to_g_background_into, FullBdpInterface::C1,
    FullBdpOrientation::U, FullBdpInterface::G,
    FullBdpOrientation::None, C_TO_G_C1_G)
WQFM_FULL_BDP_BACKGROUND_UNARY_INTO(
    full_bdp_ig_ff_background_into, FullBdpInterface::G,
    FullBdpOrientation::None, FullBdpInterface::Closed,
    FullBdpOrientation::None, IG_FF_G_CLOSED)
WQFM_FULL_BDP_BACKGROUND_UNARY_INTO(
    full_bdp_ig_tf_background_into, FullBdpInterface::G,
    FullBdpOrientation::None, FullBdpInterface::C1,
    FullBdpOrientation::U, IG_TF_G_C1)
WQFM_FULL_BDP_BACKGROUND_UNARY_INTO(
    full_bdp_ig_tt_background_into, FullBdpInterface::G,
    FullBdpOrientation::None, FullBdpInterface::C2,
    FullBdpOrientation::None, IG_TT_G_C2)
WQFM_FULL_BDP_BACKGROUND_UNARY_INTO(
    full_bdp_close_background_into, FullBdpInterface::C1,
    FullBdpOrientation::U, FullBdpInterface::Closed,
    FullBdpOrientation::None, CLOSE_C1_CLOSED)

#undef WQFM_FULL_BDP_BACKGROUND_UNARY_INTO

void full_bdp_reverse_into(int palette,
                           FullBdpConstView input,
                           FullBdpMutableView output) {
  if (input.interface == FullBdpInterface::C1) {
    if ((input.orientation != FullBdpOrientation::U &&
         input.orientation != FullBdpOrientation::D) ||
        output.orientation == input.orientation) {
      throw std::invalid_argument("REVERSE C1 view orientation mismatch");
    }
    require_full_hot_view(palette, input, FullBdpInterface::C1,
                          input.orientation, "REVERSE C1 input");
    require_full_hot_view(palette, output, FullBdpInterface::C1,
                          output.orientation, "REVERSE C1 output");
    if (palette == 2) {
      namespace n = boundary_dp_full_k2;
      generated_full_unary_into<n::eval_REV_C1_C1>(input, output,
                                                    n::REV_C1_C1);
      return;
    }
    namespace n = boundary_dp_full_k4;
    generated_full_unary_into<n::eval_REV_C1_C1>(input, output,
                                                  n::REV_C1_C1);
    return;
  }
  require_full_hot_view(palette, input, FullBdpInterface::C2,
                        FullBdpOrientation::None, "REVERSE C2 input");
  require_full_hot_view(palette, output, FullBdpInterface::C2,
                        FullBdpOrientation::None, "REVERSE C2 output");
  if (palette == 2) {
    namespace n = boundary_dp_full_k2;
    generated_full_unary_into<n::eval_REV_C2_C2>(input, output,
                                                  n::REV_C2_C2);
    return;
  }
  namespace n = boundary_dp_full_k4;
  generated_full_unary_into<n::eval_REV_C2_C2>(input, output,
                                                n::REV_C2_C2);
}

void full_bdp_reverse_background_into(int palette,
                                      FullBdpConstView input,
                                      FullBdpMutableView output) {
  if (input.interface == FullBdpInterface::C1) {
    if ((input.orientation != FullBdpOrientation::U &&
         input.orientation != FullBdpOrientation::D) ||
        output.orientation == input.orientation) {
      throw std::invalid_argument("background REVERSE C1 orientation mismatch");
    }
    require_full_hot_view(palette, input, FullBdpInterface::C1,
                          input.orientation, "background REVERSE C1 input");
    require_full_hot_view(palette, output, FullBdpInterface::C1,
                          output.orientation, "background REVERSE C1 output");
    if (palette == 2) {
      namespace n = boundary_dp_full_k2;
      generated_full_unary_into<n::eval_background_REV_C1_C1>(
          input, output, n::BACKGROUND_REV_C1_C1);
      return;
    }
    namespace n = boundary_dp_full_k4;
    generated_full_unary_into<n::eval_background_REV_C1_C1>(
        input, output, n::BACKGROUND_REV_C1_C1);
    return;
  }
  require_full_hot_view(palette, input, FullBdpInterface::C2,
                        FullBdpOrientation::None,
                        "background REVERSE C2 input");
  require_full_hot_view(palette, output, FullBdpInterface::C2,
                        FullBdpOrientation::None,
                        "background REVERSE C2 output");
  if (palette == 2) {
    namespace n = boundary_dp_full_k2;
    generated_full_unary_into<n::eval_background_REV_C2_C2>(
        input, output, n::BACKGROUND_REV_C2_C2);
    return;
  }
  namespace n = boundary_dp_full_k4;
  generated_full_unary_into<n::eval_background_REV_C2_C2>(
      input, output, n::BACKGROUND_REV_C2_C2);
}

void full_bdp_gg_into(int palette,
                      FullBdpConstView left,
                      FullBdpConstView right,
                      FullBdpMutableView output) {
  require_full_hot_view(palette, left, FullBdpInterface::G,
                        FullBdpOrientation::None, "GG_TO_G left");
  require_full_hot_view(palette, right, FullBdpInterface::G,
                        FullBdpOrientation::None, "GG_TO_G right");
  require_full_hot_view(palette, output, FullBdpInterface::G,
                        FullBdpOrientation::None, "GG_TO_G output");
  if (palette == 2) {
    namespace n = boundary_dp_full_k2;
    generated_full_binary_into<n::eval_GG_TO_G_G_G_G>(
        left, right, output, n::GG_TO_G_G_G_G);
    return;
  }
  namespace n = boundary_dp_full_k4;
  generated_full_binary_into<n::eval_GG_TO_G_G_G_G>(
      left, right, output, n::GG_TO_G_G_G_G);
}

void full_bdp_gg_background_into(int palette,
                                 FullBdpConstView left,
                                 FullBdpConstView right,
                                 FullBdpMutableView output) {
  require_full_hot_view(palette, left, FullBdpInterface::G,
                        FullBdpOrientation::None, "background GG left");
  require_full_hot_view(palette, right, FullBdpInterface::G,
                        FullBdpOrientation::None, "background GG right");
  require_full_hot_view(palette, output, FullBdpInterface::G,
                        FullBdpOrientation::None, "background GG output");
  if (palette == 2) {
    namespace n = boundary_dp_full_k2;
    generated_full_binary_into<n::eval_background_GG_TO_G_G_G_G>(
        left, right, output, n::BACKGROUND_GG_TO_G_G_G_G);
    return;
  }
  namespace n = boundary_dp_full_k4;
  generated_full_binary_into<n::eval_background_GG_TO_G_G_G_G>(
      left, right, output, n::BACKGROUND_GG_TO_G_G_G_G);
}

void full_bdp_cc_into(int palette,
                      FullBdpConstView higher,
                      FullBdpConstView lower,
                      FullBdpMutableView output) {
  const bool higher_c1 = higher.interface == FullBdpInterface::C1;
  const bool lower_c1 = lower.interface == FullBdpInterface::C1;
  require_full_hot_view(palette, higher,
                        higher_c1 ? FullBdpInterface::C1
                                  : FullBdpInterface::C2,
                        higher_c1 ? FullBdpOrientation::D
                                  : FullBdpOrientation::None,
                        "CC_TO_C higher");
  require_full_hot_view(palette, lower,
                        lower_c1 ? FullBdpInterface::C1
                                 : FullBdpInterface::C2,
                        lower_c1 ? FullBdpOrientation::U
                                 : FullBdpOrientation::None,
                        "CC_TO_C lower");
  const FullBdpInterface output_interface =
      higher_c1 && lower_c1
          ? FullBdpInterface::Closed
          : (higher_c1 || lower_c1 ? FullBdpInterface::C1
                                   : FullBdpInterface::C2);
  const FullBdpOrientation output_orientation =
      output_interface != FullBdpInterface::C1
          ? FullBdpOrientation::None
          : (higher_c1 ? FullBdpOrientation::D : FullBdpOrientation::U);
  require_full_hot_view(palette, output, output_interface,
                        output_orientation, "CC_TO_C output");
#define WQFM_CC_INTO_DISPATCH(NS)                                            \
  if (higher_c1 && lower_c1) {                                               \
    generated_full_binary_into<NS::eval_CC_TO_C_C1_C1_CLOSED>(              \
        higher, lower, output, NS::CC_TO_C_C1_C1_CLOSED);                   \
    return;                                                                  \
  }                                                                          \
  if (higher_c1) {                                                           \
    generated_full_binary_into<NS::eval_CC_TO_C_C1_C2_C1>(                  \
        higher, lower, output, NS::CC_TO_C_C1_C2_C1);                       \
    return;                                                                  \
  }                                                                          \
  if (lower_c1) {                                                            \
    generated_full_binary_into<NS::eval_CC_TO_C_C2_C1_C1>(                  \
        higher, lower, output, NS::CC_TO_C_C2_C1_C1);                       \
    return;                                                                  \
  }                                                                          \
  generated_full_binary_into<NS::eval_CC_TO_C_C2_C2_C2>(                    \
      higher, lower, output, NS::CC_TO_C_C2_C2_C2);                          \
  return
  if (palette == 2) {
    namespace n = boundary_dp_full_k2;
    WQFM_CC_INTO_DISPATCH(n);
  }
  namespace n = boundary_dp_full_k4;
  WQFM_CC_INTO_DISPATCH(n);
#undef WQFM_CC_INTO_DISPATCH
}

void full_bdp_cc_background_into(int palette,
                                 FullBdpConstView higher,
                                 FullBdpConstView lower,
                                 FullBdpMutableView output) {
  const bool higher_c1 = higher.interface == FullBdpInterface::C1;
  const bool lower_c1 = lower.interface == FullBdpInterface::C1;
  require_full_hot_view(palette, higher,
                        higher_c1 ? FullBdpInterface::C1
                                  : FullBdpInterface::C2,
                        higher_c1 ? FullBdpOrientation::D
                                  : FullBdpOrientation::None,
                        "background CC higher");
  require_full_hot_view(palette, lower,
                        lower_c1 ? FullBdpInterface::C1
                                 : FullBdpInterface::C2,
                        lower_c1 ? FullBdpOrientation::U
                                 : FullBdpOrientation::None,
                        "background CC lower");
  const FullBdpInterface output_interface =
      higher_c1 && lower_c1
          ? FullBdpInterface::Closed
          : (higher_c1 || lower_c1 ? FullBdpInterface::C1
                                   : FullBdpInterface::C2);
  const FullBdpOrientation output_orientation =
      output_interface != FullBdpInterface::C1
          ? FullBdpOrientation::None
          : (higher_c1 ? FullBdpOrientation::D : FullBdpOrientation::U);
  require_full_hot_view(palette, output, output_interface,
                        output_orientation, "background CC output");
#define WQFM_BACKGROUND_CC_DISPATCH(NS)                                     \
  if (higher_c1 && lower_c1) {                                              \
    generated_full_binary_into<                                             \
        NS::eval_background_CC_TO_C_C1_C1_CLOSED>(                          \
        higher, lower, output, NS::BACKGROUND_CC_TO_C_C1_C1_CLOSED);       \
    return;                                                                 \
  }                                                                         \
  if (higher_c1) {                                                          \
    generated_full_binary_into<NS::eval_background_CC_TO_C_C1_C2_C1>(      \
        higher, lower, output, NS::BACKGROUND_CC_TO_C_C1_C2_C1);           \
    return;                                                                 \
  }                                                                         \
  if (lower_c1) {                                                           \
    generated_full_binary_into<NS::eval_background_CC_TO_C_C2_C1_C1>(      \
        higher, lower, output, NS::BACKGROUND_CC_TO_C_C2_C1_C1);           \
    return;                                                                 \
  }                                                                         \
  generated_full_binary_into<NS::eval_background_CC_TO_C_C2_C2_C2>(       \
      higher, lower, output, NS::BACKGROUND_CC_TO_C_C2_C2_C2);             \
  return
  if (palette == 2) {
    namespace n = boundary_dp_full_k2;
    WQFM_BACKGROUND_CC_DISPATCH(n);
  }
  namespace n = boundary_dp_full_k4;
  WQFM_BACKGROUND_CC_DISPATCH(n);
#undef WQFM_BACKGROUND_CC_DISPATCH
}

namespace {

template <class Masks>
std::uint8_t value_color_mask_impl(const FullBdpVec& value,
                                   const Masks& masks) {
  std::uint8_t mask = 0;
  for (int coordinate = 0; coordinate < value.dimension; ++coordinate) {
    if (value.values[static_cast<std::size_t>(coordinate)] != 0.0) {
      mask |= masks[static_cast<std::size_t>(coordinate)];
    }
  }
  return mask;
}

template <class Terms>
std::array<double, 3> color_mask_readout_impl(
    int palette,
    std::uint8_t color_mask,
    const double* closed,
    const Terms& terms) {
  std::array<double, 3> output{0.0, 0.0, 0.0};
  const std::uint64_t active = full_bdp_active_coordinate_bits(
      palette, FullBdpInterface::Closed, color_mask);
  for (const auto& term : terms) {
    if ((active & (std::uint64_t{1} << term.in)) == 0) {
      continue;
    }
    ++g_full_bdp_scalar_fmas;
    const int packed = std::popcount(
        active & ((std::uint64_t{1} << term.in) - 1u));
    output[term.out] +=
        term.sign > 0 ? closed[packed] : -closed[packed];
  }
  return output;
}

}  // namespace

std::uint64_t full_bdp_active_coordinate_bits(
    int palette,
    FullBdpInterface interface,
    std::uint8_t color_mask) {
  if (palette != 2 && palette != 4) {
    throw std::invalid_argument("unsupported boundary-DP palette");
  }
  if (color_mask >= (1u << palette)) {
    throw std::invalid_argument("boundary-DP color mask outside palette");
  }
#define WQFM_ACTIVE_BITS(IFACE)                                             \
  (palette == 2                                                             \
       ? boundary_dp_full_k2::active_coordinate_bits_##IFACE[color_mask]   \
       : boundary_dp_full_k4::active_coordinate_bits_##IFACE[color_mask])
  switch (interface) {
    case FullBdpInterface::Closed: return WQFM_ACTIVE_BITS(Closed);
    case FullBdpInterface::C1: return WQFM_ACTIVE_BITS(C1);
    case FullBdpInterface::C2: return WQFM_ACTIVE_BITS(C2);
    case FullBdpInterface::G: return WQFM_ACTIVE_BITS(G);
  }
#undef WQFM_ACTIVE_BITS
  return 0;
}

int full_bdp_active_coordinate_count(int palette,
                                     FullBdpInterface interface,
                                     std::uint8_t color_mask) {
  return std::popcount(
      full_bdp_active_coordinate_bits(palette, interface, color_mask));
}

int full_bdp_packed_coordinate_rank(int palette,
                                    FullBdpInterface interface,
                                    std::uint8_t color_mask,
                                    int coordinate) {
  if (!full_bdp_coordinate_active(
          palette, interface, coordinate, color_mask)) {
    return -1;
  }
  const std::uint64_t active =
      full_bdp_active_coordinate_bits(palette, interface, color_mask);
  return std::popcount(
      active & ((std::uint64_t{1} << coordinate) - 1u));
}

bool full_bdp_coordinate_active(int palette,
                                FullBdpInterface interface,
                                int coordinate,
                                std::uint8_t color_mask) {
  if (coordinate < 0 || coordinate >= full_bdp_dimension(palette, interface)) {
    throw std::out_of_range("boundary-DP coordinate outside interface");
  }
  return (full_bdp_active_coordinate_bits(palette, interface, color_mask) &
          (std::uint64_t{1} << coordinate)) != 0;
}

std::uint8_t full_bdp_value_color_mask(int palette,
                                       const FullBdpVec& value) {
#define WQFM_VALUE_MASK(IFACE)                                              \
  (palette == 2                                                             \
       ? value_color_mask_impl(                                             \
             value, boundary_dp_full_k2::coordinate_color_mask_##IFACE)     \
       : value_color_mask_impl(                                             \
             value, boundary_dp_full_k4::coordinate_color_mask_##IFACE))
  switch (value.interface) {
    case FullBdpInterface::Closed: return WQFM_VALUE_MASK(Closed);
    case FullBdpInterface::C1: return WQFM_VALUE_MASK(C1);
    case FullBdpInterface::C2: return WQFM_VALUE_MASK(C2);
    case FullBdpInterface::G: return WQFM_VALUE_MASK(G);
  }
#undef WQFM_VALUE_MASK
  return 0;
}

void full_bdp_leaf_color_mask_into(int palette,
                                   int color,
                                   double weight,
                                   std::uint8_t color_mask,
                                   double* output) {
  if (color < 0 || color >= palette || output == nullptr) {
    throw std::invalid_argument("invalid color-graded leaf");
  }
  std::uint64_t active = full_bdp_active_coordinate_bits(
      palette, FullBdpInterface::C1, color_mask);
  std::fill_n(output, std::popcount(active), 0.0);
  const auto add = [&](const auto& terms, double scale) {
    for (const auto& term : terms) {
      const std::uint64_t coordinate_bit =
          std::uint64_t{1} << term.index;
      if ((active & coordinate_bit) == 0) continue;
      ++g_full_bdp_scalar_fmas;
      const int packed = std::popcount(active & (coordinate_bit - 1u));
      output[packed] +=
          term.sign > 0 ? scale : -scale;
    }
  };
  if (palette == 2) {
    namespace n = boundary_dp_full_k2;
    add(n::LEAF_ABSENT_C1, 1.0);
    add(color == 0 ? n::LEAF_PRESENT_0_C1 : n::LEAF_PRESENT_1_C1,
        weight);
    return;
  }
  namespace n = boundary_dp_full_k4;
  add(n::LEAF_ABSENT_C1, 1.0);
  switch (color) {
    case 0: add(n::LEAF_PRESENT_0_C1, weight); break;
    case 1: add(n::LEAF_PRESENT_1_C1, weight); break;
    case 2: add(n::LEAF_PRESENT_2_C1, weight); break;
    case 3: add(n::LEAF_PRESENT_3_C1, weight); break;
    default: break;
  }
}

#define WQFM_COLOR_MASK_UNARY(NAME, ID)                                     \
  void NAME(int palette, std::uint8_t color_mask, const double* input,      \
            double* output) {                                               \
    int terms = 0;                                                          \
    if (palette == 2) {                                                     \
      terms = boundary_dp_full_k2::eval_color_mask_##ID(                    \
          color_mask, input, output);                                       \
    } else {                                                                \
      terms = boundary_dp_full_k4::eval_color_mask_##ID(                    \
          color_mask, input, output);                                       \
    }                                                                       \
    g_full_bdp_scalar_fmas += terms;                                        \
  }

WQFM_COLOR_MASK_UNARY(full_bdp_c_to_g_color_mask_into, C_TO_G_C1_G)
WQFM_COLOR_MASK_UNARY(full_bdp_ig_ff_color_mask_into, IG_FF_G_CLOSED)
WQFM_COLOR_MASK_UNARY(full_bdp_ig_tf_color_mask_into, IG_TF_G_C1)
WQFM_COLOR_MASK_UNARY(full_bdp_ig_tt_color_mask_into, IG_TT_G_C2)
WQFM_COLOR_MASK_UNARY(full_bdp_close_color_mask_into, CLOSE_C1_CLOSED)
WQFM_COLOR_MASK_UNARY(full_bdp_reverse_c1_color_mask_into, REV_C1_C1)
WQFM_COLOR_MASK_UNARY(full_bdp_reverse_c2_color_mask_into, REV_C2_C2)

#undef WQFM_COLOR_MASK_UNARY

void full_bdp_gg_color_mask_into(int palette,
                                 std::uint8_t left_mask,
                                 std::uint8_t right_mask,
                                 const double* left,
                                 const double* right,
                                 double* output) {
  const int terms = palette == 2
                        ? boundary_dp_full_k2::eval_color_mask_GG_TO_G_G_G_G(
                              left_mask, right_mask, left, right, output)
                        : boundary_dp_full_k4::eval_color_mask_GG_TO_G_G_G_G(
                              left_mask, right_mask, left, right, output);
  g_full_bdp_scalar_fmas += terms;
}

void full_bdp_cc_color_mask_into(int palette,
                                 FullBdpInterface higher_interface,
                                 FullBdpInterface lower_interface,
                                 std::uint8_t higher_mask,
                                 std::uint8_t lower_mask,
                                 const double* higher,
                                 const double* lower,
                                 double* output) {
  const bool higher_c1 = higher_interface == FullBdpInterface::C1;
  const bool lower_c1 = lower_interface == FullBdpInterface::C1;
  if ((!higher_c1 && higher_interface != FullBdpInterface::C2) ||
      (!lower_c1 && lower_interface != FullBdpInterface::C2)) {
    throw std::invalid_argument("invalid color-graded CC interface");
  }
#define WQFM_COLOR_MASK_CC_DISPATCH(NS)                                     \
  (higher_c1 && lower_c1                                                    \
       ? NS::eval_color_mask_CC_TO_C_C1_C1_CLOSED(                         \
             higher_mask, lower_mask, higher, lower, output)                \
       : (higher_c1                                                         \
              ? NS::eval_color_mask_CC_TO_C_C1_C2_C1(                      \
                    higher_mask, lower_mask, higher, lower, output)         \
              : (lower_c1                                                  \
                     ? NS::eval_color_mask_CC_TO_C_C2_C1_C1(               \
                           higher_mask, lower_mask, higher, lower, output)  \
                     : NS::eval_color_mask_CC_TO_C_C2_C2_C2(               \
                           higher_mask, lower_mask, higher, lower, output))))
  const int terms = palette == 2
                        ? WQFM_COLOR_MASK_CC_DISPATCH(boundary_dp_full_k2)
                        : WQFM_COLOR_MASK_CC_DISPATCH(boundary_dp_full_k4);
#undef WQFM_COLOR_MASK_CC_DISPATCH
  g_full_bdp_scalar_fmas += terms;
}

std::array<double, 3> full_bdp_color_mask_readout(
    int palette,
    std::uint8_t color_mask,
    const double* closed) {
  if (palette == 2) {
    return color_mask_readout_impl(
        palette, color_mask, closed,
        boundary_dp_full_k2::READOUT_CLOSED);
  }
  return color_mask_readout_impl(
      palette, color_mask, closed,
      boundary_dp_full_k4::READOUT_CLOSED);
}

std::array<double, 3> full_bdp_readout(int palette,
                                       FullBdpConstView closed) {
  require_full_hot_view(palette, closed, FullBdpInterface::Closed,
                        FullBdpOrientation::None, "READOUT");
  std::array<double, 3> output{0.0, 0.0, 0.0};
  const auto read = [&](const auto& terms) {
    g_full_bdp_scalar_fmas += static_cast<long long>(terms.size());
    for (const auto& term : terms) {
      output[term.out] +=
          term.sign > 0 ? closed.values[term.in] : -closed.values[term.in];
    }
  };
  if (palette == 2) read(boundary_dp_full_k2::READOUT_CLOSED);
  else read(boundary_dp_full_k4::READOUT_CLOSED);
  return output;
}

#define WQFM_FULL_BDP_DUAL_UNARY(NAME, INPUT_INTERFACE, INPUT_ORIENTATION,   \
                                 OUTPUT_INTERFACE, OUTPUT_ORIENTATION, ID)  \
  FullBdpDual NAME(int palette, const FullBdpDual& input) {                 \
    if (input.current_equals_background) {                                  \
      FullBdpVec current = NAME(palette, input.current);                     \
      return FullBdpDual{current, current, true};                            \
    }                                                                        \
    require_full_value(palette, input.current, INPUT_INTERFACE,              \
                       INPUT_ORIENTATION, #NAME " current");                \
    require_full_value(palette, input.background, INPUT_INTERFACE,           \
                       INPUT_ORIENTATION, #NAME " background");             \
    if (palette == 2) {                                                       \
      namespace n = boundary_dp_full_k2;                                    \
      return generated_full_unary_dual<n::eval_dual_##ID>(                  \
          palette, input, OUTPUT_INTERFACE, OUTPUT_ORIENTATION, n::ID);      \
    }                                                                        \
    namespace n = boundary_dp_full_k4;                                      \
    return generated_full_unary_dual<n::eval_dual_##ID>(                    \
        palette, input, OUTPUT_INTERFACE, OUTPUT_ORIENTATION, n::ID);        \
  }

WQFM_FULL_BDP_DUAL_UNARY(
    full_bdp_c_to_g, FullBdpInterface::C1, FullBdpOrientation::U,
    FullBdpInterface::G, FullBdpOrientation::None, C_TO_G_C1_G)
WQFM_FULL_BDP_DUAL_UNARY(
    full_bdp_ig_ff, FullBdpInterface::G, FullBdpOrientation::None,
    FullBdpInterface::Closed, FullBdpOrientation::None, IG_FF_G_CLOSED)
WQFM_FULL_BDP_DUAL_UNARY(
    full_bdp_ig_tf, FullBdpInterface::G, FullBdpOrientation::None,
    FullBdpInterface::C1, FullBdpOrientation::U, IG_TF_G_C1)
WQFM_FULL_BDP_DUAL_UNARY(
    full_bdp_ig_tt, FullBdpInterface::G, FullBdpOrientation::None,
    FullBdpInterface::C2, FullBdpOrientation::None, IG_TT_G_C2)
WQFM_FULL_BDP_DUAL_UNARY(
    full_bdp_close, FullBdpInterface::C1, FullBdpOrientation::U,
    FullBdpInterface::Closed, FullBdpOrientation::None, CLOSE_C1_CLOSED)

#undef WQFM_FULL_BDP_DUAL_UNARY

FullBdpDual full_bdp_reverse(int palette, const FullBdpDual& input) {
  if (input.current_equals_background) {
    FullBdpVec current = full_bdp_reverse(palette, input.current);
    return FullBdpDual{current, current, true};
  }
  if (input.current.interface == FullBdpInterface::C1) {
    if (input.current.orientation != input.background.orientation ||
        (input.current.orientation != FullBdpOrientation::U &&
         input.current.orientation != FullBdpOrientation::D)) {
      throw std::invalid_argument("REVERSE C1 dual orientation mismatch");
    }
    require_full_value(palette, input.current, FullBdpInterface::C1,
                       input.current.orientation, "REVERSE C1 current");
    require_full_value(palette, input.background, FullBdpInterface::C1,
                       input.current.orientation, "REVERSE C1 background");
    const FullBdpOrientation output_orientation =
        input.current.orientation == FullBdpOrientation::U
            ? FullBdpOrientation::D
            : FullBdpOrientation::U;
    if (palette == 2) {
      namespace n = boundary_dp_full_k2;
      return generated_full_unary_dual<n::eval_dual_REV_C1_C1>(
          palette, input, FullBdpInterface::C1, output_orientation,
          n::REV_C1_C1);
    }
    namespace n = boundary_dp_full_k4;
    return generated_full_unary_dual<n::eval_dual_REV_C1_C1>(
        palette, input, FullBdpInterface::C1, output_orientation,
        n::REV_C1_C1);
  }
  require_full_value(palette, input.current, FullBdpInterface::C2,
                     FullBdpOrientation::None, "REVERSE C2 current");
  require_full_value(palette, input.background, FullBdpInterface::C2,
                     FullBdpOrientation::None, "REVERSE C2 background");
  if (palette == 2) {
    namespace n = boundary_dp_full_k2;
    return generated_full_unary_dual<n::eval_dual_REV_C2_C2>(
        palette, input, FullBdpInterface::C2, FullBdpOrientation::None,
        n::REV_C2_C2);
  }
  namespace n = boundary_dp_full_k4;
  return generated_full_unary_dual<n::eval_dual_REV_C2_C2>(
      palette, input, FullBdpInterface::C2, FullBdpOrientation::None,
      n::REV_C2_C2);
}

FullBdpDual full_bdp_gg(int palette,
                        const FullBdpDual& left,
                        const FullBdpDual& right) {
  const bool identical = left.current_equals_background &&
                         right.current_equals_background;
  if (identical) {
    FullBdpVec current = full_bdp_gg(palette, left.current, right.current);
    return FullBdpDual{current, current, true};
  }
  for (const FullBdpVec* value : {&left.current, &left.background,
                                  &right.current, &right.background}) {
    require_full_value(palette, *value, FullBdpInterface::G,
                       FullBdpOrientation::None, "GG_TO_G dual");
  }
  if (palette == 2) {
    namespace n = boundary_dp_full_k2;
    return generated_full_binary_dual<n::eval_dual_GG_TO_G_G_G_G>(
        palette, left, right, FullBdpInterface::G, FullBdpOrientation::None,
        n::GG_TO_G_G_G_G);
  }
  namespace n = boundary_dp_full_k4;
  return generated_full_binary_dual<n::eval_dual_GG_TO_G_G_G_G>(
      palette, left, right, FullBdpInterface::G, FullBdpOrientation::None,
      n::GG_TO_G_G_G_G);
}

FullBdpDual full_bdp_cc(int palette,
                        const FullBdpDual& higher,
                        const FullBdpDual& lower) {
  const bool identical = higher.current_equals_background &&
                         lower.current_equals_background;
  if (identical) {
    FullBdpVec current =
        full_bdp_cc(palette, higher.current, lower.current);
    return FullBdpDual{current, current, true};
  }
  const bool higher_c1 = higher.current.interface == FullBdpInterface::C1;
  const bool lower_c1 = lower.current.interface == FullBdpInterface::C1;
  const FullBdpInterface higher_interface =
      higher_c1 ? FullBdpInterface::C1 : FullBdpInterface::C2;
  const FullBdpInterface lower_interface =
      lower_c1 ? FullBdpInterface::C1 : FullBdpInterface::C2;
  const FullBdpOrientation higher_orientation =
      higher_c1 ? FullBdpOrientation::D : FullBdpOrientation::None;
  const FullBdpOrientation lower_orientation =
      lower_c1 ? FullBdpOrientation::U : FullBdpOrientation::None;
  require_full_value(palette, higher.current, higher_interface,
                     higher_orientation, "CC_TO_C higher current");
  require_full_value(palette, higher.background, higher_interface,
                     higher_orientation, "CC_TO_C higher background");
  require_full_value(palette, lower.current, lower_interface,
                     lower_orientation, "CC_TO_C lower current");
  require_full_value(palette, lower.background, lower_interface,
                     lower_orientation, "CC_TO_C lower background");
  const FullBdpInterface output_interface =
      higher_c1 && lower_c1
          ? FullBdpInterface::Closed
          : (higher_c1 || lower_c1 ? FullBdpInterface::C1
                                   : FullBdpInterface::C2);
  const FullBdpOrientation output_orientation =
      output_interface != FullBdpInterface::C1
          ? FullBdpOrientation::None
          : (higher_c1 ? FullBdpOrientation::D : FullBdpOrientation::U);
#define WQFM_DUAL_CC_DISPATCH(NS)                                            \
  if (higher_c1 && lower_c1) {                                               \
    return generated_full_binary_dual<                                      \
        NS::eval_dual_CC_TO_C_C1_C1_CLOSED>(                                \
        palette, higher, lower, output_interface, output_orientation,        \
        NS::CC_TO_C_C1_C1_CLOSED);                                          \
  }                                                                          \
  if (higher_c1) {                                                           \
    return generated_full_binary_dual<NS::eval_dual_CC_TO_C_C1_C2_C1>(      \
        palette, higher, lower, output_interface, output_orientation,        \
        NS::CC_TO_C_C1_C2_C1);                                              \
  }                                                                          \
  if (lower_c1) {                                                            \
    return generated_full_binary_dual<NS::eval_dual_CC_TO_C_C2_C1_C1>(      \
        palette, higher, lower, output_interface, output_orientation,        \
        NS::CC_TO_C_C2_C1_C1);                                              \
  }                                                                          \
  return generated_full_binary_dual<NS::eval_dual_CC_TO_C_C2_C2_C2>(       \
      palette, higher, lower, output_interface, output_orientation,          \
      NS::CC_TO_C_C2_C2_C2)
  if (palette == 2) {
    namespace n = boundary_dp_full_k2;
    WQFM_DUAL_CC_DISPATCH(n);
  }
  namespace n = boundary_dp_full_k4;
  WQFM_DUAL_CC_DISPATCH(n);
#undef WQFM_DUAL_CC_DISPATCH
}

}  // namespace wqfm
