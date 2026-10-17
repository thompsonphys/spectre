// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Elliptic/Systems/SelfForce/Scalar/AnalyticData/EccentricOrbit.hpp"

#include <algorithm>
#include <blaze/math/Vector.h>
#include <boost/iterator/is_iterator.hpp>
#include <complex.h>
#include <limits>
#include <memory>
extern "C" {
#include <effsource_equatorial.h>
}
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_sf.h>
#include <iomanip>
#include <ios>
#include <stdexcept>
#include <stdio.h>
#include <stdlib.h>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <string>
#include "DataStructures/Tensor/IndexType.hpp"
#include "Elliptic/Systems/SelfForce/Scalar/AnalyticData/SelfForceBackground.hpp"
extern "C" {
#include "korb.h"
}

#include <cmath>
#include <complex>
#include <cstddef>
#include <effsource_equatorial.hpp>
#include <gsl/gsl_errno.h>
#include <utility>

#include "DataStructures/Blaze/IntegerPow.hpp"
#include "DataStructures/ComplexDataVector.hpp"
#include "DataStructures/DataBox/Prefixes.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/EagerMath/Magnitude.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Elliptic/Systems/SelfForce/Scalar/Tags.hpp"
#include "NumericalAlgorithms/Integration/GslQuadAdaptive.hpp"
#include "NumericalAlgorithms/Interpolation/CubicSpline.hpp"
#include "PointwiseFunctions/GeneralRelativity/TortoiseCoordinates.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/Math.hpp"
#include "Utilities/Serialization/PupStlCpp17.hpp"
#include "Utilities/System/ParallelInfo.hpp"

namespace {
// Instantaneous time Jacobian dt/dlambda at Mino time lambda.
double dtdlambda(double lambda, const korb_params& o) {
  return korb_Tr(korb_psifromla(lambda, o), o) +
         korb_Tth(korb_chifromla(lambda, o), o);
}

// Newton iteration for lambda such that t(lambda) = t.
// Sets stalled = true if it fails to converge in 20 iterations.
double lambda_from_t(double t, double lambda, korb_params& orbpar,
                     bool& stalled) {
  for (int it = 0; it < 20; it++) {
    double dl = (korb_tfromla(lambda, orbpar) - t) / dtdlambda(lambda, orbpar);
    lambda -= dl;
    if (std::abs(dl) <=
        4.0 * std::numeric_limits<double>::epsilon() * fabs(lambda)) {
      return lambda;
    }
  }
  stalled = true;
  return lambda;
}
}  // namespace

namespace ScalarSelfForce::AnalyticData {

namespace {
std::pair<DataVector, DataVector> boost_function_and_deriv(
    const DataVector& r_star, const std::array<double, 4>& transition_points) {
  return {
      smoothstep<1>(transition_points[0], transition_points[1], r_star) +
          smoothstep<1>(transition_points[2], transition_points[3], r_star) -
          1.0,
      smoothstep_deriv<1>(transition_points[0], transition_points[1], r_star) +
          smoothstep_deriv<1>(transition_points[2], transition_points[3],
                              r_star)};
}
}  // namespace

namespace EccentricOrbit_detail {
size_t phi_points_for_sample(const double c1, const int m_mode_number) {
  constexpr size_t min_phi_points = 64;
  constexpr size_t max_phi_points = size_t{1} << 22;
  constexpr double convergence_constant = 16.;
  constexpr double min_mode_bound = 20.;
  if (not(c1 > 0.)) {
    ERROR("The separation parameter C1 = "
          << c1
          << " is not positive. The field point is at the particle position. "
             "The m-mode of the effective source is not finite there.");
  }
  const double abs_m = std::abs(m_mode_number);
  const double wanted_phi_points = std::max(
      {static_cast<double>(min_phi_points),
       std::max(abs_m, min_mode_bound) + convergence_constant / std::sqrt(c1),
       2.5 * abs_m});
  if (not(wanted_phi_points <= static_cast<double>(max_phi_points))) {
    ERROR("The numerical m-mode projection needs "
          << wanted_phi_points << " azimuthal samples at the separation C1 = "
          << c1 << ". The limit is " << max_phi_points
          << " samples. The field point is too close to the particle.");
  }
  size_t phi_points = min_phi_points;
  while (static_cast<double>(phi_points) < wanted_phi_points) {
    phi_points *= 2;
  }
  return phi_points;
}
}  // namespace EccentricOrbit_detail

EccentricOrbit::EccentricOrbit(
    const double black_hole_mass, const double black_hole_spin,
    const double semi_latus_rectum, const double eccentricity,
    const int m_mode_number, const int n_mode_number,
    const std::optional<std::array<double, 4>>
        hyperboloidal_slicing_transitions,
    const bool penetrating_horizon, const bool impose_equatorial_symmetry,
    const size_t time_points, const std::optional<size_t> phi_points)
    : black_hole_mass_(black_hole_mass),
      black_hole_spin_(black_hole_spin),
      semi_latus_rectum_(semi_latus_rectum),
      eccentricity_(eccentricity),
      m_mode_number_(m_mode_number),
      n_mode_number_(n_mode_number),
      hyperboloidal_slicing_transitions_(hyperboloidal_slicing_transitions),
      penetrating_horizon_(penetrating_horizon),
      impose_equatorial_symmetry_(impose_equatorial_symmetry),
      time_points_(time_points),
      phi_points_(phi_points) {
  if (penetrating_horizon_ and
      not hyperboloidal_slicing_transitions_.has_value()) {
    ERROR(
        "Hyperboloidal slicing must be enabled when penetrating_horizon is "
        "true.");
  }
  if (phi_points_.has_value() and
      *phi_points_ <= 2 * static_cast<size_t>(std::abs(m_mode_number_))) {
    ERROR("PhiPoints (" << *phi_points_
                        << ") must be larger than 2 * |MModeNumber| ("
                        << 2 * std::abs(m_mode_number_)
                        << "). Fewer samples alias the m-mode.");
  }
  compute_trajectory(time_points_);
}

EccentricOrbit::EccentricOrbit(CkMigrateMessage* m)
    : SelfForceBackground(m), elliptic::analytic_data::InitialGuess(m) {}

tnsr::I<double, 2> EccentricOrbit::puncture_position() const {
  const double M = black_hole_mass_;
  const double r_plus = M * (1. + sqrt(1. - square(black_hole_spin_)));
  const double r_0 = semi_latus_rectum_;
  if (penetrating_horizon_) {
    return tnsr::I<double, 2>{{{r_0, 0.}}};
  } else {
  const double r_star = gr::tortoise_radius_from_boyer_lindquist_minus_r_plus(
      r_0 - r_plus, M, black_hole_spin_);
  return tnsr::I<double, 2>{{{r_star, 0.}}};
  }
}

// Compute the particle trajectory and specify the number of time points used
// over one radial period
void EccentricOrbit::compute_trajectory(size_t time_points) {
  const int ecc = 1;
  const int inclined = 0;
  const double err = 1.0e-12;
  const double x_inclination = 1.0;
  korb_params orbpar;
  korb_getparams(ecc, inclined, black_hole_spin_, semi_latus_rectum_,
                 eccentricity_, x_inclination, err, &orbpar);
  if(eccentricity_ == 0)
  {
    orbpar.wr = orbpar.wphi; // For eccentricity = 0 case
  }
  energy = orbpar.E;
  angular_momentum = orbpar.Lz;
  Omega_r = orbpar.wr;
  Omega_phi = orbpar.wphi;
  r_of_t.resize(time_points);
  phi_of_t.resize(time_points);
  t_values.resize(time_points);
  u_r.resize(time_points);

  // Sample one full radial period. `variables` only evaluates the effective
  // source in the first half and gets the second half from time-reversal
  // symmetry.
  double T_max = 2*M_PI / (orbpar.wr);

  // Compute the radial four velocity at each time
  double X_four_velocity = orbpar.Lz - (black_hole_spin_ * orbpar.E);
  double lambda = 0;
  for (size_t i = 0; i < time_points_; i++) {
    double t = static_cast<double>(i) / time_points_ * T_max;
    t_values[i] = t;
    bool stalled=false;
    lambda = lambda_from_t(
        t, lambda + (T_max/time_points_)/(dtdlambda(lambda, orbpar)),
        orbpar, stalled);
    phi_of_t[i] = korb_phifromla(lambda, orbpar);
    double psi_of_t = korb_psifromla(lambda, orbpar);
    r_of_t[i] = korb_rfrompsi(psi_of_t, orbpar);
    u_r[i] = (eccentricity_ * gsl_sf_sin(psi_of_t) / semi_latus_rectum_) *
             sqrt(square(X_four_velocity) + square(black_hole_spin_) +
                  (2 * X_four_velocity * black_hole_spin_ * orbpar.E) -
                  (2 * square(X_four_velocity) *
                   (3 + eccentricity_ * gsl_sf_cos(psi_of_t)) /
                   semi_latus_rectum_));
  }

  korb_freepar(orbpar);

  create_source_coefficients();
}

void EccentricOrbit::create_source_coefficients() {
  for (auto* ctx : source_coefficients) {
    effsource_equatorial_free(ctx);
  }
  source_coefficients.clear();
  // Only the first half of the radial period, see `variables`
  source_coefficients.resize(std::min(time_points_, time_points_ / 2 + 1));
  for (size_t i = 0; i < source_coefficients.size(); ++i) {
    coordinate xp{};
    xp.t = t_values[i];
    xp.r = r_of_t[i];
    xp.theta = M_PI_2;
    // The particle is placed at phi = 0 and the m-modes are rotated to
    // phi_of_t[i] in `variables`
    xp.phi = 0.;
    source_coefficients[i] =
        effsource_equatorial_create(black_hole_mass_, black_hole_spin_);
    effsource_equatorial_ctx_set_particle(source_coefficients[i], &xp, energy,
                                          angular_momentum, u_r[i]);
  }
}

// Background
tuples::TaggedTuple<Tags::Alpha, Tags::Beta, Tags::Gamma>
EccentricOrbit::variables(
  const tnsr::I<DataVector, 2>& x,
  tmpl::list<Tags::Alpha, Tags::Beta, Tags::Gamma> /*meta*/) const {
  const double a = black_hole_spin_ * black_hole_mass_;
  const double M = black_hole_mass_;
  const double r_plus = M * (1. + sqrt(1. - square(black_hole_spin_)));
  const double r_minus = M * (1. - sqrt(1. - square(black_hole_spin_)));
  const double k = m_mode_number() * Omega_phi + n_mode_number() * Omega_r;

  // Resolve coordinates
  const auto& r_star_or_r = get<0>(x);
  DataVector r;
  DataVector r_star;
  DataVector r_minus_r_plus;
  if (penetrating_horizon_) {
    // NOLINTNEXTLINE
    r.set_data_ref(const_cast<DataVector*>(&r_star_or_r));
    r_minus_r_plus = r - r_plus;
    r_star = gr::tortoise_radius_from_boyer_lindquist_minus_r_plus(
        r_minus_r_plus, M, black_hole_spin_);
  } else {
    // NOLINTNEXTLINE
    r_star.set_data_ref(const_cast<DataVector*>(&r_star_or_r));
    r_minus_r_plus = gr::boyer_lindquist_radius_minus_r_plus_from_tortoise(
        r_star, M, black_hole_spin_);
    r = r_minus_r_plus + r_plus;
  }
  const auto& cos_theta_or_sq = get<1>(x);
  DataVector cos_theta_sq;
  if (impose_equatorial_symmetry_) {
    // NOLINTNEXTLINE
    cos_theta_sq.set_data_ref(const_cast<DataVector*>(&cos_theta_or_sq));
  } else {
    cos_theta_sq = square(cos_theta_or_sq);
  }
  const DataVector delta = r_minus_r_plus * (r - r_minus);
  const DataVector r_sq_plus_a_sq = square(r) + square(a);
  const DataVector r_sq_plus_a_sq_sq = square(r_sq_plus_a_sq);
  const DataVector sin_theta_squared = 1. - cos_theta_sq;
  const DataVector sigma_squared =
      r_sq_plus_a_sq_sq - square(a) * delta * sin_theta_squared;
  tuples::TaggedTuple<Tags::Alpha, Tags::Beta, Tags::Gamma> result{};
  // Hyperboloidal slicing
  ComplexDataVector H;
  ComplexDataVector dH;
  if (hyperboloidal_slicing_transitions_.has_value()) {
    std::tie(H, dH) = boost_function_and_deriv(
        r_star_or_r, hyperboloidal_slicing_transitions_.value());
  } else {
    H = make_with_value<ComplexDataVector>(r_star_or_r, 0.);
    dH = make_with_value<ComplexDataVector>(r_star_or_r, 0.);
  }
  auto& alpha = get<Tags::Alpha>(result);
  auto& beta = get<Tags::Beta>(result);
  auto& gamma = get<Tags::Gamma>(result);
  if (penetrating_horizon_) {
    get<0>(alpha) = delta / r_sq_plus_a_sq;
    get<1>(alpha) = 1.0 / r_sq_plus_a_sq;
  } else {
    get<0>(alpha) = make_with_value<DataVector>(r_star_or_r, 1.0);
    get<1>(alpha) = delta / r_sq_plus_a_sq_sq;
  }
  get(beta) = make_with_value<ComplexDataVector>(r_star_or_r, 0.);
  for (size_t p = 0; p < get(beta).size(); ++p) {
    if (penetrating_horizon_ and equal_within_roundoff(r[p], r_plus)) {
      // The following terms are zero at the horizon. Skip them to avoid
      // division by zero.
      continue;
    }
    get(beta)[p] =
        square(k) *
            (square(H[p]) - sigma_squared[p] / r_sq_plus_a_sq_sq[p]) +
        2. * a * m_mode_number_ * k *
            (2. * M * r[p] / r_sq_plus_a_sq[p] + H[p]) / r_sq_plus_a_sq[p];
    if (penetrating_horizon_) {
      get(beta)[p] /= get<0>(alpha)[p];
    }
  }
  get(beta) +=
      get<1>(alpha) * (m_mode_number_ * (m_mode_number_ + 1) +
                       2. * M / r * (1. - square(a) / M / r) +
                       std::complex<double>(0., 2. * a * m_mode_number_) *
                           (1. + a * Omega_phi * H) / r) -
      std::complex<double>(0., k) * dH;
  get<0>(gamma) =
      2. * square(a) * delta / (r * r_sq_plus_a_sq_sq) -
      std::complex<double>(0., 2. * a * m_mode_number_) / r_sq_plus_a_sq -
      std::complex<double>(0., 2. * k) * H;
  get<1>(gamma) = 2. * m_mode_number_ * cos_theta_or_sq * get<1>(alpha);
  if (impose_equatorial_symmetry_) {
    get<1>(gamma) += sin_theta_squared * get<1>(alpha);
    get<1>(gamma) *= 2.0;
  }
  get<1>(alpha) *= sin_theta_squared;
  if (impose_equatorial_symmetry_) {
    get<1>(alpha) *= 4. * cos_theta_sq;
  }
  return result;
}

// Initial guess
tuples::TaggedTuple<Tags::MMode> EccentricOrbit::variables(
    const tnsr::I<DataVector, 2>& x, tmpl::list<Tags::MMode> /*meta*/) {
  tuples::TaggedTuple<Tags::MMode> result{};
  auto& field = get<Tags::MMode>(result);
  get(field) = ComplexDataVector{get<0>(x).size(), 0.};
  return result;
}

// Fixed sources
tuples::TaggedTuple<
    ::Tags::FixedSource<Tags::MMode>, Tags::SingularField,
    ::Tags::deriv<Tags::SingularField, tmpl::size_t<2>, Frame::Inertial>,
    Tags::BoyerLindquistRadius>
EccentricOrbit::variables(
    const tnsr::I<DataVector, 2>& x, bool on_worldtube_boundary,
    tmpl::list<
        ::Tags::FixedSource<Tags::MMode>, Tags::SingularField,
        ::Tags::deriv<Tags::SingularField, tmpl::size_t<2>, Frame::Inertial>,
        Tags::BoyerLindquistRadius> /*meta*/) const {
  const double a = black_hole_spin_ * black_hole_mass_;
  const double M = black_hole_mass_;
  const double r_plus = M * (1. + sqrt(1. - square(black_hole_spin_)));
  const double r_minus = M * (1. - sqrt(1. - square(black_hole_spin_)));
  const auto& r_star_or_r = get<0>(x);
  if (hyperboloidal_slicing_transitions_.has_value() and
      ((min(r_star_or_r) < (*hyperboloidal_slicing_transitions_)[1] and
        not equal_within_roundoff(min(r_star_or_r),
                                  (*hyperboloidal_slicing_transitions_)[1])) or
       (max(r_star_or_r) > (*hyperboloidal_slicing_transitions_)[2] and
        not equal_within_roundoff(max(r_star_or_r),
                                  (*hyperboloidal_slicing_transitions_)[2])))) {
    ERROR(
        "The effective source is only valid where no hyperboloidal slicing is "
        "applied, which is in the radial range ["
        << (*hyperboloidal_slicing_transitions_)[1] << ", "
        << (*hyperboloidal_slicing_transitions_)[2]
        << "], but was requested in the range [" << min(r_star_or_r) << ", "
        << max(r_star_or_r) << "]");
  }
  DataVector r;
  DataVector r_star;
  DataVector r_minus_r_plus;
  if (penetrating_horizon_) {
    // NOLINTNEXTLINE
    r.set_data_ref(const_cast<DataVector*>(&r_star_or_r));
    r_minus_r_plus = r - r_plus;
    r_star = gr::tortoise_radius_from_boyer_lindquist_minus_r_plus(
        r_minus_r_plus, M, black_hole_spin_);
  } else {
    // NOLINTNEXTLINE
    r_star.set_data_ref(const_cast<DataVector*>(&r_star_or_r));
    r_minus_r_plus = gr::boyer_lindquist_radius_minus_r_plus_from_tortoise(
        r_star, M, black_hole_spin_);
    r = r_minus_r_plus + r_plus;
  }
  const auto& cos_theta_or_sq = get<1>(x);
  DataVector cos_theta;
  DataVector cos_theta_sq;
  if (impose_equatorial_symmetry_) {
    // NOLINTNEXTLINE
    cos_theta_sq.set_data_ref(const_cast<DataVector*>(&cos_theta_or_sq));
    cos_theta = sqrt(cos_theta_or_sq);
  } else {
    // NOLINTNEXTLINE
    cos_theta.set_data_ref(const_cast<DataVector*>(&cos_theta_or_sq));
    cos_theta_sq = square(cos_theta_or_sq);
  }
  const DataVector delta = r_minus_r_plus * (r - r_minus);
  const DataVector r_sq_plus_a_sq = square(r) + square(a);
  const DataVector r_sq_plus_a_sq_sq = square(r_sq_plus_a_sq);
  const DataVector sin_theta_sq = 1. - cos_theta_sq;
  const DataVector sin_theta = sqrt(sin_theta_sq);
  const DataVector sin_theta_pow_m = integer_pow(sin_theta, m_mode_number_);
  const DataVector delta_phi = m_mode_number_ * a / (r_plus - r_minus) *
                               log((r - r_plus) / (r - r_minus));
  const ComplexDataVector rotation =
      cos(delta_phi) - std::complex<double>(0., 1.) * sin(delta_phi);
  tuples::TaggedTuple<
      ::Tags::FixedSource<Tags::MMode>, Tags::SingularField,
      ::Tags::deriv<Tags::SingularField, tmpl::size_t<2>, Frame::Inertial>,
      Tags::BoyerLindquistRadius>
      result{};
  get(get<Tags::BoyerLindquistRadius>(result)) = r;
  const size_t num_points = get<0>(x).size();
  Scalar<ComplexDataVector>& effective_source =
      get<::Tags::FixedSource<Tags::MMode>>(result);
  get(effective_source).destructive_resize(num_points);
  Scalar<ComplexDataVector>& singular_field = get<Tags::SingularField>(result);
  get(singular_field).destructive_resize(num_points);
  tnsr::i<ComplexDataVector, 2>& deriv_singular_field =
      get<::Tags::deriv<Tags::SingularField, tmpl::size_t<2>, Frame::Inertial>>(
          result);
  get<0>(deriv_singular_field).destructive_resize(num_points);
  get<1>(deriv_singular_field).destructive_resize(num_points);

  // Project the effective source (and the singular field on the worldtube
  // boundary) onto the m-mode at each time point and onto the n-mode over one
  // radial period. Both projections are trapezoidal sums. Each snapshot is
  // evaluated with the particle at phi_p = 0, so its m-mode is rotated by
  // exp(-i m phi_p(t_i)) in the n-mode sum.
  //
  // The phases exp(-i m phi_k) of the finest azimuthal grid used so far. The
  // automatic grids are powers of two, so a coarser grid reads every
  // stride-th entry. A fixed grid always has the size of the table.
  std::vector<std::complex<double>> m_mode_phases{};
  const auto update_m_mode_phases = [this,
                                     &m_mode_phases](const size_t num_phi) {
    if (num_phi <= m_mode_phases.size()) {
      return;
    }
    m_mode_phases.resize(num_phi);
    for (size_t k = 0; k < num_phi; ++k) {
      m_mode_phases[k] =
          std::polar(1., -m_mode_number_ * 2. * M_PI * static_cast<double>(k) /
                             static_cast<double>(num_phi));
    }
  };
  const auto n_mode_weight = [this](const size_t i) {
    return std::polar(1. / static_cast<double>(time_points_),
                      m_mode_number_ * (Omega_phi * t_values[i] - phi_of_t[i]) +
                          2. * M_PI * n_mode_number_ * static_cast<double>(i) /
                              static_cast<double>(time_points_));
  };

  ComplexDataVector effective_source_n_mode{num_points, 0.};
  ComplexDataVector singular_field_n_mode{num_points, 0.};
  ComplexDataVector deriv_r_singular_field_n_mode{num_points, 0.};
  ComplexDataVector deriv_theta_singular_field_n_mode{num_points, 0.};
  const DataVector theta = acos(cos_theta);
  coordinate x_i{};
  double singular_field_4d = 0.;
  std::array<double, 4> deriv_singular_field_4d{};
  std::array<double, 10> second_deriv_singular_field_4d{};
  double effective_source_4d = 0.;
  // The orbit is symmetric under (t, phi) -> (-t, -phi): t, phi and the radial
  // phase are odd functions of Mino time, so r is even and u^r is odd. The
  // snapshot at t_{N-i} = T_r - t_i is the snapshot at t_i mirrored in phi
  // (both with phi_p = 0), so its m-mode is the complex conjugate. Only the
  // snapshots in the first half of the radial period are evaluated.
  for (size_t i = 0; i < source_coefficients.size(); ++i) {
    const size_t mirror_i = time_points_ - i;
    const bool has_mirror = i > 0 and mirror_i != i;
    const std::complex<double> weight = n_mode_weight(i);
    const std::complex<double> mirror_weight =
        has_mirror ? n_mode_weight(mirror_i) : std::complex<double>{};
    x_i.t = t_values[i];
    const effsource_equatorial_ctx& ctx = *source_coefficients[i];
    for (size_t j = 0; j < num_points; ++j) {
      x_i.r = r[j];
      x_i.theta = theta[j];
      const size_t num_phi =
          phi_points_.has_value()
              ? *phi_points_
              : EccentricOrbit_detail::phi_points_for_sample(
                    (ctx.alpha20 * square(r[j] - ctx.xp.r) +
                     ctx.alpha02 * square(theta[j] - ctx.xp.theta)) /
                        ctx.beta,
                    m_mode_number_);
      update_m_mode_phases(num_phi);
      const size_t stride = m_mode_phases.size() / num_phi;
      const double phi_spacing = 2. * M_PI / static_cast<double>(num_phi);
      std::complex<double> effective_source_m_mode{};
      std::complex<double> singular_field_m_mode{};
      std::complex<double> deriv_r_singular_field_m_mode{};
      std::complex<double> deriv_theta_singular_field_m_mode{};
      for (size_t k = 0; k < num_phi; ++k) {
        x_i.phi = phi_spacing * static_cast<double>(k);
        effsource_equatorial_ctx_calc(
            source_coefficients[i], &x_i, &singular_field_4d,
            deriv_singular_field_4d.data(),
            second_deriv_singular_field_4d.data(), &effective_source_4d);
        const std::complex<double>& phase = m_mode_phases[k * stride];
        effective_source_m_mode += effective_source_4d * phase;
        if (on_worldtube_boundary) {
          singular_field_m_mode += singular_field_4d * phase;
          deriv_r_singular_field_m_mode += deriv_singular_field_4d[1] * phase;
          deriv_theta_singular_field_m_mode +=
              deriv_singular_field_4d[2] * phase;
        }
      }
      effective_source_m_mode *= phi_spacing;
      singular_field_m_mode *= phi_spacing;
      deriv_r_singular_field_m_mode *= phi_spacing;
      deriv_theta_singular_field_m_mode *= phi_spacing;
      effective_source_n_mode[j] +=
          weight * effective_source_m_mode +
          mirror_weight * std::conj(effective_source_m_mode);
      if (on_worldtube_boundary) {
        singular_field_n_mode[j] +=
            weight * singular_field_m_mode +
            mirror_weight * std::conj(singular_field_m_mode);
        deriv_r_singular_field_n_mode[j] +=
            weight * deriv_r_singular_field_m_mode +
            mirror_weight * std::conj(deriv_r_singular_field_m_mode);
        deriv_theta_singular_field_n_mode[j] +=
            weight * deriv_theta_singular_field_m_mode +
            mirror_weight * std::conj(deriv_theta_singular_field_m_mode);
      }
    }
  }

  // Rotate by delta_phi and multiply by r / 2 pi
  const ComplexDataVector m_mode_factor = rotation * 0.5 * r / M_PI;
  // Factor Delta * (r^2 + a^2 cos^2(theta)) / Sigma^2
  // Factor Sigma^2 / (r^2 + a^2)^2 from first-order formulation
  // Factor 1/sin(theta)^m from change of variables
  get(effective_source) = effective_source_n_mode * m_mode_factor *
                          (square(r) + square(a) * cos_theta_sq) /
                          (r_sq_plus_a_sq * sin_theta_pow_m);
  if (not penetrating_horizon_) {
    get(effective_source) *= delta / r_sq_plus_a_sq;
  }

  // Only compute the singular field and its derivatives when crossing the
  // worldtube boundary
  if (not on_worldtube_boundary) {
    get(singular_field) = 0.;
    get<0>(deriv_singular_field) = 0.;
    get<1>(deriv_singular_field) = 0.;
    return result;
  }
  get(singular_field) = singular_field_n_mode * m_mode_factor / sin_theta_pow_m;
  get<0>(deriv_singular_field) =
      deriv_r_singular_field_n_mode * m_mode_factor / sin_theta_pow_m +
      get(singular_field) / r -
      std::complex<double>(0., a * m_mode_number_) / delta *
          get(singular_field);
  if (not penetrating_horizon_) {
    get<0>(deriv_singular_field) *= delta / r_sq_plus_a_sq;
  }
  get<1>(deriv_singular_field) = deriv_theta_singular_field_n_mode *
                                 m_mode_factor / sin_theta_pow_m / (-sin_theta);
  if (impose_equatorial_symmetry_) {
    get<1>(deriv_singular_field) /= 2. * cos_theta;
  }
  {
    ComplexDataVector add_term =
        m_mode_number_ * get(singular_field) / sin_theta_sq;
    if (impose_equatorial_symmetry_) {
      add_term *= 0.5;
    } else {
      add_term *= cos_theta;
    }
    get<1>(deriv_singular_field) += add_term;
  }
  return result;
}

void EccentricOrbit::pup(PUP::er& p) {
  SelfForceBackground::pup(p);
  elliptic::analytic_data::InitialGuess::pup(p);
  p | black_hole_mass_;
  p | black_hole_spin_;
  p | semi_latus_rectum_;
  p | eccentricity_;
  p | m_mode_number_;
  p | n_mode_number_;
  p | hyperboloidal_slicing_transitions_;
  p | penetrating_horizon_;
  p | impose_equatorial_symmetry_;
  p | time_points_;
  p | phi_points_;
  p | t_values;
  p | r_of_t;
  p | phi_of_t;
  p | u_r;
  p | Omega_phi;
  p | Omega_r;
  p | energy;
  p | angular_momentum;

  if (p.isUnpacking()) {
    create_source_coefficients();
  }
}

bool operator==(const EccentricOrbit& lhs, const EccentricOrbit& rhs) {
  return lhs.black_hole_mass_ == rhs.black_hole_mass_ and
         lhs.black_hole_spin_ == rhs.black_hole_spin_ and
         lhs.semi_latus_rectum_ == rhs.semi_latus_rectum_ and
         lhs.eccentricity_ == rhs.eccentricity_ and
         lhs.m_mode_number_ == rhs.m_mode_number_ and
         lhs.n_mode_number_ == rhs.n_mode_number_ and
         lhs.hyperboloidal_slicing_transitions_ ==
             rhs.hyperboloidal_slicing_transitions_ and
         lhs.penetrating_horizon_ == rhs.penetrating_horizon_ and
         lhs.impose_equatorial_symmetry_ == rhs.impose_equatorial_symmetry_ and
         lhs.time_points_ == rhs.time_points_ and
         lhs.phi_points_ == rhs.phi_points_;
         }

bool operator!=(const EccentricOrbit& lhs, const EccentricOrbit& rhs) {
  return not(lhs == rhs);
}

PUP::able::PUP_ID EccentricOrbit::my_PUP_ID = 0;  // NOLINT

}  // namespace ScalarSelfForce::AnalyticData
