// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "PointwiseFunctions/AnalyticData/SelfForce/Scalar/CircularOrbit.hpp"

#include <complex>
#include <cstddef>
#include <effsource.hpp>
#include <utility>

#include "DataStructures/ComplexDataVector.hpp"
#include "DataStructures/DataBox/Prefixes.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/EagerMath/Magnitude.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Elliptic/Systems/SelfForce/Scalar/Tags.hpp"
#include "PointwiseFunctions/GeneralRelativity/TortoiseCoordinates.hpp"
#include "Utilities/Gsl.hpp"

namespace ScalarSelfForce::AnalyticData {

CircularOrbit::CircularOrbit(const double scalar_charge,
                             const double black_hole_mass,
                             const double black_hole_spin,
                             const double orbital_radius,
                             const int m_mode_number)
    : scalar_charge_(scalar_charge),
      black_hole_mass_(black_hole_mass),
      black_hole_spin_(black_hole_spin),
      orbital_radius_(orbital_radius),
      m_mode_number_(m_mode_number) {}

// Background
tuples::TaggedTuple<Tags::Alpha, Tags::Beta, Tags::Gamma>
CircularOrbit::variables(
    const tnsr::I<DataVector, 2>& x,
    tmpl::list<Tags::Alpha, Tags::Beta, Tags::Gamma> /*meta*/) const {
  const double a = black_hole_spin_ * black_hole_mass_;
  const double M = black_hole_mass_;
  const double r_plus = M * (1. + sqrt(1. - square(black_hole_spin_)));
  const double r_0 = orbital_radius_;
  const double omega = 1. / (a + sqrt(cube(r_0) * M));
  const auto& r_star = get<0>(x);
  const auto& theta = get<1>(x);
  const DataVector r = gr::boyer_lindquist_radius_minus_r_plus_from_tortoise(
                     r_star, M, black_hole_spin_) +
                 r_plus;
  const DataVector delta = square(r) - 2.0 * M * r + square(a);
  const DataVector r_sq_plus_a_sq = square(r) + square(a);
  const DataVector r_sq_plus_a_sq_sq = square(r_sq_plus_a_sq);
  const DataVector sin_theta = sin(theta);
  const DataVector cos_theta = cos(theta);
  const DataVector sin_theta_squared = square(sin_theta);
  const DataVector sigma_squared =
      r_sq_plus_a_sq_sq - square(a) * delta * sin_theta_squared;
  tuples::TaggedTuple<Tags::Alpha, Tags::Beta, Tags::Gamma> result{};
  auto& alpha = get<Tags::Alpha>(result);
  auto& beta = get<Tags::Beta>(result);
  auto& gamma = get<Tags::Gamma>(result);
  get(alpha) = delta / r_sq_plus_a_sq_sq;
  const ComplexDataVector temp1 =
      1. / r * std::complex<double>(0., 2. * a * m_mode_number_);
  get(beta) = (-square(m_mode_number_ * omega) * sigma_squared +
               4. * a * square(m_mode_number_) * omega * M * r +
               delta * (square(m_mode_number_) / sin_theta_squared +
                        2. * M / r * (1. - square(a) / M / r) + temp1)) /
              r_sq_plus_a_sq_sq;
  get<0>(gamma) =
      -1. / r_sq_plus_a_sq * std::complex<double>(0., 2. * a * m_mode_number_) +
      2. * square(a) * get(alpha) / r;
  get<1>(gamma) = -get(alpha) * cos_theta / sin_theta;
  return result;
}

// Initial guess
tuples::TaggedTuple<Tags::MMode> CircularOrbit::variables(
    const tnsr::I<DataVector, 2>& x, tmpl::list<Tags::MMode> /*meta*/) const {
  tuples::TaggedTuple<Tags::MMode> result{};
  auto& field = get<Tags::MMode>(result);
  get<0>(field) = DataVector{get<0>(x).size(), 0.};
  get<1>(field) = DataVector{get<0>(x).size(), 0.};
  return result;
}

// Fixed sources
tuples::TaggedTuple<::Tags::FixedSource<Tags::MMode>> CircularOrbit::variables(
    const tnsr::I<DataVector, 2>& x,
    tmpl::list<::Tags::FixedSource<Tags::MMode>> /*meta*/) const {
  const double a = black_hole_spin_ * black_hole_mass_;
  const double M = black_hole_mass_;
  const double r_0 = orbital_radius_;
  const double r_plus = M * (1. + sqrt(1. - square(black_hole_spin_)));
  const double r_minus = M * (1. - sqrt(1. - square(black_hole_spin_)));
  {
    // Initialize effsource
    effsource_init(M, a);
    struct coordinate xp;
    xp.t = 0;
    xp.r = r_0;
    xp.theta = M_PI_2;
    xp.phi = 0;
    const double e = ((r_0 - 2.0 * M) * sqrt(M * r_0) + a * M) /
                     (sqrt(M * r_0) * sqrt(r_0 * r_0 - 3.0 * M * r_0 +
                                           2.0 * a * sqrt(M * r_0)));
    const double l = (M * (a * a + r_0 * r_0 - 2.0 * a * sqrt(M * r_0))) /
                     (sqrt(M * r_0) * sqrt(r_0 * r_0 - 3.0 * M * r_0 +
                                           2.0 * a * sqrt(M * r_0)));
    effsource_set_particle(&xp, e, l, 0.);
  }
  const auto& r_star = get<0>(x);
  const auto& theta = get<1>(x);
  const DataVector r = gr::boyer_lindquist_radius_minus_r_plus_from_tortoise(
                           r_star, M, black_hole_spin_) +
                       r_plus;
  // TODO: check sign of delta_phi
  const DataVector delta_phi = m_mode_number_ * a / (r_plus - r_minus) *
                               log((r - r_plus) / (r - r_minus));
  const ComplexDataVector rotation =
      cos(delta_phi) + std::complex<double>(0., 1.) * sin(delta_phi);
  tuples::TaggedTuple<::Tags::FixedSource<Tags::MMode>> result{};
  ComplexDataVector effective_source{get<0>(x).size()};
  struct coordinate x_i;
  double PhiS[2], dPhiS_dx[8], d2PhiS_dx2[20], src[2];
  for (size_t i = 0; i < get<0>(x).size(); ++i) {
    x_i.t = 0;
    x_i.r = r[i];
    x_i.theta = theta[i];
    x_i.phi = 0;
    effsource_calc_m(m_mode_number_, &x_i, PhiS, dPhiS_dx, d2PhiS_dx2, src);
    effective_source[i] = src[0] + std::complex<double>(0., 1.) * src[1];
  }
  // Rotate the source by delta_phi and multiply by r / 2 pi
  // TODO: is this correct for S_eff?
  effective_source *= rotation * 0.5 * r / M_PI;
  auto& fixed_source = get<::Tags::FixedSource<Tags::MMode>>(result);
  get<0>(fixed_source) = real(effective_source);
  get<1>(fixed_source) = imag(effective_source);
  return result;
}

void CircularOrbit::pup(PUP::er& p) {
  elliptic::analytic_data::Background::pup(p);
  elliptic::analytic_data::InitialGuess::pup(p);
  p | scalar_charge_;
  p | black_hole_mass_;
  p | black_hole_spin_;
  p | orbital_radius_;
  p | m_mode_number_;
}

bool operator==(const CircularOrbit& lhs, const CircularOrbit& rhs) {
  return lhs.scalar_charge_ == rhs.scalar_charge_ and
         lhs.black_hole_mass_ == rhs.black_hole_mass_ and
         lhs.black_hole_spin_ == rhs.black_hole_spin_ and
         lhs.orbital_radius_ == rhs.orbital_radius_ and
         lhs.m_mode_number_ == rhs.m_mode_number_;
}

bool operator!=(const CircularOrbit& lhs, const CircularOrbit& rhs) {
  return not(lhs == rhs);
}

PUP::able::PUP_ID CircularOrbit::my_PUP_ID = 0;  // NOLINT

}  // namespace ScalarSelfForce::AnalyticData
