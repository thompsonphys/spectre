// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <effsource_equatorial.hpp>
#include <utility>

#include "DataStructures/ComplexDataVector.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/TaggedTuple.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "DataStructures/Variables.hpp"
#include "Domain/BlockLogicalCoordinates.hpp"
#include "Domain/Creators/Rectilinear.hpp"
#include "Domain/ElementMap.hpp"
#include "Domain/Structure/ElementId.hpp"
#include "Elliptic/Systems/SelfForce/Scalar/AnalyticData/EccentricOrbit.hpp"
#include "Elliptic/Systems/SelfForce/Scalar/Equations.hpp"
#include "NumericalAlgorithms/LinearOperators/Divergence.tpp"
#include "NumericalAlgorithms/LinearOperators/PartialDerivatives.hpp"
#include "NumericalAlgorithms/Spectral/LogicalCoordinates.hpp"
#include "PointwiseFunctions/GeneralRelativity/TortoiseCoordinates.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/TMPL.hpp"

namespace ScalarSelfForce::AnalyticData {

namespace {
// The m-modes of the effective source and of the singular field and its
// derivatives by a trapezoidal sum over phi, the integrals of their absolute
// values over phi, which bound the rounding error of the sum, and the m-mode of
// the effective source of `effsource_equatorial_ctx_calc`. Near the particle
// the effective source of `effsource_equatorial_ctx_calc` loses digits to
// cancellation (about 1e-10 relative at C1 ~ 1e-4), so the reference for the
// effective source is `effsource_equatorial_ctx_point_src` with the
// double-double coefficients of `effsource_equatorial_ctx_prepare_dd`.
struct NumericalMModes {
  std::array<std::complex<double>, 4> modes{};
  std::array<double, 4> abs_integrals{};
  std::complex<double> calc_effective_source{};
};

NumericalMModes numerical_m_modes(effsource_equatorial_ctx* ctx, const double r,
                                  const double theta, const int m_mode_number,
                                  const size_t num_phi) {
  effsource_equatorial_point point{};
  effsource_equatorial_ctx_prepare_dd(ctx, r, theta, &point);
  coordinate x{};
  x.t = ctx->xp.t;
  x.r = r;
  x.theta = theta;
  double singular_field = 0.;
  std::array<double, 4> deriv_singular_field{};
  std::array<double, 10> second_deriv_singular_field{};
  double calc_effective_source = 0.;
  NumericalMModes result{};
  const double phi_spacing = 2. * M_PI / static_cast<double>(num_phi);
  for (size_t k = 0; k < num_phi; ++k) {
    x.phi = phi_spacing * static_cast<double>(k);
    effsource_equatorial_ctx_calc(
        ctx, &x, &singular_field, deriv_singular_field.data(),
        second_deriv_singular_field.data(), &calc_effective_source);
    const std::complex<double> phase =
        std::polar(phi_spacing, -m_mode_number * x.phi);
    const std::array<double, 4> values{
        effsource_equatorial_ctx_point_src(ctx, &point, x.phi), singular_field,
        deriv_singular_field[1], deriv_singular_field[2]};
    for (size_t field = 0; field < 4; ++field) {
      gsl::at(result.modes, field) += gsl::at(values, field) * phase;
      gsl::at(result.abs_integrals, field) +=
          std::abs(gsl::at(values, field)) * phi_spacing;
    }
    result.calc_effective_source += calc_effective_source * phase;
  }
  return result;
}

void test_analytic_m_modes() {
  // Snapshots at periapsis (t = 0) and on the way out (t = T_r / 8) of an orbit
  // with r_p in [7.1, 16.7]
  const EccentricOrbit orbit{
      1., 0.5, 10., 0.4, 2, 2, {{-25., -5., 20., 40.}}, false, false, 8};
  // Offsets (dr, dtheta) of the field point from the particle, from
  // C1 ~ 1e-4 to C1 ~ 1. The offset (0.6, -0.03) has C1 ~ 2e-3 and 1e-3, just
  // above EFFSOURCE_PREPARE_DD_C1. The trapezoidal sum with 4096 samples is
  // converged to roundoff at all of them.
  const std::array<std::array<double, 2>, 6> offsets{{{{0.3, -0.01}},
                                                      {{-0.05, -0.02}},
                                                      {{0., -0.05}},
                                                      {{0.6, -0.03}},
                                                      {{2., -0.3}},
                                                      {{-1.5, -0.8}}}};
  for (const size_t i : {size_t{0}, size_t{1}}) {
    effsource_equatorial_ctx* ctx = orbit.source_coefficients.at(i);
    for (const auto& offset : offsets) {
      const double r = ctx->xp.r + offset[0];
      const double theta = ctx->xp.theta + offset[1];
      for (const int m_mode_number : {0, 2, -2, 7}) {
        CAPTURE(i);
        CAPTURE(r);
        CAPTURE(theta);
        CAPTURE(m_mode_number);
        const auto analytic =
            EccentricOrbit_detail::m_modes(*ctx, r, theta, m_mode_number, true);
        const auto numerical =
            numerical_m_modes(ctx, r, theta, m_mode_number, 4096);
        const std::array<std::complex<double>, 4> analytic_modes{
            analytic.effective_source, analytic.singular_field,
            analytic.deriv_r_singular_field,
            analytic.deriv_theta_singular_field};
        // For C1 >= EFFSOURCE_PREPARE_DD_C1 = 1e-3,
        // effsource_equatorial_ctx_prepare uses double precision. The
        // effective source then has errors of up to about 4e-11 relative to
        // the integral of its absolute value near C1 = 1e-3.
        const double dr = r - ctx->xp.r;
        const double dtheta = theta - ctx->xp.theta;
        const double c1 =
            (ctx->alpha20 * dr * dr + ctx->alpha02 * dtheta * dtheta) /
            ctx->beta;
        const std::array<double, 4> tolerances{
            c1 < EFFSOURCE_PREPARE_DD_C1 ? 1.e-12 : 1.e-10, 1.e-13, 1.e-13,
            1.e-13};
        for (size_t field = 0; field < 4; ++field) {
          CAPTURE(field);
          CAPTURE(gsl::at(analytic_modes, field));
          CAPTURE(gsl::at(numerical.modes, field));
          CHECK(std::abs(gsl::at(analytic_modes, field) -
                         gsl::at(numerical.modes, field)) <
                gsl::at(tolerances, field) *
                    gsl::at(numerical.abs_integrals, field));
        }
        CAPTURE(numerical.calc_effective_source);
        CHECK(std::abs(analytic.effective_source -
                       numerical.calc_effective_source) <
              1.e-8 * numerical.abs_integrals[0]);
        const auto source_only = EccentricOrbit_detail::m_modes(
            *ctx, r, theta, m_mode_number, false);
        CHECK(source_only.effective_source == analytic.effective_source);
        CHECK(source_only.singular_field == std::complex<double>{});
        CHECK(source_only.deriv_r_singular_field == std::complex<double>{});
        CHECK(source_only.deriv_theta_singular_field == std::complex<double>{});
      }
    }
  }
#ifdef SPECTRE_DEBUG
  CHECK_THROWS_WITH(
      EccentricOrbit_detail::m_modes(*orbit.source_coefficients.at(0),
                                     orbit.source_coefficients.at(0)->xp.r,
                                     M_PI_2, 2, false),
      Catch::Matchers::ContainsSubstring("is on the particle"));
  CHECK_THROWS_WITH(
      EccentricOrbit_detail::m_modes(*orbit.source_coefficients.at(0), 12., 1.,
                                     250, false),
      Catch::Matchers::ContainsSubstring("only available for |m| <= 249"));
#endif
}

// Checks the self-force equations and the effective source together: the
// elliptic operator applied to the singular field must give the effective
// source. This is done numerically on a rectangular grid in (r_*, cos(theta)),
// or (r_*, cos(theta)^2) with equatorial symmetry, near the puncture.
void test_operator_consistency(const bool impose_equatorial_symmetry) {
  CAPTURE(impose_equatorial_symmetry);
  // Set up a domain. With equatorial symmetry the angular coordinate is
  // cos(theta)^2 instead of cos(theta).
  const double costheta_lower = 0.1;
  const double costheta_upper = 0.3;
  const double angular_lower = impose_equatorial_symmetry
                                   ? costheta_lower * costheta_lower
                                   : costheta_lower;
  const double angular_upper = impose_equatorial_symmetry
                                   ? costheta_upper * costheta_upper
                                   : costheta_upper;
  const double rstar_offset = 0.;
  const double delta_rstar = 5.;
  const size_t npoints = 20;
  const domain::creators::Rectangle domain_creator{
      {{rstar_offset, angular_lower}},
      {{rstar_offset + delta_rstar, angular_upper}},
      {{1, 1}},
      {{npoints, npoints}},
      {{false, false}}};
  const auto domain = domain_creator.create_domain();
  const auto& block = domain.blocks()[0];
  const ElementId<2> element_id{0};
  const ElementMap<2, Frame::Inertial> element_map{element_id, block};
  const Mesh<2> mesh{npoints, Spectral::Basis::Legendre,
                     Spectral::Quadrature::GaussLobatto};
  const auto xi = logical_coordinates(mesh);
  const auto x = element_map(xi);
  const auto inv_jacobian = element_map.inv_jacobian(xi);
  const auto& r_star = get<0>(x);
  const auto& angular_coordinate = get<1>(x);
  CAPTURE(min(r_star));
  CAPTURE(max(r_star));
  CAPTURE(min(angular_coordinate));
  CAPTURE(max(angular_coordinate));

  // Get the analytic fields
  for (int n_mode_number = 2; n_mode_number < 3; ++n_mode_number) {

    CAPTURE(n_mode_number);
    const auto eccentric_orbit = EccentricOrbit{1.,
                                                0.5,
                                                10.,
                                                0.4,
                                                2,
                                                n_mode_number,
                                                {{-25., -5., 20., 40.}},
                                                false,
                                                impose_equatorial_symmetry,
                                                128};
    const auto background =
        eccentric_orbit.variables(x, EccentricOrbit::background_tags{});
    const auto& alpha = get<Tags::Alpha>(background);
    const auto& beta = get<Tags::Beta>(background);
    const auto& gamma = get<Tags::Gamma>(background);
    const auto vars = eccentric_orbit.variables(
        x, true, EccentricOrbit::source_tags{});
    const auto& singular_field = get<Tags::SingularField>(vars);
    const auto& deriv_singular_field = get<
        ::Tags::deriv<Tags::SingularField, tmpl::size_t<2>, Frame::Inertial>>(
        vars);
    const auto& effective_source = get<::Tags::FixedSource<Tags::MMode>>(vars);
    const auto numeric_deriv_singular_field =
        partial_derivative(singular_field, mesh, inv_jacobian);
    const Approx custom_approx = Approx::custom().epsilon(1.e-10).scale(1.);
    for (size_t i = 0; i < deriv_singular_field.size(); ++i) {
      CAPTURE(i);
      CHECK_ITERABLE_CUSTOM_APPROX(numeric_deriv_singular_field[i],
                                   deriv_singular_field[i], custom_approx);
    }

    tnsr::I<ComplexDataVector, 2> flux_singular_field{};
    ScalarSelfForce::Fluxes::apply(make_not_null(&flux_singular_field), alpha,
                                   {}, deriv_singular_field);
    auto scalar_eqn = divergence(flux_singular_field, mesh, inv_jacobian);
    get(scalar_eqn) *= -1.;
    ScalarSelfForce::Sources::apply(make_not_null(&scalar_eqn), beta, gamma,
                                    singular_field, deriv_singular_field,
                                    flux_singular_field);

    // Minus sign is from the definition of the effective source:
    //   \psi = \psi_R + \psi_P = 0
    // where \psi_R is the regular part and \psi_P is the singular part
    //   => -\Delta \psi_R = \Delta \psi_P = S_eff
    // where -Delta represents the elliptic operator. So the effective source
    // for the regular part is is the negative of the elliptic operator
    // acting on the singular part.
    CHECK_ITERABLE_CUSTOM_APPROX(get(scalar_eqn), -get(effective_source),
                                 custom_approx);
  }
}

// Outside the hyperboloidal window the slicing t' = t + h(r) with
// dh/dr = -H (r^2 + a^2) / Delta changes the mode function by exp(i psi),
// psi = omega h, with omega = m Omega_phi + n Omega_r. The operator there must
// therefore act on exp(i psi) f as the window operator acts on f. Comparing
// two orbits that differ only in where the window ends checks every term of
// the operator that depends on H, in the inner (H = -1) or the outer (H = +1)
// null-slicing region.
void test_slicing_covariance(const int n_mode_number, const bool inner) {
  CAPTURE(n_mode_number);
  CAPTURE(inner);
  const double r_lower = inner ? 2.3 : 30.;
  const double r_upper = inner ? 2.9 : 35.;
  const size_t npoints = 20;
  const domain::creators::Rectangle domain_creator{{{r_lower, 0.1}},
                                                   {{r_upper, 0.3}},
                                                   {{0, 0}},
                                                   {{npoints, npoints}},
                                                   {{false, false}}};
  const auto domain = domain_creator.create_domain();
  const ElementId<2> element_id{0};
  const ElementMap<2, Frame::Inertial> element_map{element_id,
                                                   domain.blocks()[0]};
  const Mesh<2> mesh{npoints, Spectral::Basis::Legendre,
                     Spectral::Quadrature::GaussLobatto};
  const auto xi = logical_coordinates(mesh);
  const auto x = element_map(xi);
  const auto inv_jacobian = element_map.inv_jacobian(xi);
  const DataVector& r = get<0>(x);
  const DataVector& z = get<1>(x);

  // The region lies in the null slicing of `null_slicing` and in the window of
  // `window`. The background needs only the orbital frequencies, so a few time
  // points suffice.
  const EccentricOrbit null_slicing{
      1.,   0.5,   10., 0.2, 2, n_mode_number, {{3., 3., 25., 25.}},
      true, false, 4};
  const EccentricOrbit window{1.,
                              0.5,
                              10.,
                              0.2,
                              2,
                              n_mode_number,
                              inner ? std::array<double, 4>{{2., 2., 25., 25.}}
                                    : std::array<double, 4>{{3., 3., 40., 40.}},
                              true,
                              false,
                              4};
  const double H = inner ? -1. : 1.;
  const double omega =
      2. * window.omega_phi() + n_mode_number * window.omega_r();
  const double M = 1.;
  const double a = 0.5;
  const double r_plus = M * (1. + sqrt(1. - square(a)));
  const DataVector r_minus_r_plus = r - r_plus;
  const DataVector r_star =
      gr::tortoise_radius_from_boyer_lindquist_minus_r_plus(r_minus_r_plus, M,
                                                            a);
  const DataVector psi = -H * omega * r_star;
  const ComplexDataVector phase =
      cos(psi) + std::complex<double>(0., 1.) * sin(psi);

  const auto apply_operator = [&mesh, &inv_jacobian, &x](
                                  const EccentricOrbit& orbit,
                                  const Scalar<ComplexDataVector>& field) {
    const auto background =
        orbit.variables(x, EccentricOrbit::background_tags{});
    const auto deriv_field = partial_derivative(field, mesh, inv_jacobian);
    tnsr::I<ComplexDataVector, 2> flux{};
    ScalarSelfForce::Fluxes::apply(
        make_not_null(&flux), get<Tags::Alpha>(background), {}, deriv_field);
    auto result = divergence(flux, mesh, inv_jacobian);
    get(result) *= -1.;
    ScalarSelfForce::Sources::apply(
        make_not_null(&result), get<Tags::Beta>(background),
        get<Tags::Gamma>(background), field, deriv_field, flux);
    return result;
  };
  // A low-order polynomial, so that its spectral derivatives are exact
  const Scalar<ComplexDataVector> f{
      ComplexDataVector{(1. + 0.3 * r + 0.1 * square(r)) * (1. + square(z)) +
                        std::complex<double>(0., 1.) * r * z}};
  const Scalar<ComplexDataVector> phase_times_f{
      ComplexDataVector{phase * get(f)}};
  const auto expected =
      ComplexDataVector{phase * get(apply_operator(window, f))};
  const Approx custom_approx =
      Approx::custom().epsilon(1.e-9).scale(max(abs(expected)));
  CHECK_ITERABLE_CUSTOM_APPROX(get(apply_operator(null_slicing, phase_times_f)),
                               expected, custom_approx);
}
}  // namespace

SPECTRE_TEST_CASE("Unit.PointwiseFunctions.ScalarSelfForce.EccentricOrbit",
                  "[PointwiseFunctions][Unit]") {
  test_analytic_m_modes();

  test_operator_consistency(false);
  test_operator_consistency(true);
  for (const bool inner : {true, false}) {
    test_slicing_covariance(0, inner);
    test_slicing_covariance(2, inner);
  }
}

}// namespace ScalarSelfForce::AnalyticData
