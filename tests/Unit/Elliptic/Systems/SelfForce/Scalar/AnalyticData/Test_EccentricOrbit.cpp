// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <array>
#include <complex>
#include <cstddef>
#include <optional>

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
#include "Utilities/TMPL.hpp"

namespace ScalarSelfForce::AnalyticData {

namespace {
void test_phi_points_rule() {
  using EccentricOrbit_detail::phi_points_for_sample;
  // Far from the particle the floor of 64 samples applies
  CHECK(phi_points_for_sample(1., 2) == 64);
  // max(m, 20) + 16 / sqrt(C1), rounded up to a power of two
  CHECK(phi_points_for_sample(1.e-2, 2) == 256);
  CHECK(phi_points_for_sample(1.e-6, 2) == 16384);
  CHECK(phi_points_for_sample(1.e-6, -2) == 16384);
  // High modes: m + 16 / sqrt(C1) = 72 and 5 m / 2 = 100
  CHECK(phi_points_for_sample(0.25, 40) == 128);
  CHECK(phi_points_for_sample(0.25, -40) == 128);
  CHECK_THROWS_WITH(phi_points_for_sample(0., 2),
                    Catch::Matchers::ContainsSubstring("is not positive"));
  CHECK_THROWS_WITH(phi_points_for_sample(-1., 2),
                    Catch::Matchers::ContainsSubstring("is not positive"));
  CHECK_THROWS_WITH(phi_points_for_sample(1.e-12, 2),
                    Catch::Matchers::ContainsSubstring("The limit is 4194304"));
}

double max_relative_difference(const ComplexDataVector& value,
                               const ComplexDataVector& reference) {
  return max(abs(value - reference)) / max(abs(reference));
}

void test_automatic_phi_points() {
  // Field points close to the sweep of the particle, where the automatic rule
  // uses many more samples near closest approach than far from it. The fixed
  // reference grid is converged at these points.
  const tnsr::I<DataVector, 2> x{{{DataVector{11.5, 13., 11.5, 13.},
                                   DataVector{1.e-2, 1.e-2, 3.e-2, 3.e-2}}}};
  const auto source_vars = [&x](const std::optional<size_t> phi_points) {
    const EccentricOrbit orbit{
        1.,    0.5,   10., 0.4,       2, 2, {{-25., -5., 20., 40.}},
        false, false, 32,  phi_points};
    return orbit.variables(x, true, EccentricOrbit::source_tags{});
  };
  const auto automatic = source_vars(std::nullopt);
  const auto reference = source_vars(8192);
  using deriv_tag =
      ::Tags::deriv<Tags::SingularField, tmpl::size_t<2>, Frame::Inertial>;
  const double source_error = max_relative_difference(
      get(get<::Tags::FixedSource<Tags::MMode>>(automatic)),
      get(get<::Tags::FixedSource<Tags::MMode>>(reference)));
  const double singular_field_error =
      max_relative_difference(get(get<Tags::SingularField>(automatic)),
                              get(get<Tags::SingularField>(reference)));
  const double deriv_r_error = max_relative_difference(
      get<0>(get<deriv_tag>(automatic)), get<0>(get<deriv_tag>(reference)));
  const double deriv_theta_error = max_relative_difference(
      get<1>(get<deriv_tag>(automatic)), get<1>(get<deriv_tag>(reference)));
  CAPTURE(source_error);
  CAPTURE(singular_field_error);
  CAPTURE(deriv_r_error);
  CAPTURE(deriv_theta_error);
  // The rule targets 1e-8 relative to the largest mode at a point. At these
  // points the error is at the roundoff level of the pointwise source.
  CHECK(source_error < 1.e-9);
  CHECK(singular_field_error < 1.e-9);
  CHECK(deriv_r_error < 1.e-9);
  CHECK(deriv_theta_error < 1.e-9);
}
}  // namespace

SPECTRE_TEST_CASE("Unit.PointwiseFunctions.ScalarSelfForce.EccentricOrbit",
                  "[PointwiseFunctions][Unit]") {
  test_phi_points_rule();
  test_automatic_phi_points();

  // This test checks both the self-force equations and the effective source
  // computation in a very robust way: it ensures that the elliptic operator
  // applied to the singular field gives the effective source.
  // This is done numerically on a rectangular grid in (r_*, cos(theta)) near
  // the puncture.

  // Set up a domain
  const double costheta_offset = 0.1;
  const double delta_costheta = 0.2;
  const double rstar_offset = 0.;
  const double delta_rstar = 5.;
  const size_t npoints = 20;
  const domain::creators::Rectangle domain_creator{
      {{rstar_offset, costheta_offset}},
      {{rstar_offset + delta_rstar, costheta_offset + delta_costheta}},
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
  const auto& cos_theta = get<1>(x);
  CAPTURE(min(r_star));
  CAPTURE(max(r_star));
  CAPTURE(min(cos_theta));
  CAPTURE(max(cos_theta));

  // Get the analytic fields
  for (int n_mode_number = 2; n_mode_number < 3; ++n_mode_number) {

    CAPTURE(n_mode_number);
    const auto eccentric_orbit = EccentricOrbit{
        1.,    0.5,   10., 0.4, 2, n_mode_number, {{-25., -5., 20., 40.}},
        false, false, 128, 128};
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
}// namespace ScalarSelfForce::AnalyticData
