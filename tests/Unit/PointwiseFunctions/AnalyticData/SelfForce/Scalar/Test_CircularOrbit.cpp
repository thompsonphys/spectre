// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <array>
#include <complex>
#include <cstddef>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "DataStructures/Variables.hpp"
#include "Domain/Creators/Rectangle.hpp"
#include "Domain/ElementMap.hpp"
#include "Domain/Structure/ElementId.hpp"
#include "Elliptic/Systems/SelfForce/Scalar/Equations.hpp"
#include "NumericalAlgorithms/LinearOperators/Divergence.tpp"
#include "NumericalAlgorithms/LinearOperators/PartialDerivatives.hpp"
#include "NumericalAlgorithms/Spectral/LogicalCoordinates.hpp"
#include "PointwiseFunctions/AnalyticData/SelfForce/Scalar/CircularOrbit.hpp"
#include "Utilities/TMPL.hpp"
#include "Utilities/TaggedTuple.hpp"

namespace ScalarSelfForce::AnalyticData {

SPECTRE_TEST_CASE("Unit.PointwiseFunctions.ScalarSelfForce.CircularOrbit",
                  "[PointwiseFunctions][Unit]") {
  // Set up a rectangular grid
  const double theta_offset = M_PI / 10.;
  const double delta_theta = M_PI / 20.;
  const double rstar_offset = 0.;
  const double delta_rstar = 5.;
  const size_t npoints = 20;
  const domain::creators::Rectangle domain_creator{
      {{rstar_offset, M_PI_2 + theta_offset}},
      {{rstar_offset + delta_rstar, M_PI_2 + theta_offset + delta_theta}},
      {{0, 0}},
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
  CAPTURE(x);

  // Get the analytic fields
  const auto circular_orbit = CircularOrbit{1., 1., 0.9, 6., 0};
  CAPTURE(circular_orbit.puncture_position());
  const auto background =
      circular_orbit.variables(x, CircularOrbit::background_tags{});
  const auto& alpha = get<Tags::Alpha>(background);
  const auto& beta = get<Tags::Beta>(background);
  const auto& gamma = get<Tags::Gamma>(background);
  const auto vars = circular_orbit.variables(x, CircularOrbit::source_tags{});
  const auto& singular_field_complex = get<Tags::SingularField>(vars);
  tnsr::I<DataVector, 2> singular_field{};
  get<0>(singular_field) = real(get(singular_field_complex));
  get<1>(singular_field) = imag(get(singular_field_complex));
  const auto& deriv_singular_field_complex =
      get<::Tags::deriv<Tags::SingularField, tmpl::size_t<2>, Frame::Inertial>>(
          vars);
  tnsr::iJ<DataVector, 2> deriv_singular_field{};
  for (size_t i = 0; i < 2; ++i) {
    deriv_singular_field.get(i, 0) = real(deriv_singular_field_complex.get(i));
    deriv_singular_field.get(i, 1) = imag(deriv_singular_field_complex.get(i));
  }
  const auto& effective_source = get<::Tags::FixedSource<Tags::MMode>>(vars);

  // Take numeric derivative
  const auto numeric_deriv_singular_field =
      partial_derivative(singular_field, mesh, inv_jacobian);
  Approx custom_approx = Approx::custom().epsilon(1.e-10).scale(1.);
  for (size_t i = 0; i < deriv_singular_field.size(); ++i) {
    CAPTURE(i);
    CHECK_ITERABLE_CUSTOM_APPROX(numeric_deriv_singular_field[i],
                                 deriv_singular_field[i], custom_approx);
  }

  Variables<
      tmpl::list<::Tags::Flux<Tags::MMode, tmpl::size_t<2>, Frame::Inertial>>>
      fluxes{mesh.number_of_grid_points()};
  auto& flux_singular_field =
      get<::Tags::Flux<Tags::MMode, tmpl::size_t<2>, Frame::Inertial>>(fluxes);
  ScalarSelfForce::Fluxes::apply(make_not_null(&flux_singular_field), alpha, {},
                                 deriv_singular_field);
  auto divs = divergence(fluxes, mesh, inv_jacobian);
  auto& scalar_eqn = get<
      ::Tags::div<::Tags::Flux<Tags::MMode, tmpl::size_t<2>, Frame::Inertial>>>(
      divs);
  for (size_t i = 0; i < scalar_eqn.size(); ++i) {
    scalar_eqn[i] *= -1.;
  }
  ScalarSelfForce::Sources::apply(make_not_null(&scalar_eqn), beta, gamma,
                                  singular_field, flux_singular_field);
  for (size_t i = 0; i < scalar_eqn.size(); ++i) {
    CAPTURE(i);
    CHECK_ITERABLE_CUSTOM_APPROX(scalar_eqn[i], effective_source[i],
                                 custom_approx);
  }
}

}  // namespace ScalarSelfForce::AnalyticData
