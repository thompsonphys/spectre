// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <complex>
#include <cstddef>
#include <string>
#include <utility>

#include "DataStructures/ComplexDataVector.hpp"
#include "DataStructures/DataBox/DataBox.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Domain/Structure/Direction.hpp"
#include "Domain/Structure/DirectionMap.hpp"
#include "Domain/Tags.hpp"
#include "Domain/Tags/Faces.hpp"
#include "Elliptic/BoundaryConditions/ApplyBoundaryCondition.hpp"
#include "Elliptic/BoundaryConditions/BoundaryCondition.hpp"
#include "Elliptic/BoundaryConditions/BoundaryConditionType.hpp"
#include "Elliptic/Systems/SelfForce/Scalar/BoundaryConditions/Factory.hpp"
#include "Elliptic/Systems/SelfForce/Scalar/BoundaryConditions/Sommerfeld.hpp"
#include "Elliptic/Systems/SelfForce/Scalar/Tags.hpp"
#include "Framework/TestCreation.hpp"
#include "Framework/TestHelpers.hpp"
#include "NumericalAlgorithms/Spectral/Basis.hpp"
#include "NumericalAlgorithms/Spectral/CollocationPoints.hpp"
#include "NumericalAlgorithms/Spectral/Mesh.hpp"
#include "NumericalAlgorithms/Spectral/Quadrature.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/MakeWithValue.hpp"

namespace ScalarSelfForce::BoundaryConditions {

namespace {
using BC = elliptic::BoundaryConditions::BoundaryCondition<2>;

constexpr const char* base_sommerfeld =
    "Sommerfeld:\n"
    "  BlackHoleMass: 1.0\n"
    "  BlackHoleSpin: 0.5\n"
    "  SemiLatusRectum: 10.0\n"
    "  Eccentricity: 0.0\n"
    "  MModeNumber: 2\n"
    "  NModeNumber: 0\n";

// The face of `mesh` in `direction` at r = 1000, where cos(theta) = 0.5 + 0.3
// eta in the face's logical coordinate eta
tnsr::I<DataVector, 2> face_coordinates(const Mesh<2>& mesh,
                                        const Direction<2>& direction) {
  const auto& eta =
      Spectral::collocation_points(mesh.slice_away(direction.dimension()));
  return tnsr::I<DataVector, 2>{
      {{DataVector(eta.size(), 1000.), DataVector{0.5 + 0.3 * eta}}}};
}

// Apply the `boundary_condition` on the face of `mesh` in `direction` with the
// background on that face and the field and its derivative there. Returns the
// normal flux it imposes.
template <bool Linearized>
ComplexDataVector apply(const Sommerfeld& boundary_condition,
                        const Mesh<2>& mesh, const Direction<2>& direction,
                        tnsr::I<ComplexDataVector, 2> alpha,
                        Scalar<ComplexDataVector> beta,
                        tnsr::i<ComplexDataVector, 2> gamma,
                        Scalar<ComplexDataVector> field,
                        const tnsr::i<ComplexDataVector, 2>& deriv_field) {
  const auto box = db::create<db::AddSimpleTags<
      domain::Tags::Faces<2, Tags::Alpha>, domain::Tags::Faces<2, Tags::Beta>,
      domain::Tags::Faces<2, Tags::Gamma>,
      domain::Tags::Faces<2, domain::Tags::Direction<2>>,
      domain::Tags::Faces<2, domain::Tags::Coordinates<2, Frame::Inertial>>,
      domain::Tags::Mesh<2>>>(
      DirectionMap<2, tnsr::I<ComplexDataVector, 2>>{
          {direction, std::move(alpha)}},
      DirectionMap<2, Scalar<ComplexDataVector>>{{direction, std::move(beta)}},
      DirectionMap<2, tnsr::i<ComplexDataVector, 2>>{
          {direction, std::move(gamma)}},
      DirectionMap<2, Direction<2>>{{direction, direction}},
      DirectionMap<2, tnsr::I<DataVector, 2>>{
          {direction, face_coordinates(mesh, direction)}},
      mesh);
  auto n_dot_flux = make_with_value<Scalar<ComplexDataVector>>(field, 0.);
  elliptic::apply_boundary_condition<Linearized, void,
                                     standard_boundary_conditions>(
      boundary_condition, box, direction, make_not_null(&field),
      make_not_null(&n_dot_flux), deriv_field);
  return get(n_dot_flux);
}

void test_order_1_and_2() {
  const Mesh<2> mesh{
      {{3, 5}}, Spectral::Basis::Legendre, Spectral::Quadrature::Gauss};
  const Direction<2> direction = Direction<2>::lower_xi();
  const size_t num_pts = 5;
  const auto zero_deriv =
      make_with_value<tnsr::i<ComplexDataVector, 2>>(num_pts, 0.);
  {
    INFO("Apply boundary condition: Order 1");
    const auto created = TestHelpers::test_factory_creation<BC, Sommerfeld>(
        std::string{base_sommerfeld} +
        "  HyperboloidalSlicing: False\n  Order: 1");
    const auto& boundary_condition = dynamic_cast<const Sommerfeld&>(*created);
    {
      INFO("Semantics");
      test_serialization(boundary_condition);
      test_copy_semantics(boundary_condition);
      auto move_boundary_condition = boundary_condition;
      test_move_semantics(std::move(move_boundary_condition),
                          boundary_condition);
    }
    {
      INFO("Properties");
      CHECK(boundary_condition.black_hole_mass() == 1.0);
      CHECK(boundary_condition.black_hole_spin() == 0.5);
      CHECK(boundary_condition.semi_latus_rectum() == 10.0);
      CHECK(boundary_condition.eccentricity() == 0.0);
      CHECK(boundary_condition.m_mode_number() == 2);
      CHECK(boundary_condition.n_mode_number() == 0);
      CHECK(boundary_condition.hyperboloidal_slicing() == false);
      CHECK(boundary_condition.boundary_condition_types() ==
            std::vector<elliptic::BoundaryConditionType>{
                elliptic::BoundaryConditionType::Neumann});
    }
    const auto n_dot_flux = apply<false>(
        boundary_condition, mesh, direction,
        make_with_value<tnsr::I<ComplexDataVector, 2>>(num_pts, 1.),
        make_with_value<Scalar<ComplexDataVector>>(num_pts, 0.),
        make_with_value<tnsr::i<ComplexDataVector, 2>>(num_pts, 0.),
        make_with_value<Scalar<ComplexDataVector>>(num_pts, 1.), zero_deriv);
    const ComplexDataVector expected{
        num_pts, std::complex<double>(0., 0.06226111848298833)};
    CHECK_ITERABLE_APPROX(n_dot_flux, expected);
  }
  {
    INFO(
        "Apply boundary condition: Order 2 for a field without angular "
        "dependence");
    const auto created = TestHelpers::test_factory_creation<BC, Sommerfeld>(
        std::string{base_sommerfeld} +
        "  HyperboloidalSlicing: True\n  Order: 2");
    const auto& boundary_condition = dynamic_cast<const Sommerfeld&>(*created);
    auto gamma = make_with_value<tnsr::i<ComplexDataVector, 2>>(num_pts, 0.);
    get<0>(gamma) = ComplexDataVector{num_pts, {2.0, 0.0}};
    const auto n_dot_flux = apply<false>(
        boundary_condition, mesh, direction,
        make_with_value<tnsr::I<ComplexDataVector, 2>>(num_pts, 1.),
        make_with_value<Scalar<ComplexDataVector>>(
            num_pts, std::complex<double>(0.5, 0.0)),
        std::move(gamma),
        make_with_value<Scalar<ComplexDataVector>>(num_pts, 1.), zero_deriv);
    const ComplexDataVector expected(num_pts, std::complex<double>(-0.25, 0.0));
    CHECK_ITERABLE_APPROX(n_dot_flux, expected);
  }
}

void test_order_2_and_3() {
  INFO("Apply boundary condition with hyperboloidal slicing: Order 2 and 3");
  const auto created = TestHelpers::test_factory_creation<BC, Sommerfeld>(
      std::string{base_sommerfeld} +
      "  HyperboloidalSlicing: True\n  Order: 2");
  const auto& boundary_condition = dynamic_cast<const Sommerfeld&>(*created);
  CHECK(boundary_condition.order() == 2);
  // Outer radial face with 6 points in cos(theta) = z = 0.5 + 0.3 eta
  const Mesh<2> mesh{
      {{4, 6}}, Spectral::Basis::Legendre, Spectral::Quadrature::Gauss};
  const Direction<2> direction = Direction<2>::upper_xi();
  const ComplexDataVector z{get<1>(face_coordinates(mesh, direction))};
  const size_t num_pts = z.size();
  const std::complex<double> i{0., 1.};
  // A background with the structure of the null-slicing operator: the angular
  // flux coefficient vanishes at the poles
  tnsr::I<ComplexDataVector, 2> alpha{num_pts};
  get<0>(alpha) = ComplexDataVector{num_pts, 0.9};
  get<1>(alpha) = 0.01 * (1. - square(z));
  const auto beta = make_with_value<Scalar<ComplexDataVector>>(
      num_pts, std::complex<double>(0.3, 0.1));
  tnsr::i<ComplexDataVector, 2> gamma{num_pts};
  get<0>(gamma) = ComplexDataVector{num_pts, std::complex<double>(0.004, -0.1)};
  get<1>(gamma) = 0.02 * z;
  // A field u = z^3 + 2 i z, whose angular flux alpha^z du/dz is a polynomial
  // of degree 4 that the face differentiates exactly
  const Scalar<ComplexDataVector> field{
      ComplexDataVector{cube(z) + 2. * i * z}};
  tnsr::i<ComplexDataVector, 2> deriv_field{num_pts};
  get<0>(deriv_field) = ComplexDataVector{num_pts, 0.7};
  get<1>(deriv_field) = 3. * square(z) + 2. * i;
  // Order 2: the large-r equation without the radial flux divergence,
  // gamma_r du/dr = -beta u - gamma_z du/dz + d/dz(alpha^z du/dz),
  // with d/dz(alpha^z du/dz) = 0.01 (6 z - 12 z^3 - 4 i z)
  const ComplexDataVector angular_flux_divergence =
      0.01 * (6. * z - 12. * cube(z) - 4. * i * z);
  const ComplexDataVector right_hand_side =
      -get(beta) * get(field) - get<1>(gamma) * get<1>(deriv_field) +
      angular_flux_divergence;
  const ComplexDataVector expected =
      get<0>(alpha) * right_hand_side / get<0>(gamma);
  CHECK_ITERABLE_APPROX((apply<false>(boundary_condition, mesh, direction,
                                      alpha, beta, gamma, field, deriv_field)),
                        expected);
  CHECK_ITERABLE_APPROX((apply<true>(boundary_condition, mesh, direction, alpha,
                                     beta, gamma, field, deriv_field)),
                        expected);
  // Order 3 also includes the radial flux divergence to leading order,
  // -d/dr(alpha^r du/dr) = (2 / r) alpha^r du/dr, on the face at r = 1000
  const auto order_3 = TestHelpers::test_factory_creation<BC, Sommerfeld>(
      std::string{base_sommerfeld} +
      "  HyperboloidalSlicing: True\n  Order: 3");
  const ComplexDataVector expected_order_3 =
      get<0>(alpha) * right_hand_side /
      (get<0>(gamma) + 2. * get<0>(alpha) / 1000.);
  CHECK_ITERABLE_APPROX(
      (apply<false>(dynamic_cast<const Sommerfeld&>(*order_3), mesh, direction,
                    alpha, beta, gamma, field, deriv_field)),
      expected_order_3);
  CHECK_ITERABLE_APPROX(
      (apply<true>(dynamic_cast<const Sommerfeld&>(*order_3), mesh, direction,
                   alpha, beta, gamma, field, deriv_field)),
      expected_order_3);

  INFO("Order 3 needs hyperboloidal slicing, and order 4 does not exist");
  for (const auto& [order, slicing] :
       {std::pair{3, "False"}, std::pair{4, "False"}, std::pair{4, "True"}}) {
    CAPTURE(order);
    CAPTURE(slicing);
    const auto created_unsupported =
        TestHelpers::test_factory_creation<BC, Sommerfeld>(
            std::string{base_sommerfeld} + "  HyperboloidalSlicing: " +
            slicing + "\n  Order: " + std::to_string(order));
    CHECK_THROWS_WITH(
        (apply<false>(dynamic_cast<const Sommerfeld&>(*created_unsupported),
                      mesh, direction, alpha, beta, gamma, field, deriv_field)),
        Catch::Matchers::ContainsSubstring("Order " + std::to_string(order)));
  }
}
}  // namespace

SPECTRE_TEST_CASE("Unit.ScalarSelfForce.BoundaryConditions.Sommerfeld",
                  "[Unit][Elliptic]") {
  test_order_1_and_2();
  test_order_2_and_3();
}

}  // namespace ScalarSelfForce::BoundaryConditions
