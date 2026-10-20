// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <complex>
#include <type_traits>

#include "DataStructures/ComplexDataVector.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Elliptic/Systems/SelfForce/Scalar/Equations.hpp"
#include "Elliptic/Systems/SelfForce/Scalar/FirstOrderSystem.hpp"
#include "Elliptic/Systems/SelfForce/Scalar/Tags.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/TMPL.hpp"

SPECTRE_TEST_CASE("Unit.Elliptic.Systems.ScalarSelfForce.Equations",
                  "[Unit][Elliptic]") {
  static_assert(
      std::is_same_v<ScalarSelfForce::FirstOrderSystem::complex_shift_weight,
                     ScalarSelfForce::ComplexShiftWeight>);
  static_assert(
      std::is_same_v<ScalarSelfForce::ComplexShiftWeight::argument_tags,
                     tmpl::list<ScalarSelfForce::Tags::Beta>>);
  // The weight is the negative real part of beta where it is positive, i.e.
  // where the operator is indefinite, and zero elsewhere
  const Scalar<ComplexDataVector> beta{ComplexDataVector{
      std::complex<double>{-2., 1.}, std::complex<double>{0.5, -3.},
      std::complex<double>{0., 4.}, std::complex<double>{-1.e-3, 0.}}};
  Scalar<DataVector> weight{4_st};
  ScalarSelfForce::ComplexShiftWeight::apply(make_not_null(&weight), beta);
  CHECK_ITERABLE_APPROX(get(weight), (DataVector{2., 0., 0., 1.e-3}));
}
