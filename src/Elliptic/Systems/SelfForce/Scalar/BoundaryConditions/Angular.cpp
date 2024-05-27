// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Elliptic/Systems/SelfForce/Scalar/BoundaryConditions/Angular.hpp"

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/EagerMath/Magnitude.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Utilities/Gsl.hpp"

namespace ScalarSelfForce::BoundaryConditions {

Angular::Angular(int m_mode_number) : m_mode_number_(m_mode_number) {}

void Angular::apply(
    const gsl::not_null<tnsr::I<DataVector, 2>*> field,
    const gsl::not_null<tnsr::I<DataVector, 2>*> n_dot_field_gradient,
    const tnsr::iJ<DataVector, 2>& /*deriv_field*/) const {
  // Assert that theta = 0 or pi and that normal is in angular direction?
  if (m_mode_number_ == 0) {
    get<0>(*n_dot_field_gradient) = 0.;
    get<1>(*n_dot_field_gradient) = 0.;
  } else {
    get<0>(*field) = 0.;
    get<1>(*field) = 0.;
  }
}

void Angular::apply_linearized(
    const gsl::not_null<tnsr::I<DataVector, 2>*> field_correction,
    const gsl::not_null<tnsr::I<DataVector, 2>*>
        n_dot_field_gradient_correction,
    const tnsr::iJ<DataVector, 2>& deriv_field_correction) const {
  apply(field_correction, n_dot_field_gradient_correction,
        deriv_field_correction);
}

bool operator==(const Angular& lhs, const Angular& rhs) {
  return lhs.m_mode_number_ == rhs.m_mode_number_;
}

bool operator!=(const Angular& lhs, const Angular& rhs) {
  return not(lhs == rhs);
}

PUP::able::PUP_ID Angular::my_PUP_ID = 0;  // NOLINT

}  // namespace ScalarSelfForce::BoundaryConditions
