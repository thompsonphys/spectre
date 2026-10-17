// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Elliptic/Systems/SelfForce/Scalar/BoundaryConditions/Neumann.hpp"

#include "DataStructures/Tensor/Tensor.hpp"

namespace ScalarSelfForce::BoundaryConditions {

Neumann::Neumann(CkMigrateMessage* m) : Base(m) {}

std::unique_ptr<domain::BoundaryConditions::BoundaryCondition>
Neumann::get_clone() const {
  return std::make_unique<Neumann>(*this);
}

void Neumann::apply(
    const gsl::not_null<Scalar<ComplexDataVector>*> /*field*/,
    const gsl::not_null<Scalar<ComplexDataVector>*> n_dot_field_gradient,
    const tnsr::i<ComplexDataVector, 2>& /*deriv_field*/) const {
  get(*n_dot_field_gradient) = 0.;
}

void Neumann::apply_linearized(
    const gsl::not_null<Scalar<ComplexDataVector>*> /*field_correction*/,
    const gsl::not_null<Scalar<ComplexDataVector>*>
        n_dot_field_gradient_correction,
    const tnsr::i<ComplexDataVector, 2>& /*deriv_field_correction*/) const {
  get(*n_dot_field_gradient_correction) = 0.;
}

bool operator==(const Neumann& /*lhs*/, const Neumann& /*rhs*/) { return true; }

bool operator!=(const Neumann& lhs, const Neumann& rhs) {
  return not(lhs == rhs);
}

PUP::able::PUP_ID Neumann::my_PUP_ID = 0;  // NOLINT

}  // namespace ScalarSelfForce::BoundaryConditions
