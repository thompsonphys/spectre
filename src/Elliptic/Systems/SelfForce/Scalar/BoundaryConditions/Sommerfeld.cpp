// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Elliptic/Systems/SelfForce/Scalar/BoundaryConditions/Sommerfeld.hpp"

#include "DataStructures/ComplexDataVector.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/EagerMath/Magnitude.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Utilities/Gsl.hpp"

namespace ScalarSelfForce::BoundaryConditions {

Sommerfeld::Sommerfeld(const double black_hole_mass,
                       const double black_hole_spin,
                       const double orbital_radius, const int m_mode_number,
                       const bool outgoing)
    : black_hole_mass_(black_hole_mass),
      black_hole_spin_(black_hole_spin),
      orbital_radius_(orbital_radius),
      m_mode_number_(m_mode_number),
      outgoing_(outgoing) {}

void Sommerfeld::apply(
    const gsl::not_null<tnsr::I<DataVector, 2>*> field,
    const gsl::not_null<tnsr::I<DataVector, 2>*> n_dot_field_gradient,
    const tnsr::iJ<DataVector, 2>& /*deriv_field*/) const {
  // Assuming circular orbit here!
  const double a = black_hole_spin_ * black_hole_mass_;
  const double M = black_hole_mass_;
  const double r_0 = orbital_radius_;
  const double omega = 1. / (a + sqrt(cube(r_0) * M));
  Scalar<ComplexDataVector> field_complex{field->begin()->size()};
  get(field_complex) =
      get<0>(*field) + std::complex<double>(0.0, 1.0) * get<1>(*field);
  const ComplexDataVector n_dot_field_gradient_complex =
      std::complex<double>(0.0,
                           (outgoing_ ? 1.0 : -1.0) * m_mode_number_ * omega) *
      get(field_complex);
  get<0>(*n_dot_field_gradient) = real(n_dot_field_gradient_complex);
  get<1>(*n_dot_field_gradient) = imag(n_dot_field_gradient_complex);
}

void Sommerfeld::apply_linearized(
    const gsl::not_null<tnsr::I<DataVector, 2>*> field_correction,
    const gsl::not_null<tnsr::I<DataVector, 2>*>
        n_dot_field_gradient_correction,
    const tnsr::iJ<DataVector, 2>& deriv_field_correction) const {
  apply(field_correction, n_dot_field_gradient_correction,
        deriv_field_correction);
}

void Sommerfeld::pup(PUP::er& p) {
  p | black_hole_mass_;
  p | black_hole_spin_;
  p | orbital_radius_;
  p | m_mode_number_;
  p | outgoing_;
}

bool operator==(const Sommerfeld& lhs, const Sommerfeld& rhs) {
  return lhs.black_hole_mass_ == rhs.black_hole_mass_ and
         lhs.black_hole_spin_ == rhs.black_hole_spin_ and
         lhs.orbital_radius_ == rhs.orbital_radius_ and
         lhs.m_mode_number_ == rhs.m_mode_number_ and
         lhs.outgoing_ == rhs.outgoing_;
}

bool operator!=(const Sommerfeld& lhs, const Sommerfeld& rhs) {
  return not(lhs == rhs);
}

PUP::able::PUP_ID Sommerfeld::my_PUP_ID = 0;  // NOLINT

}  // namespace ScalarSelfForce::BoundaryConditions
