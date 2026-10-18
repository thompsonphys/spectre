// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Elliptic/Systems/SelfForce/Scalar/BoundaryConditions/Sommerfeld.hpp"

#include <cmath>
#include <complex>
#include <cstddef>

#include "DataStructures/ComplexDataVector.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/EagerMath/Magnitude.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Utilities/ErrorHandling/Assert.hpp"
#include "Utilities/ErrorHandling/Error.hpp"
#include "Utilities/Gsl.hpp"

extern "C" {
#include "korb.h"
}

namespace ScalarSelfForce::BoundaryConditions {

namespace {
// The derivative at the distinct `points` of the polynomial that interpolates
// `values` there, in barycentric form
ComplexDataVector interpolating_derivative(const DataVector& points,
                                           const ComplexDataVector& values) {
  const size_t num_points = points.size();
  DataVector weights(num_points, 1.);
  for (size_t i = 0; i < num_points; ++i) {
    for (size_t j = 0; j < num_points; ++j) {
      if (j != i) {
        weights[i] /= points[i] - points[j];
      }
    }
  }
  ComplexDataVector result(num_points, 0.);
  for (size_t i = 0; i < num_points; ++i) {
    for (size_t j = 0; j < num_points; ++j) {
      if (j != i) {
        result[i] += weights[j] / weights[i] / (points[i] - points[j]) *
                     (values[j] - values[i]);
      }
    }
  }
  return result;
}
}  // namespace

Sommerfeld::Sommerfeld(const double black_hole_mass,
                       const double black_hole_spin,
                       const double semi_latus_rectum,
                       const double eccentricity,
                       const int m_mode_number, const int n_mode_number,
                       const bool hyperboloidal_slicing, const int order)
    : black_hole_mass_(black_hole_mass),
      black_hole_spin_(black_hole_spin),
      semi_latus_rectum_(semi_latus_rectum),
      eccentricity_(eccentricity),
      m_mode_number_(m_mode_number),
      n_mode_number_(n_mode_number),
      hyperboloidal_slicing_(hyperboloidal_slicing),
      order_(order) {}

Sommerfeld::Sommerfeld(CkMigrateMessage* m) : Base(m) {}

std::unique_ptr<domain::BoundaryConditions::BoundaryCondition>
Sommerfeld::get_clone() const {
  return std::make_unique<Sommerfeld>(*this);
}

void Sommerfeld::apply(
    const gsl::not_null<Scalar<ComplexDataVector>*> field,
    const gsl::not_null<Scalar<ComplexDataVector>*> n_dot_flux,
    const tnsr::i<ComplexDataVector, 2>& deriv_field,
    const tnsr::I<ComplexDataVector, 2>& alpha,
    const Scalar<ComplexDataVector>& beta,
    const tnsr::i<ComplexDataVector, 2>& gamma,
    const tnsr::I<DataVector, 2>& face_coordinates) const {
  if (hyperboloidal_slicing_) {
    if (order_ == 1) {
      get(*n_dot_flux) = 0.;
    } else if (order_ == 2 or order_ == 3) {
      ASSERT(max(abs(get<0>(face_coordinates) - get<0>(face_coordinates)[0])) <
                 1.e-10 * std::abs(get<0>(face_coordinates)[0]),
             "The Sommerfeld boundary condition of order "
                 << order_
                 << " takes angular derivatives along the boundary face, so "
                    "the face must lie at constant radius, but it spans r = "
                 << min(get<0>(face_coordinates)) << " to "
                 << max(get<0>(face_coordinates)) << ".");
      const ComplexDataVector angular_flux_divergence =
          interpolating_derivative(
              get<1>(face_coordinates),
              ComplexDataVector{get<1>(alpha) * get<1>(deriv_field)});
      // Order 3 includes the radial flux divergence to leading order in 1/r,
      // -d_r(alpha^r d_r u) = (2 / r) alpha^r d_r u
      ComplexDataVector radial_coefficient = get<0>(gamma);
      if (order_ == 3) {
        radial_coefficient += 2. * get<0>(alpha) / get<0>(face_coordinates);
      }
      get(*n_dot_flux) =
          get<0>(alpha) *
          (-get(beta) * get(*field) - get<1>(gamma) * get<1>(deriv_field) +
           angular_flux_divergence) /
          radial_coefficient;
    } else {
      ERROR("Order " << order_
                     << " not implemented for Sommerfeld boundary condition "
                        "with hyperboloidal slicing.");
    }
    return;
  }

  // Calculate frequencies associated with the orbit
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
  double Omega_r = orbpar.wr;
  double Omega_phi = orbpar.wphi;
  const double a = black_hole_spin_ * black_hole_mass_;
  const double M = black_hole_mass_;
  const double k = (m_mode_number_ * Omega_phi) + (n_mode_number_ * Omega_r);
  if (order_ == 1) {
    get(*n_dot_flux) = std::complex<double>(0.0, k) * get(*field);
  } else if (order_ == 2) {
    get(*n_dot_flux) =
        (square(k) - get(beta)) /
        (get<0>(gamma) - std::complex<double>(0.0, 2. * k)) * get(*field);
  } else {
    ERROR("Order " << order_
                   << " not implemented for Sommerfeld boundary condition.");
  }
  korb_freepar(orbpar);
}

void Sommerfeld::apply_linearized(
    const gsl::not_null<Scalar<ComplexDataVector>*> field_correction,
    const gsl::not_null<Scalar<ComplexDataVector>*> n_dot_flux_correction,
    const tnsr::i<ComplexDataVector, 2>& deriv_field_correction,
    const tnsr::I<ComplexDataVector, 2>& alpha,
    const Scalar<ComplexDataVector>& beta,
    const tnsr::i<ComplexDataVector, 2>& gamma,
    const tnsr::I<DataVector, 2>& face_coordinates) const {
  apply(field_correction, n_dot_flux_correction, deriv_field_correction, alpha,
        beta, gamma, face_coordinates);
}

void Sommerfeld::pup(PUP::er& p) {
  Base::pup(p);
  p | black_hole_mass_;
  p | black_hole_spin_;
  p | semi_latus_rectum_;
  p | eccentricity_;
  p | m_mode_number_;
  p | n_mode_number_;
  p | hyperboloidal_slicing_;
  p | order_;
}

bool operator==(const Sommerfeld& lhs, const Sommerfeld& rhs) {
  return lhs.black_hole_mass_ == rhs.black_hole_mass_ and
         lhs.black_hole_spin_ == rhs.black_hole_spin_ and
         lhs.semi_latus_rectum_ == rhs.semi_latus_rectum_ and
         lhs.eccentricity_ == rhs.eccentricity_ and
         lhs.m_mode_number_ == rhs.m_mode_number_ and
         lhs.n_mode_number_ == rhs.n_mode_number_ and
         lhs.hyperboloidal_slicing_ == rhs.hyperboloidal_slicing_ and
         lhs.order_ == rhs.order_;
}

bool operator!=(const Sommerfeld& lhs, const Sommerfeld& rhs) {
  return not(lhs == rhs);
}

PUP::able::PUP_ID Sommerfeld::my_PUP_ID = 0;  // NOLINT

}  // namespace ScalarSelfForce::BoundaryConditions
