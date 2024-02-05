// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Elliptic/Systems/SelfForce/Scalar/Equations.hpp"

#include <cstddef>

#include "DataStructures/ComplexDataVector.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"

namespace ScalarSelfForce {

void fluxes(const gsl::not_null<tnsr::I<ComplexDataVector, 2>*> flux,
            const Scalar<ComplexDataVector>& alpha,
            const tnsr::i<ComplexDataVector, 2>& field_gradient) {
  get<0>(*flux) = get<0>(field_gradient);
  get<1>(*flux) = get(alpha) * get<1>(field_gradient);
}

void fluxes_on_face(const gsl::not_null<tnsr::I<ComplexDataVector, 2>*> flux,
                    const Scalar<ComplexDataVector>& alpha,
                    const tnsr::I<DataVector, 2>& face_normal_vector,
                    const Scalar<ComplexDataVector>& field) {
  get<0>(*flux) = get<0>(face_normal_vector) * get(field);
  get<1>(*flux) = get(alpha) * get<1>(face_normal_vector) * get(field);
}

void add_sources(const gsl::not_null<Scalar<ComplexDataVector>*> source,
                 const Scalar<ComplexDataVector>& beta,
                 const tnsr::i<ComplexDataVector, 2>& gamma,
                 const Scalar<ComplexDataVector>& field,
                 const tnsr::I<ComplexDataVector, 2>& flux) {
  // question: isn't just the Poisson equation on a curved background, where
  // the metric is the Jacobian of the coordinate transformation?
  get(*source) += get(beta) * get(field) + get<0>(gamma) * get<0>(flux) +
                  get<1>(gamma) * get<1>(flux);
}

void Fluxes::apply(const gsl::not_null<tnsr::IJ<DataVector, 2>*> flux,
                   const Scalar<ComplexDataVector>& alpha,
                   const tnsr::I<DataVector, 2>& /*field*/,
                   const tnsr::iJ<DataVector, 2>& field_gradient) {
  tnsr::I<ComplexDataVector, 2> flux_complex{field_gradient.begin()->size()};
  tnsr::i<ComplexDataVector, 2> field_gradient_complex{
      field_gradient.begin()->size()};
  for (size_t d = 0; d < 2; ++d) {
    field_gradient_complex.get(d) =
        field_gradient.get(d, 0) +
        std::complex<double>(0.0, 1.0) * field_gradient.get(d, 1);
  }
  fluxes(make_not_null(&flux_complex), alpha, field_gradient_complex);
  get<0, 0>(*flux) = real(get<0>(flux_complex));
  get<0, 1>(*flux) = imag(get<0>(flux_complex));
  get<1, 0>(*flux) = real(get<1>(flux_complex));
  get<1, 1>(*flux) = imag(get<1>(flux_complex));
}

void Fluxes::apply(const gsl::not_null<tnsr::IJ<DataVector, 2>*> flux,
                   const Scalar<ComplexDataVector>& alpha,
                   const tnsr::i<DataVector, 2>& /*face_normal*/,
                   const tnsr::I<DataVector, 2>& face_normal_vector,
                   const tnsr::I<DataVector, 2>& field) {
  tnsr::I<ComplexDataVector, 2> flux_complex{field.begin()->size()};
  Scalar<ComplexDataVector> field_complex{field.begin()->size()};
  get(field_complex) =
      get<0>(field) + std::complex<double>(0.0, 1.0) * get<1>(field);
  fluxes_on_face(make_not_null(&flux_complex), alpha, face_normal_vector,
                 field_complex);
  get<0, 0>(*flux) = real(get<0>(flux_complex));
  get<0, 1>(*flux) = imag(get<0>(flux_complex));
  get<1, 0>(*flux) = real(get<1>(flux_complex));
  get<1, 1>(*flux) = imag(get<1>(flux_complex));
}

void Sources::apply(
    const gsl::not_null<tnsr::I<DataVector, 2>*> scalar_equation,
    const Scalar<ComplexDataVector>& beta,
    const tnsr::i<ComplexDataVector, 2>& gamma,
    const tnsr::I<DataVector, 2>& field, const tnsr::IJ<DataVector, 2>& flux) {
  Scalar<ComplexDataVector> scalar_equation_complex{field.begin()->size()};
  get(scalar_equation_complex) =
      get<0>(*scalar_equation) +
      std::complex<double>(0.0, 1.0) * get<1>(*scalar_equation);
  Scalar<ComplexDataVector> field_complex{field.begin()->size()};
  get(field_complex) =
      get<0>(field) + std::complex<double>(0.0, 1.0) * get<1>(field);
  tnsr::I<ComplexDataVector, 2> flux_complex{flux.begin()->size()};
  for (size_t d = 0; d < 2; ++d) {
    flux_complex.get(d) =
        flux.get(d, 0) + std::complex<double>(0.0, 1.0) * flux.get(d, 1);
  }
  add_sources(make_not_null(&scalar_equation_complex), beta, gamma,
              field_complex, flux_complex);
  get<0>(*scalar_equation) = real(get(scalar_equation_complex));
  get<1>(*scalar_equation) = imag(get(scalar_equation_complex));
}

}  // namespace ScalarSelfForce
