// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Elliptic/Systems/SelfForce/GeneralRelativity/Equations.hpp"

#include <cstddef>

#include "DataStructures/ComplexDataVector.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"

namespace GrSelfForce {

void fluxes(const gsl::not_null<tnsr::Ijj<ComplexDataVector, 2>*> flux,
            const Scalar<ComplexDataVector>& alpha,
            const tnsr::ijj<ComplexDataVector, 2>& field_gradient) {
  // get<0>(*flux) = get<0>(field_gradient);
  // get<1>(*flux) = get(alpha) * get<1>(field_gradient);
}

void fluxes_on_face(const gsl::not_null<tnsr::Ijj<ComplexDataVector, 2>*> flux,
                    const Scalar<ComplexDataVector>& alpha,
                    const tnsr::I<DataVector, 2>& face_normal_vector,
                    const tnsr::ii<ComplexDataVector, 2>& field) {
  // get<0>(*flux) = get<0>(face_normal_vector) * get(field);
  // get<1>(*flux) = get(alpha) * get<1>(face_normal_vector) * get(field);
}

void add_sources(const gsl::not_null<tnsr::ii<ComplexDataVector, 2>*> source,
                 const Scalar<ComplexDataVector>& beta,
                 const tnsr::i<ComplexDataVector, 2>& gamma,
                 const tnsr::ii<ComplexDataVector, 2>& field,
                 const tnsr::Ijj<ComplexDataVector, 2>& flux) {
  // get(*source) += get(beta) * get(field) + get<0>(gamma) * get<0>(flux) +
  //                 get<1>(gamma) * get<1>(flux);
}

void Fluxes::apply(const gsl::not_null<tnsr::Ijj<DataVector, 2>*> flux_re,
                   const gsl::not_null<tnsr::Ijj<DataVector, 2>*> flux_im,
                   const Scalar<ComplexDataVector>& alpha,
                   const tnsr::ii<DataVector, 2>& /*field_re*/,
                   const tnsr::ii<DataVector, 2>& /*field_im*/,
                   const tnsr::ijj<DataVector, 2>& field_gradient_re,
                   const tnsr::ijj<DataVector, 2>& field_gradient_im) {
  tnsr::Ijj<ComplexDataVector, 2> flux_complex{
      field_gradient_re.begin()->size()};
  tnsr::ijj<ComplexDataVector, 2> field_gradient_complex{
      field_gradient_re.begin()->size()};
  for (size_t i = 0; i < field_gradient_re.size(); ++i) {
    field_gradient_complex[i] =
        field_gradient_re[i] +
        std::complex<double>(0.0, 1.0) * field_gradient_im[i];
  }
  fluxes(make_not_null(&flux_complex), alpha, field_gradient_complex);
  for (size_t i = 0; i < flux_complex.size(); ++i) {
    (*flux_re)[i] = real(flux_complex[i]);
    (*flux_im)[i] = imag(flux_complex[i]);
  }
}

void Fluxes::apply(const gsl::not_null<tnsr::Ijj<DataVector, 2>*> flux_re,
                   const gsl::not_null<tnsr::Ijj<DataVector, 2>*> flux_im,
                   const Scalar<ComplexDataVector>& alpha,
                   const tnsr::i<DataVector, 2>& /*face_normal*/,
                   const tnsr::I<DataVector, 2>& face_normal_vector,
                   const tnsr::ii<DataVector, 2>& field_re,
                   const tnsr::ii<DataVector, 2>& field_im) {
  tnsr::Ijj<ComplexDataVector, 2> flux_complex{field_re.begin()->size()};
  tnsr::ii<ComplexDataVector, 2> field_complex{field_re.begin()->size()};
  for (size_t i = 0; i < field_re.size(); ++i) {
    field_complex[i] =
        field_re[i] + std::complex<double>(0.0, 1.0) * field_im[i];
  }
  fluxes_on_face(make_not_null(&flux_complex), alpha, face_normal_vector,
                 field_complex);
  for (size_t i = 0; i < flux_complex.size(); ++i) {
    (*flux_re)[i] = real(flux_complex[i]);
    (*flux_im)[i] = imag(flux_complex[i]);
  }
}

void Sources::apply(const gsl::not_null<tnsr::ii<DataVector, 2>*> sources_re,
                    const gsl::not_null<tnsr::ii<DataVector, 2>*> sources_im,
                    const Scalar<ComplexDataVector>& beta,
                    const tnsr::i<ComplexDataVector, 2>& gamma,
                    const tnsr::ii<DataVector, 2>& field_re,
                    const tnsr::ii<DataVector, 2>& field_im,
                    const tnsr::Ijj<DataVector, 2>& flux_re,
                    const tnsr::Ijj<DataVector, 2>& flux_im) {
  tnsr::ii<ComplexDataVector, 2> sources_complex{field_re.begin()->size()};
  tnsr::ii<ComplexDataVector, 2> field_complex{field_re.begin()->size()};
  for (size_t i = 0; i < field_re.size(); ++i) {
    sources_complex[i] =
        (*sources_re)[i] + std::complex<double>(0.0, 1.0) * (*sources_im)[i];
    field_complex[i] =
        field_re[i] + std::complex<double>(0.0, 1.0) * field_im[i];
  }
  tnsr::Ijj<ComplexDataVector, 2> flux_complex{flux_re.begin()->size()};
  for (size_t i = 0; i < flux_re.size(); ++i) {
    flux_complex[i] = flux_re[i] + std::complex<double>(0.0, 1.0) * flux_im[i];
  }
  add_sources(make_not_null(&sources_complex), beta, gamma, field_complex,
              flux_complex);
  for (size_t i = 0; i < sources_complex.size(); ++i) {
    (*sources_re)[i] = real(sources_complex[i]);
    (*sources_im)[i] = imag(sources_complex[i]);
  }
}

}  // namespace ScalarSelfForce
