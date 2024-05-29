// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Elliptic/Systems/SelfForce/GeneralRelativity/Equations.hpp"

#include <cstddef>

#include "DataStructures/ComplexDataVector.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"

namespace GrSelfForce {

void fluxes(const gsl::not_null<tnsr::Iaa<ComplexDataVector, 3>*> flux,
            const Scalar<ComplexDataVector>& alpha,
            const tnsr::iaa<ComplexDataVector, 3>& field_gradient) {
  // get<0>(*flux) = get<0>(field_gradient);
  // get<1>(*flux) = get(alpha) * get<1>(field_gradient);
}

void fluxes_on_face(const gsl::not_null<tnsr::Iaa<ComplexDataVector, 3>*> flux,
                    const Scalar<ComplexDataVector>& alpha,
                    const tnsr::I<DataVector, 3>& face_normal_vector,
                    const tnsr::aa<ComplexDataVector, 3>& field) {
  // get<0>(*flux) = get<0>(face_normal_vector) * get(field);
  // get<1>(*flux) = get(alpha) * get<1>(face_normal_vector) * get(field);
}

void add_sources(const gsl::not_null<tnsr::aa<ComplexDataVector, 3>*> source,
                 const tnsr::aaBB<ComplexDataVector, 3>& beta,
                 const tnsr::aaBB<ComplexDataVector, 3>& gamma_rstar,
                 const tnsr::aaBB<ComplexDataVector, 3>& gamma_theta,
                 const tnsr::aa<ComplexDataVector, 3>& field,
                 const tnsr::Iaa<ComplexDataVector, 3>& flux) {
  // get(*source) += get(beta) * get(field) + get<0>(gamma) * get<0>(flux) +
  //                 get<1>(gamma) * get<1>(flux);
}

void Fluxes::apply(const gsl::not_null<tnsr::Iaa<DataVector, 3>*> flux_re,
                   const gsl::not_null<tnsr::Iaa<DataVector, 3>*> flux_im,
                   const Scalar<ComplexDataVector>& alpha,
                   const tnsr::aa<DataVector, 3>& /*field_re*/,
                   const tnsr::aa<DataVector, 3>& /*field_im*/,
                   const tnsr::iaa<DataVector, 3>& field_gradient_re,
                   const tnsr::iaa<DataVector, 3>& field_gradient_im) {
  tnsr::Iaa<ComplexDataVector, 3> flux_complex{
      field_gradient_re.begin()->size()};
  tnsr::iaa<ComplexDataVector, 3> field_gradient_complex{
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

void Fluxes::apply(const gsl::not_null<tnsr::Iaa<DataVector, 3>*> flux_re,
                   const gsl::not_null<tnsr::Iaa<DataVector, 3>*> flux_im,
                   const Scalar<ComplexDataVector>& alpha,
                   const tnsr::i<DataVector, 3>& /*face_normal*/,
                   const tnsr::I<DataVector, 3>& face_normal_vector,
                   const tnsr::aa<DataVector, 3>& field_re,
                   const tnsr::aa<DataVector, 3>& field_im) {
  tnsr::Iaa<ComplexDataVector, 3> flux_complex{field_re.begin()->size()};
  tnsr::aa<ComplexDataVector, 3> field_complex{field_re.begin()->size()};
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

void Sources::apply(const gsl::not_null<tnsr::aa<DataVector, 3>*> sources_re,
                    const gsl::not_null<tnsr::aa<DataVector, 3>*> sources_im,
                    const tnsr::aaBB<ComplexDataVector, 3>& beta,
                    const tnsr::aaBB<ComplexDataVector, 3>& gamma_rstar,
                    const tnsr::aaBB<ComplexDataVector, 3>& gamma_theta,
                    const tnsr::aa<DataVector, 3>& field_re,
                    const tnsr::aa<DataVector, 3>& field_im,
                    const tnsr::Iaa<DataVector, 3>& flux_re,
                    const tnsr::Iaa<DataVector, 3>& flux_im) {
  tnsr::aa<ComplexDataVector, 3> sources_complex{field_re.begin()->size()};
  tnsr::aa<ComplexDataVector, 3> field_complex{field_re.begin()->size()};
  for (size_t i = 0; i < field_re.size(); ++i) {
    sources_complex[i] =
        (*sources_re)[i] + std::complex<double>(0.0, 1.0) * (*sources_im)[i];
    field_complex[i] =
        field_re[i] + std::complex<double>(0.0, 1.0) * field_im[i];
  }
  tnsr::Iaa<ComplexDataVector, 3> flux_complex{flux_re.begin()->size()};
  for (size_t i = 0; i < flux_re.size(); ++i) {
    flux_complex[i] = flux_re[i] + std::complex<double>(0.0, 1.0) * flux_im[i];
  }
  add_sources(make_not_null(&sources_complex), beta, gamma_rstar, gamma_theta,
              field_complex, flux_complex);
  for (size_t i = 0; i < sources_complex.size(); ++i) {
    (*sources_re)[i] = real(sources_complex[i]);
    (*sources_im)[i] = imag(sources_complex[i]);
  }
}

}  // namespace GrSelfForce
