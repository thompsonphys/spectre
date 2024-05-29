// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <pup.h>

#include "DataStructures/ComplexDataVector.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Elliptic/Systems/SelfForce/GeneralRelativity/Tags.hpp"
#include "ParallelAlgorithms/LinearSolver/Schwarz/Tags.hpp"
#include "PointwiseFunctions/GeneralRelativity/Tags.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/MakeWithValue.hpp"
#include "Utilities/TMPL.hpp"

namespace elliptic::OptionTags {
class SchwarzSmootherGroup;
}

namespace GrSelfForce {

void fluxes(gsl::not_null<tnsr::Ijj<ComplexDataVector, 2>*> flux,
            const Scalar<ComplexDataVector>& alpha,
            const tnsr::ijj<ComplexDataVector, 2>& field_gradient);

void fluxes_on_face(gsl::not_null<tnsr::Ijj<ComplexDataVector, 2>*> flux,
                    const Scalar<ComplexDataVector>& alpha,
                    const tnsr::I<DataVector, 2>& face_normal_vector,
                    const tnsr::ii<ComplexDataVector, 2>& field);

void add_sources(gsl::not_null<tnsr::ii<ComplexDataVector, 2>*> source,
                 const Scalar<ComplexDataVector>& beta,
                 const tnsr::i<ComplexDataVector, 2>& gamma,
                 const tnsr::ii<ComplexDataVector, 2>& field,
                 const tnsr::Ijj<ComplexDataVector, 2>& flux);

struct Fluxes {
  using argument_tags = tmpl::list<Tags::Alpha>;
  using volume_tags = tmpl::list<>;
  using const_global_cache_tags = tmpl::list<>;
  static constexpr bool is_trivial = false;
  static constexpr bool is_discontinuous = false;
  static void apply(gsl::not_null<tnsr::Ijj<DataVector, 2>*> flux_re,
                    gsl::not_null<tnsr::Ijj<DataVector, 2>*> flux_im,
                    const Scalar<ComplexDataVector>& alpha,
                    const tnsr::ii<DataVector, 2>& /*field_re*/,
                    const tnsr::ii<DataVector, 2>& /*field_im*/,
                    const tnsr::ijj<DataVector, 2>& field_gradient_re,
                    const tnsr::ijj<DataVector, 2>& field_gradient_im);
  static void apply(const gsl::not_null<tnsr::Ijj<DataVector, 2>*> flux_re,
                    const gsl::not_null<tnsr::Ijj<DataVector, 2>*> flux_im,
                    const Scalar<ComplexDataVector>& alpha,
                    const tnsr::i<DataVector, 2>& /*face_normal*/,
                    const tnsr::I<DataVector, 2>& face_normal_vector,
                    const tnsr::ii<DataVector, 2>& field_re,
                    const tnsr::ii<DataVector, 2>& field_im);
};

struct Sources {
  using argument_tags = tmpl::list<Tags::Beta, Tags::Gamma>;
  static void apply(gsl::not_null<tnsr::ii<DataVector, 2>*> scalar_equation_re,
                    gsl::not_null<tnsr::ii<DataVector, 2>*> scalar_equation_im,
                    const Scalar<ComplexDataVector>& beta,
                    const tnsr::i<ComplexDataVector, 2>& gamma,
                    const tnsr::ii<DataVector, 2>& field_re,
                    const tnsr::ii<DataVector, 2>& field_im,
                    const tnsr::Ijj<DataVector, 2>& flux_re,
                    const tnsr::Ijj<DataVector, 2>& flux_im);
};

struct ModifyBoundaryData {
 private:
  static constexpr size_t Dim = 2;
  using SchwarzOptionsGroup = elliptic::OptionTags::SchwarzSmootherGroup;
  template <typename Tag>
  using overlaps_tag =
      LinearSolver::Schwarz::Tags::Overlaps<Tag, Dim, SchwarzOptionsGroup>;

 public:
  using argument_tags =
      tmpl::list<Tags::FieldIsRegularized,
                 overlaps_tag<Tags::FieldIsRegularized>, Tags::SingularField,
                 ::Tags::NormalDotFlux<Tags::SingularField>>;
  using volume_tags = tmpl::list<Tags::FieldIsRegularized,
                                 overlaps_tag<Tags::FieldIsRegularized>>;
  static void apply(
      gsl::not_null<tnsr::ii<DataVector, 2>*> field_re,
      gsl::not_null<tnsr::ii<DataVector, 2>*> field_im,
      gsl::not_null<tnsr::ii<DataVector, 2>*> n_dot_flux_re,
      gsl::not_null<tnsr::ii<DataVector, 2>*> n_dot_flux_im,
      const DirectionalId<Dim>& mortar_id, const bool sending,
      const tnsr::i<DataVector, Dim>& /*face_normal*/,
      const bool field_is_regularized,
      const DirectionalIdMap<Dim, bool>& neighbors_field_is_regularized,
      const tnsr::ii<ComplexDataVector, 2>& singular_field,
      const tnsr::ii<ComplexDataVector, 2>& singular_field_n_dot_flux) {
    if (field_is_regularized == neighbors_field_is_regularized.at(mortar_id)) {
      // Both elements solve for the same field. Nothing to do.
      return;
    }
    if (field_is_regularized) {
      // We're on an element that's regularized and sending to or receiving from
      // an element that's not regularized. We have to add or subtract the
      // singular field.
      const double sign = sending ? 1. : -1.;
      for (size_t i = 0; i < singular_field.size(); ++i) {
        (*field_re)[i] += sign * real(singular_field[i]);
        (*field_im)[i] += sign * imag(singular_field[i]);
        (*n_dot_flux_re)[i] += real(singular_field_n_dot_flux[i]);
        (*n_dot_flux_im)[i] += imag(singular_field_n_dot_flux[i]);
      }
    }
  }
};

}  // namespace ScalarSelfForce
