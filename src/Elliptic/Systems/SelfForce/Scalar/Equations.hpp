// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <pup.h>

#include "DataStructures/ComplexDataVector.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Elliptic/Systems/SelfForce/Scalar/Tags.hpp"
#include "ParallelAlgorithms/LinearSolver/Schwarz/Tags.hpp"
#include "PointwiseFunctions/GeneralRelativity/Tags.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/MakeWithValue.hpp"
#include "Utilities/TMPL.hpp"

namespace elliptic::OptionTags {
class SchwarzSmootherGroup;
}

namespace ScalarSelfForce {

void fluxes(gsl::not_null<tnsr::I<ComplexDataVector, 2>*> flux,
            const Scalar<ComplexDataVector>& alpha,
            const tnsr::i<ComplexDataVector, 2>& field_gradient);

void fluxes_on_face(gsl::not_null<tnsr::I<ComplexDataVector, 2>*> flux,
                    const Scalar<ComplexDataVector>& alpha,
                    const tnsr::I<DataVector, 2>& face_normal_vector,
                    const Scalar<ComplexDataVector>& field);

void add_sources(gsl::not_null<Scalar<ComplexDataVector>*> source,
                 const Scalar<ComplexDataVector>& beta,
                 const tnsr::i<ComplexDataVector, 2>& gamma,
                 const Scalar<ComplexDataVector>& field,
                 const tnsr::I<ComplexDataVector, 2>& flux);

struct Fluxes {
  using argument_tags = tmpl::list<Tags::Alpha>;
  using volume_tags = tmpl::list<>;
  using const_global_cache_tags = tmpl::list<>;
  static constexpr bool is_trivial = false;
  static constexpr bool is_discontinuous = false;
  static void apply(gsl::not_null<tnsr::IJ<DataVector, 2>*> flux,
                    const Scalar<ComplexDataVector>& alpha,
                    const tnsr::I<DataVector, 2>& /*field*/,
                    const tnsr::iJ<DataVector, 2>& field_gradient);
  static void apply(const gsl::not_null<tnsr::IJ<DataVector, 2>*> flux,
                    const Scalar<ComplexDataVector>& alpha,
                    const tnsr::i<DataVector, 2>& /*face_normal*/,
                    const tnsr::I<DataVector, 2>& face_normal_vector,
                    const tnsr::I<DataVector, 2>& field);
};

struct Sources {
  using argument_tags = tmpl::list<Tags::Beta, Tags::Gamma>;
  static void apply(gsl::not_null<tnsr::I<DataVector, 2>*> scalar_equation,
                    const Scalar<ComplexDataVector>& beta,
                    const tnsr::i<ComplexDataVector, 2>& gamma,
                    const tnsr::I<DataVector, 2>& field,
                    const tnsr::IJ<DataVector, 2>& flux);
};

struct ModifyBoundaryData {
 private:
  static constexpr size_t Dim = 2;
  using SchwarzOptionsGroup = elliptic::OptionTags::SchwarzSmootherGroup;
  template <typename Tag>
  using overlaps_tag =
      LinearSolver::Schwarz::Tags::Overlaps<Tag, Dim, SchwarzOptionsGroup>;

 public:
  using argument_tags = tmpl::list<Tags::FieldIsRegularized,
                                   overlaps_tag<Tags::FieldIsRegularized>>;
  static void apply(
      gsl::not_null<tnsr::I<DataVector, 2>*> field,
      gsl::not_null<tnsr::I<DataVector, 2>*> n_dot_flux,
      const DirectionalId<Dim>& mortar_id, const bool sending,
      const bool field_is_regularized,
      const DirectionalIdMap<Dim, bool>& neighbors_field_is_regularized) {
    if (field_is_regularized == neighbors_field_is_regularized.at(mortar_id)) {
      // Both elements solve for the same field. Nothing to do.
      return;
    }
    if (field_is_regularized) {
      // We're on an element that's regularized and sending to or receiving from
      // an element that's not regularized. We have to add or subtract the
      // singular field.
      const double sign = sending ? 1. : -1.;
      const size_t num_points = field->begin()->size();
      const DataVector singular_field{num_points, 1.};  // TODO
      const DataVector singular_field_n_dot_flux{num_points, 0.};
      get<0>(*field) += sign * singular_field;
      get<0>(*n_dot_flux) += sign * singular_field_n_dot_flux;
    }
  }
};

}  // namespace ScalarSelfForce
