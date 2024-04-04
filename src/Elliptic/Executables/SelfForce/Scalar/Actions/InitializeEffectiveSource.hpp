// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <optional>
#include <utility>

#include "DataStructures/DataBox/DataBox.hpp"
#include "DataStructures/DataBox/PrefixHelpers.hpp"
#include "DataStructures/DataBox/Prefixes.hpp"
#include "DataStructures/Variables.hpp"
#include "Domain/BlockLogicalCoordinates.hpp"
#include "Domain/Creators/Tags/Domain.hpp"
#include "Domain/Domain.hpp"
#include "Domain/ElementLogicalCoordinates.hpp"
#include "Domain/Structure/ElementId.hpp"
#include "Domain/Tags.hpp"
#include "Elliptic/DiscontinuousGalerkin/Tags.hpp"
#include "NumericalAlgorithms/DiscontinuousGalerkin/ApplyMassMatrix.hpp"
#include "Parallel/AlgorithmExecution.hpp"
#include "Parallel/Tags/Metavariables.hpp"
#include "ParallelAlgorithms/Amr/Protocols/Projector.hpp"
#include "ParallelAlgorithms/Initialization/MutateAssign.hpp"
#include "ParallelAlgorithms/LinearSolver/Schwarz/Tags.hpp"
#include "PointwiseFunctions/AnalyticData/SelfForce/Scalar/CircularOrbit.hpp"
#include "Utilities/CallWithDynamicType.hpp"
#include "Utilities/MakeWithValue.hpp"
#include "Utilities/TMPL.hpp"
#include "Utilities/TaggedTuple.hpp"

/// \cond
namespace Parallel {
template <typename Metavariables>
struct GlobalCache;
}  // namespace Parallel
/// \endcond

namespace ScalarSelfForce::Actions {

/*!
 * \brief Initialize the "fixed sources" of the elliptic equations, i.e. their
 * variable-independent source term \f$f(x)\f$
 *
 * This action initializes \f$f(x)\f$ in an elliptic system of PDEs \f$-div(F) +
 * S = f(x)\f$.
 *
 * Uses:
 * - System:
 *   - `primal_fields`
 * - DataBox:
 *   - `BackgroundTag`
 *   - `Tags::Coordinates<Dim, Frame::Inertial>`
 *
 * DataBox:
 * - Adds:
 *   - `db::wrap_tags_in<::Tags::FixedSource, primal_fields>`
 */
template <typename System, typename BackgroundTag, typename SchwarzOptionsGroup,
          size_t Dim = System::volume_dim>
struct InitializeEffectiveSource : tt::ConformsTo<::amr::protocols::Projector> {
 private:
  using fixed_sources_tag = ::Tags::Variables<
      db::wrap_tags_in<::Tags::FixedSource, typename System::primal_fields>>;
  template <typename Tag>
  using overlaps_tag =
      LinearSolver::Schwarz::Tags::Overlaps<Tag, Dim, SchwarzOptionsGroup>;

 public:  // Iterable action
  using const_global_cache_tags =
      tmpl::list<elliptic::dg::Tags::Massive, BackgroundTag>;
  using simple_tags = tmpl::list<fixed_sources_tag, Tags::FieldIsRegularized,
                                 overlaps_tag<Tags::FieldIsRegularized>>;
  using compute_tags = tmpl::list<>;

  template <typename DbTagsList, typename... InboxTags, typename Metavariables,
            typename ActionList, typename ParallelComponent>
  static Parallel::iterable_action_return_t apply(
      db::DataBox<DbTagsList>& box,
      const tuples::TaggedTuple<InboxTags...>& /*inboxes*/,
      const Parallel::GlobalCache<Metavariables>& /*cache*/,
      const ElementId<Dim>& /*array_index*/, const ActionList /*meta*/,
      const ParallelComponent* const /*meta*/) {
    db::mutate_apply<InitializeEffectiveSource>(make_not_null(&box));
    return {Parallel::AlgorithmExecution::Continue, std::nullopt};
  }

 public:  // DataBox mutator, amr::protocols::Projector
  using return_tags = simple_tags;
  using argument_tags = tmpl::list<
      domain::Tags::Coordinates<Dim, Frame::Inertial>,
      domain::Tags::Domain<Dim>, domain::Tags::Element<Dim>, BackgroundTag,
      elliptic::dg::Tags::Massive, domain::Tags::Mesh<Dim>,
      domain::Tags::DetInvJacobian<Frame::ElementLogical, Frame::Inertial>,
      Parallel::Tags::Metavariables>;

  template <typename Background, typename Metavariables, typename... AmrData>
  static void apply(
      const gsl::not_null<typename fixed_sources_tag::type*> fixed_sources,
      const gsl::not_null<bool*> field_is_regularized,
      const gsl::not_null<DirectionalIdMap<Dim, bool>*>
          neighbors_field_is_regularized,
      const tnsr::I<DataVector, Dim>& inertial_coords,
      const Domain<Dim>& domain, const Element<Dim>& element,
      const Background& background, const bool massive, const Mesh<Dim>& mesh,
      const Scalar<DataVector>& det_inv_jacobian, const Metavariables& /*meta*/,
      const AmrData&... /*amr_data*/) {
    const auto& circular_orbit =
        dynamic_cast<const ScalarSelfForce::AnalyticData::CircularOrbit&>(
            background);

    // Check if this element and its neighbors solve for the regular field or
    // the full field
    const tnsr::I<double, 2> puncture_pos = circular_orbit.puncture_position();
    const auto puncture_in_element =
        [&puncture_pos, &domain](const ElementId<Dim>& element_id) -> bool {
      const auto& block = domain.blocks()[element_id.block_id()];
      const auto block_logical_coords =
          block_logical_coordinates_single_point(puncture_pos, block);
      if (not block_logical_coords.has_value()) {
        return false;
      }
      return element_logical_coordinates(*block_logical_coords, element_id)
          .has_value();
    };
    *field_is_regularized = puncture_in_element(element.id());
    neighbors_field_is_regularized->clear();
    for (const auto& [direction, neighbors] : element.neighbors()) {
      for (const auto& neighbor_id : neighbors) {
        const LinearSolver::Schwarz::OverlapId<Dim> overlap_id{direction,
                                                               neighbor_id};
        neighbors_field_is_regularized->emplace(
            overlap_id, puncture_in_element(neighbor_id));
      }
    }

    // Only set the effective source if solving for the regular field
    if (*field_is_regularized) {
      *fixed_sources = variables_from_tagged_tuple(circular_orbit.variables(
          inertial_coords, typename fixed_sources_tag::tags_list{}));
      // Apply DG mass matrix to the fixed sources if the DG operator is massive
      if (massive) {
        *fixed_sources /= get(det_inv_jacobian);
        ::dg::apply_mass_matrix(fixed_sources, mesh);
      }
    } else {
      *fixed_sources = Variables<typename fixed_sources_tag::tags_list>{
          mesh.number_of_grid_points(), 0.};
    }
  }
};

}  // namespace ScalarSelfForce::Actions
