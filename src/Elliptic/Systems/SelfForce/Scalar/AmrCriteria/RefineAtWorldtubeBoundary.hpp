// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>
#include <pup.h>
#include <string>

#include "Domain/Amr/Flag.hpp"
#include "Domain/Creators/Tags/Domain.hpp"
#include "Domain/Domain.hpp"
#include "Domain/Structure/ElementId.hpp"
#include "Elliptic/Tags.hpp"
#include "Options/String.hpp"
#include "Parallel/GlobalCache.hpp"
#include "ParallelAlgorithms/Amr/Criteria/Criterion.hpp"
#include "ParallelAlgorithms/Amr/Criteria/Type.hpp"
#include "PointwiseFunctions/InitialDataUtilities/Background.hpp"
#include "Utilities/Serialization/CharmPupable.hpp"
#include "Utilities/TMPL.hpp"

namespace ScalarSelfForce::AmrCriteria {

/*!
 * \brief h-refine (split) elements along the polar boundaries of the worldtube
 *
 * The worldtube is the block that contains the puncture, in which the system is
 * solved for the regular field instead of the full field (see
 * `ScalarSelfForce::Actions::InitializeEffectiveSource`). Across its boundary
 * the solution jumps by the singular field. This criterion splits the elements
 * that have a face on a polar boundary of the worldtube, on both sides of the
 * boundary and along the radial extent of the worldtube. It splits them in the
 * polar direction and, if `RefineRadially` is set, also in the radial
 * direction. Polar boundaries of the worldtube that are external boundaries of
 * the domain, such as the equator with equatorial symmetry, are not refined.
 */
class RefineAtWorldtubeBoundary : public amr::Criterion {
 public:
  struct RefineRadially {
    using type = bool;
    static constexpr Options::String help = {
        "Also split the elements in the radial direction. Otherwise only the "
        "polar direction is split."};
  };

  using options = tmpl::list<RefineRadially>;

  static constexpr Options::String help = {
      "h-refine (split) elements on the polar boundaries of the worldtube, the "
      "block that contains the puncture."};

  RefineAtWorldtubeBoundary() = default;

  explicit RefineAtWorldtubeBoundary(bool refine_radially);

  /// \cond
  explicit RefineAtWorldtubeBoundary(CkMigrateMessage* msg) : Criterion(msg) {}
  using PUP::able::register_constructor;
  WRAPPED_PUPable_decl_template(RefineAtWorldtubeBoundary);  // NOLINT
  /// \endcond

  amr::Criteria::Type type() override { return amr::Criteria::Type::h; }

  std::string observation_name() override {
    return "RefineAtWorldtubeBoundary";
  }

  using argument_tags = tmpl::list<
      elliptic::Tags::Background<elliptic::analytic_data::Background>,
      domain::Tags::Domain<2>>;
  using compute_tags_for_observation_box = tmpl::list<>;

  template <typename Metavariables>
  std::array<amr::Flag, 2> operator()(
      const elliptic::analytic_data::Background& background,
      const Domain<2>& domain, Parallel::GlobalCache<Metavariables>& /*cache*/,
      const ElementId<2>& element_id) const {
    return impl(background, domain, element_id);
  }

  // NOLINTNEXTLINE(google-runtime-references)
  void pup(PUP::er& p) override;

 private:
  std::array<amr::Flag, 2> impl(
      const elliptic::analytic_data::Background& background,
      const Domain<2>& domain, const ElementId<2>& element_id) const;

  bool refine_radially_{false};
};

}  // namespace ScalarSelfForce::AmrCriteria
