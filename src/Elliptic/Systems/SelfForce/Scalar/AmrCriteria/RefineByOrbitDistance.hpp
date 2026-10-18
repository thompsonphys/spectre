// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>
#include <cstddef>
#include <optional>
#include <pup.h>
#include <string>

#include "Domain/Amr/Flag.hpp"
#include "Domain/Creators/Tags/Domain.hpp"
#include "Domain/Domain.hpp"
#include "Domain/Structure/ElementId.hpp"
#include "Elliptic/Tags.hpp"
#include "Options/Auto.hpp"
#include "Options/Context.hpp"
#include "Options/String.hpp"
#include "Parallel/GlobalCache.hpp"
#include "ParallelAlgorithms/Amr/Criteria/Criterion.hpp"
#include "ParallelAlgorithms/Amr/Criteria/Type.hpp"
#include "PointwiseFunctions/InitialDataUtilities/Background.hpp"
#include "Utilities/Serialization/CharmPupable.hpp"
#include "Utilities/TMPL.hpp"

/// \cond
template <size_t VolumeDim>
class Block;
/// \endcond

namespace ScalarSelfForce::AmrCriteria {

namespace RefineByOrbitDistance_detail {
/*!
 * \brief The flags of `RefineByOrbitDistance` for the element `element_id` of
 * `block`, with the orbit segment at \f$\cos\theta = 0\f$ and
 * `orbit_r_min` \f$\le r \le\f$ `orbit_r_max`,
 * \f$\kappa\f$ = `max_size_over_distance` and \f$h_\mathrm{min}\f$ =
 * `min_size`, where `min_size` = 0 sets no minimum size
 */
std::array<amr::Flag, 2> flags(const Block<2>& block,
                               const ElementId<2>& element_id,
                               double orbit_r_min, double orbit_r_max,
                               double max_size_over_distance, double min_size);
}  // namespace RefineByOrbitDistance_detail

/*!
 * \brief h-refine (split) elements that are large compared to their distance
 * from the orbit
 *
 * Over a radial period the orbit sweeps the segment
 * \f$\Gamma = \{\cos\theta = 0,\ r_\mathrm{min} \le r \le r_\mathrm{max}\}\f$
 * with \f$r_\mathrm{min} = p / (1 + e)\f$ and \f$r_\mathrm{max} = p / (1 -
 * e)\f$, on which the solution is not smooth. The rate at which the spectral
 * expansion of an element converges is set by the size of the element relative
 * to its distance from \f$\Gamma\f$. This criterion splits an element in the
 * direction \f$d\f$ if
 *
 * \f[
 * h_d > \max(\kappa D, h_\mathrm{min}),
 * \f]
 *
 * where \f$\kappa\f$ is `MaxSizeOverDistance`, \f$h_\mathrm{min}\f$ is
 * `MinSize` (zero if `None`), \f$h_d\f$ is the size of the element in
 * direction \f$d\f$ and \f$D\f$ is its distance from \f$\Gamma\f$.
 * Both are measured at the point of the element nearest to \f$\Gamma\f$, with
 * the block map linearized there. The element covers the block-logical
 * coordinates \f$[\xi_\mathrm{lo}, \xi_\mathrm{hi}] \times [\eta_\mathrm{lo},
 * \eta_\mathrm{hi}]\f$, and the inverse block map takes \f$\Gamma\f$ to
 * \f$[\xi_a, \xi_b] \times \{\eta_\mathrm{eq}\}\f$, also when \f$\Gamma\f$ lies
 * outside the block. With the gaps
 *
 * \f{align}
 * g_\xi &= \max(0, \xi_a - \xi_\mathrm{hi}, \xi_\mathrm{lo} - \xi_b), &
 * g_\eta &= \max(0, \eta_\mathrm{eq} - \eta_\mathrm{hi}, \eta_\mathrm{lo} -
 * \eta_\mathrm{eq}),
 * \f}
 *
 * the nearest point \f$(\xi_n, \eta_n)\f$ of the element to \f$(\xi_a,
 * \eta_\mathrm{eq})\f$ at \f$(r, \cos\theta) = (r_n, \mu_n)\f$ and the local
 * scales
 *
 * \f{align}
 * s_\xi &= \left|\frac{\partial r}{\partial \xi}\right|, &
 * s_\eta &= \frac{r_n}{\sqrt{1 - \mu_n^2}}
 * \left|\frac{\partial \cos\theta}{\partial \eta}\right|
 * \f}
 *
 * there, the sizes and the distance are
 *
 * \f{align}
 * h_\xi &= s_\xi (\xi_\mathrm{hi} - \xi_\mathrm{lo}), &
 * h_\eta &= s_\eta (\eta_\mathrm{hi} - \eta_\mathrm{lo}), &
 * D &= \sqrt{(s_\xi g_\xi)^2 + (s_\eta g_\eta)^2}.
 * \f}
 *
 * Without a minimum size, elements that touch \f$\Gamma\f$ (\f$D = 0\f$) split
 * in both directions at every refinement step, as with `RefineAtPuncture`,
 * which is the limit \f$\kappa \to \infty\f$ of this criterion. The number of
 * elements along \f$\Gamma\f$ then doubles with every step. With a minimum
 * size they split only while they are larger than \f$h_\mathrm{min}\f$, and
 * further resolution near \f$\Gamma\f$ is left to p-refinement. In a block with
 * an affine map, \f$D\f$ is to first order in the angular offset the flat
 * distance in the metric \f$dr^2 + r^2 d\theta^2\f$. Measuring the gaps in
 * block-logical coordinates avoids splitting elements of blocks with nonlinear
 * radial maps, such as an `Inverse` distribution, that are wide in \f$r\f$ but
 * far from \f$\Gamma\f$ in the mapped coordinate.
 *
 * The criterion requires the coordinates \f$(r, \cos\theta)\f$, i.e.
 * `PenetratingHorizon` without `ImposeEquatorialSymmetry`, on the full domain
 * or on the upper half with a boundary at the equator. With the coordinate
 * \f$\cos^2\theta\f$ the distance from \f$\Gamma\f$ enters quadratically and
 * needs a different estimate.
 */
class RefineByOrbitDistance : public amr::Criterion {
 public:
  struct MaxSizeOverDistance {
    using type = double;
    static constexpr Options::String help = {
        "Split an element in a direction if its size in that direction "
        "exceeds this factor times its distance from the orbit. Must be "
        "positive."};
  };

  struct MinSize {
    using type = Options::Auto<double, Options::AutoLabel::None>;
    static constexpr Options::String help = {
        "Don't split an element in a direction if its size in that direction "
        "is at most this length, also when it touches the orbit. Must be "
        "positive. 'None' sets no minimum size, so that elements that touch "
        "the orbit split at every refinement step."};
  };

  using options = tmpl::list<MaxSizeOverDistance, MinSize>;

  static constexpr Options::String help = {
      "h-refine (split) elements that are large compared to their distance "
      "from the orbit."};

  RefineByOrbitDistance() = default;

  RefineByOrbitDistance(double max_size_over_distance,
                        std::optional<double> min_size,
                        const Options::Context& context = {});

  /// \cond
  explicit RefineByOrbitDistance(CkMigrateMessage* msg) : Criterion(msg) {}
  using PUP::able::register_constructor;
  WRAPPED_PUPable_decl_template(RefineByOrbitDistance);  // NOLINT
  /// \endcond

  amr::Criteria::Type type() override { return amr::Criteria::Type::h; }

  std::string observation_name() override { return "RefineByOrbitDistance"; }

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

  double max_size_over_distance_{};
  std::optional<double> min_size_{};
};

}  // namespace ScalarSelfForce::AmrCriteria
