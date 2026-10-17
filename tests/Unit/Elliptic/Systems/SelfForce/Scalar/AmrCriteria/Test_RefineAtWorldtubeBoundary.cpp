// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <array>
#include <cstddef>
#include <memory>
#include <optional>
#include <string>
#include <vector>

#include "DataStructures/DataBox/DataBox.hpp"
#include "DataStructures/DataBox/ObservationBox.hpp"
#include "Domain/Amr/Flag.hpp"
#include "Domain/Block.hpp"
#include "Domain/Creators/AlignedLattice.hpp"
#include "Domain/Creators/Tags/Domain.hpp"
#include "Domain/Domain.hpp"
#include "Domain/Structure/ElementId.hpp"
#include "Domain/Structure/SegmentId.hpp"
#include "Elliptic/Systems/SelfForce/Scalar/AmrCriteria/RefineAtWorldtubeBoundary.hpp"
#include "Elliptic/Systems/SelfForce/Scalar/AnalyticData/CircularOrbit.hpp"
#include "Framework/TestCreation.hpp"
#include "Framework/TestHelpers.hpp"
#include "Options/Protocols/FactoryCreation.hpp"
#include "Parallel/GlobalCache.hpp"
#include "ParallelAlgorithms/Amr/Criteria/Criterion.hpp"
#include "Utilities/EqualWithinRoundoff.hpp"
#include "Utilities/ErrorHandling/Error.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/Serialization/RegisterDerivedClassesWithCharm.hpp"
#include "Utilities/Serialization/Serialize.hpp"
#include "Utilities/TMPL.hpp"

namespace ScalarSelfForce::AmrCriteria {

namespace {

struct Metavariables {
  static constexpr size_t volume_dim = 2;
  using component_list = tmpl::list<>;
  using const_global_cache_tags = tmpl::list<>;
  struct factory_creation
      : tt::ConformsTo<Options::protocols::FactoryCreation> {
    using factory_classes = tmpl::map<
        tmpl::pair<amr::Criterion, tmpl::list<RefineAtWorldtubeBoundary>>>;
  };
};

using background_tag =
    elliptic::Tags::Background<elliptic::analytic_data::Background>;
using box_tags = tmpl::list<background_tag, domain::Tags::Domain<2>>;

// The id of the block whose lower corner is at (r, cos_theta)
size_t block_id(const Domain<2>& domain, const double r,
                const double cos_theta) {
  for (const auto& block : domain.blocks()) {
    const auto corner = block.stationary_map()(
        tnsr::I<double, 2, Frame::BlockLogical>{{-1., -1.}});
    if (equal_within_roundoff(get<0>(corner), r) and
        equal_within_roundoff(get<1>(corner), cos_theta)) {
      return block.id();
    }
  }
  ERROR("No block has its lower corner at (" << r << ", " << cos_theta << ")");
}

// The element with refinement level 2 in both dimensions and the given
// segment indices in the block whose lower corner is at (r, cos_theta)
ElementId<2> element(const Domain<2>& domain, const double r,
                     const double cos_theta, const size_t radial_index,
                     const size_t polar_index) {
  return ElementId<2>{
      block_id(domain, r, cos_theta),
      {{SegmentId{2, radial_index}, SegmentId{2, polar_index}}}};
}

void test_criterion(const std::vector<double>& polar_block_bounds,
                    const bool refine_radially) {
  CAPTURE(polar_block_bounds);
  CAPTURE(refine_radially);
  const auto created =
      TestHelpers::test_creation<std::unique_ptr<amr::Criterion>,
                                 Metavariables>(
          "RefineAtWorldtubeBoundary:\n  RefineRadially: " +
          std::string(refine_radially ? "True" : "False"));
  REQUIRE(dynamic_cast<const RefineAtWorldtubeBoundary*>(created.get()) !=
          nullptr);
  const auto& criterion = serialize_and_deserialize(
      dynamic_cast<const RefineAtWorldtubeBoundary&>(*created));

  // Blocks in r: [2, 4], [4, 20], [20, 30]. The puncture at r = 10 is in the
  // middle column, in the block that contains the equator.
  const domain::creators::AlignedLattice<2> domain_creator{
      {{{2., 4., 20., 30.}, polar_block_bounds}},
      {},
      {},
      {{2, 2}},
      {{3, 3}},
      {},
      {},
      {}};
  auto domain = domain_creator.create_domain();
  auto databox = db::create<box_tags>(
      std::unique_ptr<elliptic::analytic_data::Background>{
          std::make_unique<ScalarSelfForce::AnalyticData::CircularOrbit>(
              1., 0.5, 10., 2,
              std::make_optional(std::array<double, 4>{{3., 3., 25., 25.}}),
              true, false)},
      std::move(domain));
  const ObservationBox<tmpl::list<>, db::DataBox<box_tags>> box{
      make_not_null(&databox)};
  Parallel::GlobalCache<Metavariables> empty_cache{};
  const auto& the_domain = db::get<domain::Tags::Domain<2>>(databox);

  const std::array<amr::Flag, 2> refine{
      {refine_radially ? amr::Flag::Split : amr::Flag::DoNothing,
       amr::Flag::Split}};
  const std::array<amr::Flag, 2> do_nothing{
      {amr::Flag::DoNothing, amr::Flag::DoNothing}};
  const auto check = [&](const ElementId<2>& element_id,
                         const std::array<amr::Flag, 2>& expected) {
    CAPTURE(element_id);
    CHECK(criterion.evaluate(box, empty_cache, element_id) == expected);
  };

  const bool has_lower_boundary = polar_block_bounds.size() == 4;
  const double worldtube_bottom =
      has_lower_boundary ? polar_block_bounds[1] : polar_block_bounds[0];
  const double worldtube_top =
      has_lower_boundary ? polar_block_bounds[2] : polar_block_bounds[1];
  {
    INFO("Elements on the upper polar boundary, inside and outside");
    for (size_t radial_index = 0; radial_index < 4; ++radial_index) {
      check(element(the_domain, 4., worldtube_bottom, radial_index, 3), refine);
      check(element(the_domain, 4., worldtube_top, radial_index, 0), refine);
    }
  }
  {
    INFO("Elements away from the polar boundaries");
    check(element(the_domain, 4., worldtube_top, 1, 1), do_nothing);
    check(element(the_domain, 4., worldtube_top, 2, 3), do_nothing);
    check(element(the_domain, 4., worldtube_bottom, 1, 1), do_nothing);
    check(element(the_domain, 4., worldtube_bottom, 2, 2), do_nothing);
  }
  {
    INFO("Elements outside the radial extent of the worldtube");
    check(element(the_domain, 2., worldtube_top, 3, 0), do_nothing);
    check(element(the_domain, 2., worldtube_bottom, 3, 3), do_nothing);
    check(element(the_domain, 20., worldtube_top, 0, 0), do_nothing);
    check(element(the_domain, 20., worldtube_bottom, 0, 3), do_nothing);
  }
  if (has_lower_boundary) {
    INFO("Elements on the lower polar boundary, inside and outside");
    check(element(the_domain, 4., worldtube_bottom, 2, 0), refine);
    check(element(the_domain, 4., polar_block_bounds[0], 2, 3), refine);
  } else {
    INFO("The lower boundary of the worldtube is the equator");
    check(element(the_domain, 4., worldtube_bottom, 2, 0), do_nothing);
  }
}

}  // namespace

SPECTRE_TEST_CASE("Unit.ScalarSelfForce.AmrCriteria.RefineAtWorldtubeBoundary",
                  "[Unit][ParallelAlgorithms]") {
  register_factory_classes_with_charm<Metavariables>();
  for (const bool refine_radially : {false, true}) {
    // Full domain in cos(theta)
    test_criterion({-1., -0.5, 0.5, 1.}, refine_radially);
    // Equatorial symmetry, where the angular coordinate is cos(theta)^2 and
    // the equator is an external boundary of the domain
    test_criterion({0., 0.25, 1.}, refine_radially);
  }
}

}  // namespace ScalarSelfForce::AmrCriteria
