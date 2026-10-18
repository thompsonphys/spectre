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
#include "Domain/CoordinateMaps/Distribution.hpp"
#include "Domain/Creators/AlignedLattice.hpp"
#include "Domain/Creators/Tags/Domain.hpp"
#include "Domain/Domain.hpp"
#include "Domain/Structure/ElementId.hpp"
#include "Domain/Structure/SegmentId.hpp"
#include "Elliptic/Systems/SelfForce/Scalar/AmrCriteria/RefineAtPuncture.hpp"
#include "Elliptic/Systems/SelfForce/Scalar/AmrCriteria/RefineByOrbitDistance.hpp"
#include "Elliptic/Systems/SelfForce/Scalar/AnalyticData/CircularOrbit.hpp"
#include "Elliptic/Systems/SelfForce/Scalar/AnalyticData/EccentricOrbit.hpp"
#include "Framework/TestCreation.hpp"
#include "Framework/TestHelpers.hpp"
#include "Options/Protocols/FactoryCreation.hpp"
#include "Parallel/GlobalCache.hpp"
#include "ParallelAlgorithms/Amr/Criteria/Criterion.hpp"
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
    using factory_classes = tmpl::map<tmpl::pair<
        amr::Criterion, tmpl::list<RefineByOrbitDistance, RefineAtPuncture>>>;
  };
};

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-function"
struct OtherBackground : elliptic::analytic_data::Background,
                         elliptic::analytic_data::InitialGuess {
  OtherBackground() = default;
  explicit OtherBackground(CkMigrateMessage* /*m*/) {}
  using PUP::able::register_constructor;
  WRAPPED_PUPable_decl_template(OtherBackground);
};

PUP::able::PUP_ID OtherBackground::my_PUP_ID = 0;  // NOLINT
#pragma GCC diagnostic pop

using background_tag =
    elliptic::Tags::Background<elliptic::analytic_data::Background>;
using box_tags = tmpl::list<background_tag, domain::Tags::Domain<2>>;
using domain::CoordinateMaps::Distribution;

constexpr amr::Flag S = amr::Flag::Split;
constexpr amr::Flag N = amr::Flag::DoNothing;

// A lattice of blocks with the given bounds in r and cos(theta). The
// singularity positions must lie outside the blocks even for Linear maps.
Domain<2> lattice(const std::vector<double>& radial_bounds,
                  const std::vector<double>& polar_bounds,
                  std::array<std::vector<Distribution>, 2> distributions = {},
                  std::array<std::vector<double>, 2> singularities = {}) {
  return domain::creators::AlignedLattice<2>{{{radial_bounds, polar_bounds}},
                                             std::move(distributions),
                                             std::move(singularities),
                                             {{0, 0}},
                                             {{3, 3}},
                                             {},
                                             {},
                                             {}}
      .create_domain();
}

// The element of `block_id` with the given refinement levels and indices
ElementId<2> element(const std::array<size_t, 2>& levels,
                     const std::array<size_t, 2>& indices,
                     const size_t block_id = 0) {
  return ElementId<2>{
      block_id,
      {{SegmentId{levels[0], indices[0]}, SegmentId{levels[1], indices[1]}}}};
}

void check_flags(const Domain<2>& domain, const ElementId<2>& element_id,
                 const double kappa, const std::array<amr::Flag, 2>& expected,
                 const double min_size = 0.) {
  CAPTURE(element_id);
  CAPTURE(kappa);
  CAPTURE(min_size);
  CHECK(RefineByOrbitDistance_detail::flags(
            domain.blocks()[element_id.block_id()], element_id, 8., 12., kappa,
            min_size) == expected);
}

std::unique_ptr<elliptic::analytic_data::Background> circular_orbit(
    const bool penetrating_horizon = true,
    const bool impose_equatorial_symmetry = false) {
  return std::make_unique<AnalyticData::CircularOrbit>(
      1., 0.5, 12., 2, std::array<double, 4>{{3., 3., 25., 25.}},
      penetrating_horizon, impose_equatorial_symmetry);
}

std::unique_ptr<elliptic::analytic_data::Background> eccentric_orbit() {
  // r_min = 8.333..., r_max = 12.5
  return std::make_unique<AnalyticData::EccentricOrbit>(
      1., 0.5, 10., 0.2, 2, 0, std::array<double, 4>{{3., 3., 25., 25.}}, true,
      false, 32, std::nullopt);
}

// The flags of `criterion` for `element_id` with the given background
std::array<amr::Flag, 2> evaluate(
    const amr::Criterion& criterion,
    std::unique_ptr<elliptic::analytic_data::Background> background,
    Domain<2> domain, const ElementId<2>& element_id) {
  auto databox = db::create<box_tags>(std::move(background), std::move(domain));
  const ObservationBox<tmpl::list<>, db::DataBox<box_tags>> box{
      make_not_null(&databox)};
  Parallel::GlobalCache<Metavariables> empty_cache{};
  return criterion.evaluate(box, empty_cache, element_id);
}

void test_flags() {
  INFO("Elements near the orbit segment r in [8, 12]");
  const auto full = lattice({4., 20.}, {-0.5, 0.5});
  const auto reflected = lattice({4., 20.}, {0., 1.});
  const std::array<double, 3> kappas{{0.9, 1.5, 2.5}};
  const auto check_all =
      [&kappas](const Domain<2>& domain, const ElementId<2>& element_id,
                const std::array<std::array<amr::Flag, 2>, 3>& expected) {
        for (size_t i = 0; i < 3; ++i) {
          check_flags(domain, element_id, kappas[i], expected[i]);
        }
      };
  // r in [8, 12], cos(theta) in [0, 0.25]: on the segment
  check_all(full, element({{2, 2}}, {{1, 2}}), {{{S, S}, {S, S}, {S, S}}});
  // r in [4, 8], cos(theta) in [0, 0.25]: touches the end of the segment
  check_all(full, element({{2, 2}}, {{0, 2}}), {{{S, S}, {S, S}, {S, S}}});
  // r in [8, 12], cos(theta) in [0.25, 0.5]: h / D = 1.936 (r), 1.000 (theta)
  check_all(full, element({{2, 2}}, {{1, 3}}), {{{S, S}, {S, N}, {N, N}}});
  // r in [16, 20], cos(theta) in [0, 0.5]: h / D = 1.0 (r), 2.0 (theta)
  check_all(full, element({{2, 1}}, {{3, 1}}), {{{S, S}, {N, S}, {N, N}}});
  // The same element in the reflected domain, where the equator is the lower
  // boundary of the block
  check_all(reflected, element({{2, 1}}, {{3, 0}}), {{{S, S}, {N, S}, {N, N}}});
  // r in [16, 20], cos(theta) in [0.25, 0.5]: the nearest point is the corner
  // (16, 0.25), h / D = 0.696 (r), 0.718 (theta). The larger of the two
  // distances would give 0.968, 1.0 and their sum 0.492, 0.508. The scales at
  // the element center would give 0.636, 0.772.
  check_flags(full, element({{2, 2}}, {{3, 3}}), 0.6, {{S, S}});
  check_flags(full, element({{2, 2}}, {{3, 3}}), 0.75, {{N, N}});
  check_flags(full, element({{2, 2}}, {{3, 3}}), 0.8, {{N, N}});
  // r in [8, 12], cos(theta) in [0.5, 1]: h / D = 0.866 (r), 1.0 (theta).
  // Without the factor 1 / sin(theta) in the polar scale, h / D = 1.0 in r.
  check_flags(lattice({4., 20.}, {-1., 1.}), element({{2, 2}}, {{1, 3}}), 0.9,
              {{N, S}});
  // r in [2, 6], cos(theta) in [0, 0.25], below the segment in r: h / D = 2.0
  // (r), 0.75 (theta)
  check_all(lattice({2., 18.}, {-0.5, 0.5}), element({{2, 2}}, {{0, 2}}),
            {{{S, N}, {S, N}, {N, N}}});
}

void test_inverse_block() {
  INFO("Outer block with an Inverse distribution");
  const auto domain = lattice(
      {4., 20., 1000.}, {-0.5, 0.5},
      {{{Distribution::Linear, Distribution::Inverse}, {Distribution::Linear}}},
      {{{-51., 0.}, {-2.}}});
  // r in [20, 1000], cos(theta) in [0, 0.5]: h / D = 1.47 (r), 0.75 (theta).
  // The physical extents would give h / D = 73.5 in r.
  check_flags(domain, element({{0, 1}}, {{0, 1}}, 1), 2., {{N, N}});
  check_flags(domain, element({{0, 1}}, {{0, 1}}, 1), 1., {{S, N}});
  // r in [39.2, 1000], cos(theta) in [0, 0.5]: h / D = 0.42 (r), 0.22 (theta)
  check_flags(domain, element({{1, 1}}, {{1, 1}}, 1), 1., {{N, N}});
  // r in [20, 26.5], cos(theta) in [0.5, 1], with gaps in both directions:
  // h / D = 0.278 (r), 0.655 (theta). The scales at the element center would
  // give 0.705 in theta.
  const auto polar_domain = lattice(
      {4., 20., 1000.}, {-1., 1.},
      {{{Distribution::Linear, Distribution::Inverse}, {Distribution::Linear}}},
      {{{-51., 0.}, {-2.}}});
  check_flags(polar_domain, element({{2, 2}}, {{0, 3}}, 1), 0.62, {{N, S}});
  check_flags(polar_domain, element({{2, 2}}, {{0, 3}}, 1), 0.68, {{N, N}});
}

void test_min_size() {
  INFO("Minimum size");
  const auto full = lattice({4., 20.}, {-0.5, 0.5});
  // r in [8, 12], cos(theta) in [0, 0.25]: on the segment, D = 0, with
  // h = 4 (r) and 2 (theta) at the nearest point (8, 0)
  const auto on_segment = element({{2, 2}}, {{1, 2}});
  check_flags(full, on_segment, 0.9, {{S, S}});
  check_flags(full, on_segment, 0.9, {{S, S}}, 1.);
  check_flags(full, on_segment, 0.9, {{S, N}}, 3.);
  check_flags(full, on_segment, 0.9, {{N, N}}, 5.);
  // r in [4, 8], cos(theta) in [0, 0.25]: touches the end of the segment, with
  // h = 4 (r) and 2 (theta) at the nearest point (8, 0)
  check_flags(full, element({{2, 2}}, {{0, 2}}), 0.9, {{S, N}}, 3.);
  // r in [8, 12], cos(theta) in [0.25, 0.5]: h / D = 1.936 (r), 1.000
  // (theta) with h = 4 (r) and 2.066 (theta). The minimum size stops the
  // angular split that kappa D alone would request.
  const auto above_segment = element({{2, 2}}, {{1, 3}});
  check_flags(full, above_segment, 0.9, {{S, S}}, 1.);
  check_flags(full, above_segment, 0.9, {{S, N}}, 3.);
  // r in [16, 20], cos(theta) in [0, 0.5]: h / D = 1.0 (r), 2.0 (theta). A
  // minimum size below kappa D leaves the flags unchanged.
  check_flags(full, element({{2, 1}}, {{3, 1}}), 1.5, {{N, S}}, 0.1);
}

void test_large_kappa() {
  INFO("Large kappa splits exactly the elements that RefineAtPuncture splits");
  const auto make_domain = []() { return lattice({4., 20.}, {-0.5, 0.5}); };
  const auto element_ids = []() {
    std::vector<ElementId<2>> result{};
    for (size_t radial_index = 0; radial_index < 4; ++radial_index) {
      for (size_t polar_index = 0; polar_index < 4; ++polar_index) {
        result.push_back(element({{2, 2}}, {{radial_index, polar_index}}));
      }
    }
    return result;
  }();
  const auto domain = make_domain();
  for (const auto& element_id : element_ids) {
    const bool touches_segment = element_id.segment_id(0).index() < 3 and
                                 (element_id.segment_id(1).index() == 1 or
                                  element_id.segment_id(1).index() == 2);
    check_flags(domain, element_id, 1.e3,
                touches_segment ? std::array{S, S} : std::array{N, N});
  }
  const RefineByOrbitDistance criterion{1.e3, std::nullopt};
  const RefineAtPuncture refine_at_puncture{};
  for (const bool eccentric : {false, true}) {
    CAPTURE(eccentric);
    auto databox = db::create<box_tags>(
        eccentric ? eccentric_orbit() : circular_orbit(), make_domain());
    const ObservationBox<tmpl::list<>, db::DataBox<box_tags>> box{
        make_not_null(&databox)};
    Parallel::GlobalCache<Metavariables> empty_cache{};
    size_t num_split = 0;
    for (const auto& element_id : element_ids) {
      CAPTURE(element_id);
      const auto flags = criterion.evaluate(box, empty_cache, element_id);
      CHECK(flags == refine_at_puncture.evaluate(box, empty_cache, element_id));
      if (flags == std::array{S, S}) {
        ++num_split;
      }
    }
    CHECK(num_split == 4);
  }
}

void test_criterion() {
  INFO("Creation and evaluation through the background");
  const auto created =
      TestHelpers::test_creation<std::unique_ptr<amr::Criterion>,
                                 Metavariables>(
          "RefineByOrbitDistance:\n  MaxSizeOverDistance: 1.5\n"
          "  MinSize: None");
  REQUIRE(dynamic_cast<const RefineByOrbitDistance*>(created.get()) != nullptr);
  const auto criterion = serialize_and_deserialize(
      dynamic_cast<const RefineByOrbitDistance&>(*created));
  // r in [8, 12], cos(theta) in [0.25, 0.5]: h / D = 1.859 (r), 1.000 (theta)
  CHECK(evaluate(criterion, eccentric_orbit(), lattice({4., 20.}, {-0.5, 0.5}),
                 element({{2, 2}}, {{1, 3}})) == std::array{S, N});
  // r in [8, 12], cos(theta) in [0, 0.25]: on the segment, which starts at
  // r_min = 8.33, with h = 4 (r) and 2.08 (theta) at the nearest point
  // (8.33, 0)
  CHECK(evaluate(criterion, eccentric_orbit(), lattice({4., 20.}, {-0.5, 0.5}),
                 element({{2, 2}}, {{1, 2}})) == std::array{S, S});
  const auto created_with_min_size =
      TestHelpers::test_creation<std::unique_ptr<amr::Criterion>,
                                 Metavariables>(
          "RefineByOrbitDistance:\n  MaxSizeOverDistance: 1.5\n"
          "  MinSize: 3.");
  REQUIRE(dynamic_cast<const RefineByOrbitDistance*>(
              created_with_min_size.get()) != nullptr);
  const auto criterion_with_min_size = serialize_and_deserialize(
      dynamic_cast<const RefineByOrbitDistance&>(*created_with_min_size));
  CHECK(evaluate(criterion_with_min_size, eccentric_orbit(),
                 lattice({4., 20.}, {-0.5, 0.5}),
                 element({{2, 2}}, {{1, 2}})) == std::array{S, N});
}

void test_errors() {
  INFO("Errors");
  const RefineByOrbitDistance criterion{1., std::nullopt};
  const auto evaluate_with =
      [&criterion](
          std::unique_ptr<elliptic::analytic_data::Background> background) {
        return evaluate(criterion, std::move(background),
                        lattice({4., 20.}, {-0.5, 0.5}),
                        element({{2, 2}}, {{1, 3}}));
      };
  CHECK_THROWS_WITH(
      evaluate_with(std::make_unique<OtherBackground>()),
      Catch::Matchers::ContainsSubstring(
          "RefineByOrbitDistance requires a SelfForceBackground"));
  CHECK_THROWS_WITH(
      evaluate_with(circular_orbit(true, true)),
      Catch::Matchers::ContainsSubstring("ImposeEquatorialSymmetry"));
  CHECK_THROWS_WITH(evaluate_with(circular_orbit(false, false)),
                    Catch::Matchers::ContainsSubstring("PenetratingHorizon"));
  CHECK_THROWS_WITH((TestHelpers::test_creation<std::unique_ptr<amr::Criterion>,
                                                Metavariables>(
                        "RefineByOrbitDistance:\n  MaxSizeOverDistance: 0.\n"
                        "  MinSize: None")),
                    Catch::Matchers::ContainsSubstring(
                        "MaxSizeOverDistance must be positive"));
  CHECK_THROWS_WITH(
      (TestHelpers::test_creation<std::unique_ptr<amr::Criterion>,
                                  Metavariables>(
          "RefineByOrbitDistance:\n  MaxSizeOverDistance: 1.\n"
          "  MinSize: 0.")),
      Catch::Matchers::ContainsSubstring("MinSize must be positive"));
  // The logarithmic map of the outer block is undefined at the orbit, which
  // lies beyond its singularity at r = 15
  const auto logarithmic =
      lattice({4., 20., 1000.}, {-0.5, 0.5},
              {{{Distribution::Linear, Distribution::Logarithmic},
                {Distribution::Linear}}},
              {{{-51., 15.}, {-2.}}});
  CHECK_THROWS_WITH(RefineByOrbitDistance_detail::flags(
                        logarithmic.blocks()[1], element({{0, 0}}, {{0, 0}}, 1),
                        8., 12., 1., 0.),
                    Catch::Matchers::ContainsSubstring("not finite"));
  // The inverse map of an Inverse distribution with its singularity at r = 15,
  // between the orbit and the block, is finite at the orbit but places it
  // beyond the outer boundary of the block
  const auto inverse = lattice(
      {4., 20., 1000.}, {-0.5, 0.5},
      {{{Distribution::Linear, Distribution::Inverse}, {Distribution::Linear}}},
      {{{-51., 15.}, {-2.}}});
  CHECK_THROWS_WITH(
      RefineByOrbitDistance_detail::flags(
          inverse.blocks()[1], element({{0, 0}}, {{0, 0}}, 1), 8., 12., 1., 0.),
      Catch::Matchers::ContainsSubstring("on the wrong side"));
}

}  // namespace

SPECTRE_TEST_CASE("Unit.ScalarSelfForce.AmrCriteria.RefineByOrbitDistance",
                  "[Unit][ParallelAlgorithms]") {
  register_factory_classes_with_charm<Metavariables>();
  test_flags();
  test_inverse_block();
  test_min_size();
  test_large_kappa();
  test_criterion();
  test_errors();
}

}  // namespace ScalarSelfForce::AmrCriteria
