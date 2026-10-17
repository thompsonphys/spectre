// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Elliptic/Systems/SelfForce/Scalar/AmrCriteria/RefineAtWorldtubeBoundary.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <pup.h>

#include "DataStructures/Tensor/IndexType.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Domain/Amr/Flag.hpp"
#include "Domain/Block.hpp"
#include "Domain/BlockLogicalCoordinates.hpp"
#include "Domain/Domain.hpp"
#include "Domain/Structure/Direction.hpp"
#include "Domain/Structure/ElementId.hpp"
#include "Domain/Structure/Side.hpp"
#include "Elliptic/Systems/SelfForce/Scalar/AnalyticData/SelfForceBackground.hpp"
#include "PointwiseFunctions/InitialDataUtilities/Background.hpp"
#include "Utilities/EqualWithinRoundoff.hpp"
#include "Utilities/Gsl.hpp"

namespace ScalarSelfForce::AmrCriteria {

namespace {
// The inertial lower and upper corners of the region of `block` that spans the
// given block-logical coordinates. The maps of the domain are assumed to be
// aligned with the coordinate axes, so the region is a box.
std::array<std::array<double, 2>, 2> inertial_box(
    const Block<2>& block, const std::array<double, 2>& logical_lower,
    const std::array<double, 2>& logical_upper) {
  const auto& map = block.stationary_map();
  const auto lower = map(tnsr::I<double, 2, Frame::BlockLogical>{
      {logical_lower[0], logical_lower[1]}});
  const auto upper = map(tnsr::I<double, 2, Frame::BlockLogical>{
      {logical_upper[0], logical_upper[1]}});
  std::array<std::array<double, 2>, 2> box{};
  for (size_t d = 0; d < 2; ++d) {
    gsl::at(box[0], d) = std::min(lower.get(d), upper.get(d));
    gsl::at(box[1], d) = std::max(lower.get(d), upper.get(d));
  }
  return box;
}
}  // namespace

RefineAtWorldtubeBoundary::RefineAtWorldtubeBoundary(const bool refine_radially)
    : refine_radially_(refine_radially) {}

std::array<amr::Flag, 2> RefineAtWorldtubeBoundary::impl(
    const elliptic::analytic_data::Background& background,
    const Domain<2>& domain, const ElementId<2>& element_id) const {
  const auto& orbit =
      dynamic_cast<const AnalyticData::SelfForceBackground&>(background);
  const tnsr::I<double, 2> puncture_position = orbit.puncture_position();
  const auto element_box =
      inertial_box(domain.blocks()[element_id.block_id()],
                   {{element_id.segment_id(0).endpoint(Side::Lower),
                     element_id.segment_id(1).endpoint(Side::Lower)}},
                   {{element_id.segment_id(0).endpoint(Side::Upper),
                     element_id.segment_id(1).endpoint(Side::Upper)}});
  for (const auto& block : domain.blocks()) {
    // The worldtube is the block that contains the puncture, as in
    // InitializeEffectiveSource
    if (not block_logical_coordinates_single_point(puncture_position, block)
                .has_value()) {
      continue;
    }
    const auto worldtube_box = inertial_box(block, {{-1., -1.}}, {{1., 1.}});
    const double radial_overlap =
        std::min(element_box[1][0], worldtube_box[1][0]) -
        std::max(element_box[0][0], worldtube_box[0][0]);
    if (radial_overlap <=
        1.e-10 * (worldtube_box[1][0] - worldtube_box[0][0])) {
      continue;
    }
    for (const Side side : {Side::Lower, Side::Upper}) {
      if (block.external_boundaries().contains(Direction<2>{1, side})) {
        continue;
      }
      const double polar_boundary =
          side == Side::Lower ? worldtube_box[0][1] : worldtube_box[1][1];
      if (equal_within_roundoff(element_box[0][1], polar_boundary) or
          equal_within_roundoff(element_box[1][1], polar_boundary)) {
        return {{refine_radially_ ? amr::Flag::Split : amr::Flag::DoNothing,
                 amr::Flag::Split}};
      }
    }
  }
  return {{amr::Flag::DoNothing, amr::Flag::DoNothing}};
}

void RefineAtWorldtubeBoundary::pup(PUP::er& p) {
  amr::Criterion::pup(p);
  p | refine_radially_;
}

PUP::able::PUP_ID RefineAtWorldtubeBoundary::my_PUP_ID = 0;  // NOLINT

}  // namespace ScalarSelfForce::AmrCriteria
