// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Elliptic/Systems/SelfForce/Scalar/AmrCriteria/RefineByOrbitDistance.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <optional>
#include <pup.h>

#include "DataStructures/Tensor/Tensor.hpp"
#include "Domain/Amr/Flag.hpp"
#include "Domain/Block.hpp"
#include "Domain/Domain.hpp"
#include "Domain/Structure/ElementId.hpp"
#include "Domain/Structure/Side.hpp"
#include "Elliptic/Systems/SelfForce/Scalar/AnalyticData/SelfForceBackground.hpp"
#include "Options/ParseError.hpp"
#include "PointwiseFunctions/InitialDataUtilities/Background.hpp"
#include "Utilities/ConstantExpressions.hpp"
#include "Utilities/ErrorHandling/Error.hpp"
#include "Utilities/ErrorHandling/FloatingPointExceptions.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/Serialization/PupStlCpp17.hpp"

namespace ScalarSelfForce::AmrCriteria {

namespace {
// Whether the block-logical coordinate `logical` of a point lies on the same
// side of each block boundary (at logical -1 and 1) as the point's physical
// coordinate `physical` lies of the boundary's physical coordinate
bool on_same_side_of_boundaries(const double physical, const double logical,
                                const double physical_lower,
                                const double physical_upper) {
  constexpr double tolerance = 1.e-10;
  const double orientation = physical_upper > physical_lower ? 1. : -1.;
  const bool beyond_lower = orientation * (physical - physical_lower) < 0.;
  const bool beyond_upper = orientation * (physical - physical_upper) > 0.;
  return (beyond_lower ? logical < -1. + tolerance
                       : logical > -1. - tolerance) and
         (beyond_upper ? logical > 1. - tolerance : logical < 1. + tolerance);
}
}  // namespace

namespace RefineByOrbitDistance_detail {
std::array<amr::Flag, 2> flags(const Block<2>& block,
                               const ElementId<2>& element_id,
                               const double orbit_r_min,
                               const double orbit_r_max,
                               const double max_size_over_distance,
                               const double min_size) {
  const auto& map = block.stationary_map();
  // The block-logical coordinates of the ends of the orbit segment. The orbit
  // usually lies outside the block, where the inverse map is extended
  // analytically.
  std::array<std::array<double, 2>, 2> orbit_logical{};
  for (size_t i = 0; i < 2; ++i) {
    const double r = i == 0 ? orbit_r_min : orbit_r_max;
    const auto logical = [&map, &r]() {
      // The inverse map can be undefined at the orbit, e.g. beyond the
      // singularity of a logarithmic map. Floating point exceptions are
      // disabled so that this shows up as a non-finite value below.
      const ScopedFpeState disable_fpes(false);
      return map.inverse(tnsr::I<double, 2, Frame::Inertial>{{r, 0.}});
    }();
    if (not logical.has_value() or not std::isfinite(get<0>(*logical)) or
        not std::isfinite(get<1>(*logical))) {
      ERROR("The block-logical coordinates of the orbit at r = "
            << r << ", cos(theta) = 0 are not finite in block " << block.name()
            << ". RefineByOrbitDistance needs a block map that can be "
               "inverted at the orbit.");
    }
    // The inverse map can also be finite at the orbit but beyond a singularity
    // of the map, e.g. of an Inverse distribution whose singularity lies
    // between the orbit and the block. The logical coordinate then lies beyond
    // the block boundary that is farther from the orbit.
    for (size_t d = 0; d < 2; ++d) {
      const auto boundary = [&map, &d](const double side) {
        tnsr::I<double, 2, Frame::BlockLogical> point{{0., 0.}};
        point.get(d) = side;
        return map(point).get(d);
      };
      const double physical = d == 0 ? r : 0.;
      if (not on_same_side_of_boundaries(physical, logical->get(d),
                                         boundary(-1.), boundary(1.))) {
        ERROR("The block-logical coordinate "
              << d << " of the orbit at r = " << r << ", cos(theta) = 0 is "
              << logical->get(d) << " in block " << block.name()
              << ", which spans " << boundary(-1.) << " to " << boundary(1.)
              << " in this direction. The logical coordinate lies on the "
                 "wrong side of the block, beyond a singularity of the map. "
                 "RefineByOrbitDistance needs a block map that can be "
                 "inverted at the orbit.");
      }
    }
    gsl::at(orbit_logical, i) = {{get<0>(*logical), get<1>(*logical)}};
  }
  const double xi_a = std::min(orbit_logical[0][0], orbit_logical[1][0]);
  const double xi_b = std::max(orbit_logical[0][0], orbit_logical[1][0]);
  const double eta_equator = orbit_logical[0][1];

  const double xi_lo = element_id.segment_id(0).endpoint(Side::Lower);
  const double xi_hi = element_id.segment_id(0).endpoint(Side::Upper);
  const double eta_lo = element_id.segment_id(1).endpoint(Side::Lower);
  const double eta_hi = element_id.segment_id(1).endpoint(Side::Upper);
  const double gap_xi = std::max({0., xi_a - xi_hi, xi_lo - xi_b});
  const double gap_eta =
      std::max({0., eta_equator - eta_hi, eta_lo - eta_equator});

  const tnsr::I<double, 2, Frame::BlockLogical> nearest{
      {std::clamp(xi_a, xi_lo, xi_hi),
       std::clamp(eta_equator, eta_lo, eta_hi)}};
  const auto x = map(nearest);
  const auto jacobian = map.jacobian(nearest);
  // Local lengths per unit block-logical coordinate: dr and the arc r dtheta
  const double scale_xi = std::abs(get<0, 0>(jacobian));
  const double scale_eta = get<0>(x) * std::abs(get<1, 1>(jacobian)) /
                           std::sqrt(1. - square(get<1>(x)));
  // Elements that touch the orbit segment have zero distance, so without a
  // minimum size they always split, as in RefineAtPuncture
  const double distance = std::hypot(scale_xi * gap_xi, scale_eta * gap_eta);
  const double threshold =
      std::max(max_size_over_distance * distance, min_size);
  const auto flag = [&threshold](const double size) {
    return size > threshold ? amr::Flag::Split : amr::Flag::DoNothing;
  };
  return {
      {flag(scale_xi * (xi_hi - xi_lo)), flag(scale_eta * (eta_hi - eta_lo))}};
}
}  // namespace RefineByOrbitDistance_detail

RefineByOrbitDistance::RefineByOrbitDistance(
    const double max_size_over_distance, const std::optional<double> min_size,
    const Options::Context& context)
    : max_size_over_distance_(max_size_over_distance), min_size_(min_size) {
  if (not(max_size_over_distance_ > 0.)) {
    PARSE_ERROR(context, "MaxSizeOverDistance must be positive, but is "
                             << max_size_over_distance_ << ".");
  }
  if (min_size_.has_value() and not(*min_size_ > 0.)) {
    PARSE_ERROR(context, "MinSize must be positive or 'None', but is "
                             << *min_size_ << ".");
  }
}

std::array<amr::Flag, 2> RefineByOrbitDistance::impl(
    const elliptic::analytic_data::Background& background,
    const Domain<2>& domain, const ElementId<2>& element_id) const {
  const auto* const orbit =
      dynamic_cast<const AnalyticData::SelfForceBackground*>(&background);
  if (orbit == nullptr) {
    ERROR(
        "RefineByOrbitDistance requires a SelfForceBackground, such as "
        "CircularOrbit or EccentricOrbit.");
  }
  if (not orbit->penetrating_horizon()) {
    ERROR(
        "RefineByOrbitDistance requires PenetratingHorizon, where the radial "
        "coordinate is r and not the tortoise radius.");
  }
  if (orbit->impose_equatorial_symmetry()) {
    ERROR(
        "RefineByOrbitDistance does not support ImposeEquatorialSymmetry, "
        "where the angular coordinate is cos(theta)^2.");
  }
  const double semi_latus_rectum = orbit->semi_latus_rectum();
  const double eccentricity = orbit->eccentricity();
  return RefineByOrbitDistance_detail::flags(
      domain.blocks()[element_id.block_id()], element_id,
      semi_latus_rectum / (1. + eccentricity),
      semi_latus_rectum / (1. - eccentricity), max_size_over_distance_,
      min_size_.value_or(0.));
}

void RefineByOrbitDistance::pup(PUP::er& p) {
  amr::Criterion::pup(p);
  p | max_size_over_distance_;
  p | min_size_;
}

PUP::able::PUP_ID RefineByOrbitDistance::my_PUP_ID = 0;  // NOLINT

}  // namespace ScalarSelfForce::AmrCriteria
