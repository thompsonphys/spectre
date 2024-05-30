// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>
#include <cstddef>
#include <iosfwd>
#include <limits>
#include <memory>
#include <vector>

#include "DataStructures/Index.hpp"
#include "Domain/BoundaryConditions/BoundaryCondition.hpp"
#include "Domain/BoundaryConditions/GetBoundaryConditionsBase.hpp"
#include "Domain/Creators/DomainCreator.hpp"
#include "Domain/Domain.hpp"
#include "Domain/Structure/DirectionMap.hpp"
#include "Options/Context.hpp"
#include "Options/String.hpp"
#include "Utilities/MakeArray.hpp"
#include "Utilities/TMPL.hpp"

/// \cond
namespace domain {
namespace CoordinateMaps {
class Affine;
template <typename Map1, typename Map2>
class ProductOf2Maps;
template <typename Map1, typename Map2, typename Map3>
class ProductOf3Maps;
}  // namespace CoordinateMaps

template <typename SourceFrame, typename TargetFrame, typename... Maps>
class CoordinateMap;
}  // namespace domain
/// \endcond

namespace domain {
namespace creators {

template <size_t VolumeDim>
struct RefinementRegion {
  std::array<size_t, VolumeDim> lower_corner_index;
  std::array<size_t, VolumeDim> upper_corner_index;
  std::array<size_t, VolumeDim> refinement;

  struct LowerCornerIndex {
    using type = std::array<size_t, VolumeDim>;
    static constexpr Options::String help = {"Lower bound of refined region."};
  };

  struct UpperCornerIndex {
    using type = std::array<size_t, VolumeDim>;
    static constexpr Options::String help = {"Upper bound of refined region."};
  };

  struct Refinement {
    using type = std::array<size_t, VolumeDim>;
    static constexpr Options::String help = {"Refinement inside region."};
  };

  static constexpr Options::String help = {
      "A region to be refined differently from the default for the lattice.\n"
      "The region is a box between the block boundaries indexed by the\n"
      "Lower- and UpperCornerIndex options."};
  using options = tmpl::list<LowerCornerIndex, UpperCornerIndex, Refinement>;
  RefinementRegion(const std::array<size_t, VolumeDim>& lower_corner_index_in,
                   const std::array<size_t, VolumeDim>& upper_corner_index_in,
                   const std::array<size_t, VolumeDim>& refinement_in)
      : lower_corner_index(lower_corner_index_in),
        upper_corner_index(upper_corner_index_in),
        refinement(refinement_in) {}
  RefinementRegion() = default;
};

/// \cond
// This is needed to print the default value for the RefinedGridPoints
// option.  Since the default value is an empty vector, this function
// is never actually called.
template <size_t VolumeDim>
[[noreturn]] std::ostream& operator<<(
    std::ostream& /*s*/, const RefinementRegion<VolumeDim>& /*unused*/);
/// \endcond

/// \brief Create a Domain consisting of multiple aligned Blocks arrayed in a
/// lattice.
///
/// This is useful for setting up problems with piecewise smooth initial data,
/// problems that specify different boundary conditions on distinct parts of
/// the boundary, or problems that need different length scales initially.
///
/// \note Adaptive mesh refinement can never join Block%s, so use the fewest
/// number of Block%s that your problem needs.  More initial Element%s can be
/// created by specifying a larger `InitialRefinement`.
template <size_t VolumeDim>
class AlignedLattice : public DomainCreator<VolumeDim> {
 public:
  using maps_list = tmpl::list<
      domain::CoordinateMap<Frame::BlockLogical, Frame::Inertial,
                            CoordinateMaps::Affine>,
      domain::CoordinateMap<
          Frame::BlockLogical, Frame::Inertial,
          CoordinateMaps::ProductOf2Maps<CoordinateMaps::Affine,
                                         CoordinateMaps::Affine>>,
      domain::CoordinateMap<Frame::BlockLogical, Frame::Inertial,
                            CoordinateMaps::ProductOf3Maps<
                                CoordinateMaps::Affine, CoordinateMaps::Affine,
                                CoordinateMaps::Affine>>>;

  struct BlockBounds {
    using type = std::array<std::vector<double>, VolumeDim>;
    static constexpr Options::String help = {
        "Coordinates of block boundaries in each dimension."};
  };

  struct IsPeriodicIn {
    using type = std::array<bool, VolumeDim>;
    static constexpr Options::String help = {
        "Whether the domain is periodic in each dimension."};
  };

  struct InitialLevels {
    using type = std::array<size_t, VolumeDim>;
    static constexpr Options::String help = {
        "Initial refinement level in each dimension."};
  };

  struct InitialGridPoints {
    using type = std::array<size_t, VolumeDim>;
    static constexpr Options::String help = {
        "Initial number of grid points in each dimension."};
  };

  struct RefinedLevels {
    using type = std::vector<RefinementRegion<VolumeDim>>;
    static constexpr Options::String help = {
        "h-refined regions.  Later entries take priority."};
  };

  struct RefinedGridPoints {
    using type = std::vector<RefinementRegion<VolumeDim>>;
    static constexpr Options::String help = {
        "p-refined regions.  Later entries take priority."};
  };

  struct BlocksToExclude {
    using type = std::vector<std::array<size_t, VolumeDim>>;
    static constexpr Options::String help = {
        "List of Block indices to exclude, if any."};
  };

  template <typename BoundaryConditionsBase>
  struct LowerUpperBoundaryCondition {
    static constexpr Options::String help =
        "Lower and upper Boundary Conditions";
    struct LowerBC {
      using type = std::unique_ptr<BoundaryConditionsBase>;
      static constexpr Options::String help = "Lower Boundary Condition";
      static std::string name() { return "LowerBoundaryCondition"; };
    };
    struct UpperBC {
      using type = std::unique_ptr<BoundaryConditionsBase>;
      static constexpr Options::String help = "Upper Boundary Condition";
      static std::string name() { return "UpperBoundaryCondition"; };
    };
    LowerUpperBoundaryCondition(typename LowerBC::type lower_bc,
                                typename UpperBC::type upper_bc)
        : lower(std::move(lower_bc)), upper(std::move(upper_bc)) {};
    LowerUpperBoundaryCondition() = default;
    std::unique_ptr<BoundaryConditionsBase> lower;
    std::unique_ptr<BoundaryConditionsBase> upper;
    using options = tmpl::list<LowerBC, UpperBC>;
  };

  template <typename BoundaryConditionsBase, size_t Dim>
  struct BoundaryCondition {
    static std::string name() {
      return "BoundaryConditionIn" +
             std::string{Dim == 0 ? 'X' : (Dim == 1 ? 'Y' : 'Z')};
    };
    static constexpr Options::String help = {
        "The boundary condition to be imposed for boundaries at "
        "the given direction. Either specify one B.C. to be imposed for both "
        "lower and upper boundary or a pair of B.C. for lower and upper "
        "boundary respectively."};
    using type =
        std::variant<std::unique_ptr<BoundaryConditionsBase>,
                     LowerUpperBoundaryCondition<BoundaryConditionsBase>>;
  };

  using common_options =
      tmpl::list<BlockBounds, InitialLevels, InitialGridPoints, RefinedLevels,
                 RefinedGridPoints, BlocksToExclude>;
  using options_periodic = tmpl::list<IsPeriodicIn>;

  template <typename Metavariables>
  using options = tmpl::append<
      common_options,
      tmpl::conditional_t<
          domain::BoundaryConditions::has_boundary_conditions_base_v<
              typename Metavariables::system>,
          tmpl::list<BoundaryCondition<domain::BoundaryConditions::
                                           get_boundary_conditions_base<
                                               typename Metavariables::system>,
                                       0>,
                     BoundaryCondition<domain::BoundaryConditions::
                                           get_boundary_conditions_base<
                                               typename Metavariables::system>,
                                       1>>,
          options_periodic>>;

  static constexpr Options::String help = {
      "AlignedLattice creates a regular lattice of blocks whose corners are\n"
      "given by tensor products of the specified BlockBounds. Each Block in\n"
      "the lattice is identified by a VolumeDim-tuple of zero-based indices\n"
      "Supplying a list of these tuples to BlocksToExclude will result in\n"
      "the domain having the corresponding Blocks excluded. See the Domain\n"
      "Creation tutorial in the documentation for more information on Block\n"
      "numberings in rectilinear domains. Note that if any Blocks are\n"
      "excluded, setting the option IsPeriodicIn to `true` in any dimension\n"
      "will trigger an error, as periodic boundary\n"
      "conditions for this domain with holes is not supported."};

  AlignedLattice(typename BlockBounds::type block_bounds,
                 typename InitialLevels::type initial_refinement_levels,
                 typename InitialGridPoints::type initial_number_of_grid_points,
                 typename RefinedLevels::type refined_refinement,
                 typename RefinedGridPoints::type refined_grid_points,
                 typename BlocksToExclude::type blocks_to_exclude,
                 typename IsPeriodicIn::type is_periodic_in,
                 const Options::Context& context = {});

  AlignedLattice(typename BlockBounds::type block_bounds,
                 typename InitialLevels::type initial_refinement_levels,
                 typename InitialGridPoints::type initial_number_of_grid_points,
                 typename RefinedLevels::type refined_refinement,
                 typename RefinedGridPoints::type refined_grid_points,
                 typename BlocksToExclude::type blocks_to_exclude,
                 std::unique_ptr<domain::BoundaryConditions::BoundaryCondition>
                     boundary_condition_in_lower_x,
                 std::unique_ptr<domain::BoundaryConditions::BoundaryCondition>
                     boundary_condition_in_upper_x,
                 std::unique_ptr<domain::BoundaryConditions::BoundaryCondition>
                     boundary_condition_in_lower_y,
                 std::unique_ptr<domain::BoundaryConditions::BoundaryCondition>
                     boundary_condition_in_upper_y,
                 const Options::Context& context = {});

  template <typename BoundaryConditionsBase>
  AlignedLattice(
      typename BlockBounds::type block_bounds,
      typename InitialLevels::type initial_refinement_levels,
      typename InitialGridPoints::type initial_number_of_grid_points,
      typename RefinedLevels::type refined_refinement,
      typename RefinedGridPoints::type refined_grid_points,
      typename BlocksToExclude::type blocks_to_exclude,
      std::variant<std::unique_ptr<BoundaryConditionsBase>,
                   LowerUpperBoundaryCondition<BoundaryConditionsBase>>
          boundary_conditions_in_x,
      std::variant<std::unique_ptr<BoundaryConditionsBase>,
                   LowerUpperBoundaryCondition<BoundaryConditionsBase>>
          boundary_conditions_in_y,
      const Options::Context& context = {});

  AlignedLattice() = default;
  AlignedLattice(const AlignedLattice&) = delete;
  AlignedLattice(AlignedLattice&&) = default;
  AlignedLattice& operator=(const AlignedLattice&) = delete;
  AlignedLattice& operator=(AlignedLattice&&) = default;
  ~AlignedLattice() override = default;

  Domain<VolumeDim> create_domain() const override;

  std::vector<DirectionMap<
      VolumeDim,
      std::unique_ptr<domain::BoundaryConditions::BoundaryCondition>>>
  external_boundary_conditions() const override;

  std::vector<std::array<size_t, VolumeDim>> initial_extents() const override;

  std::vector<std::array<size_t, VolumeDim>> initial_refinement_levels()
      const override;

 private:
  typename BlockBounds::type block_bounds_{
      make_array<VolumeDim, std::vector<double>>({})};
  typename IsPeriodicIn::type is_periodic_in_{make_array<VolumeDim>(false)};
  typename InitialLevels::type initial_refinement_levels_{
      make_array<VolumeDim>(std::numeric_limits<size_t>::max())};
  typename InitialGridPoints::type initial_number_of_grid_points_{
      make_array<VolumeDim>(std::numeric_limits<size_t>::max())};
  typename RefinedLevels::type refined_refinement_{};
  typename RefinedGridPoints::type refined_grid_points_{};
  typename BlocksToExclude::type blocks_to_exclude_{};
  Index<VolumeDim> number_of_blocks_by_dim_{};
  std::unique_ptr<domain::BoundaryConditions::BoundaryCondition>
      boundary_condition_in_lower_x_;
  std::unique_ptr<domain::BoundaryConditions::BoundaryCondition>
      boundary_condition_in_upper_x_;
  std::unique_ptr<domain::BoundaryConditions::BoundaryCondition>
      boundary_condition_in_lower_y_;
  std::unique_ptr<domain::BoundaryConditions::BoundaryCondition>
      boundary_condition_in_upper_y_;
};

template <size_t Dim>
template <typename BoundaryConditionsBase>
AlignedLattice<Dim>::AlignedLattice(
    typename BlockBounds::type block_bounds,
    typename InitialLevels::type initial_refinement_levels,
    typename InitialGridPoints::type initial_number_of_grid_points,
    typename RefinedLevels::type refined_refinement,
    typename RefinedGridPoints::type refined_grid_points,
    typename BlocksToExclude::type blocks_to_exclude,
    std::variant<std::unique_ptr<BoundaryConditionsBase>,
                 LowerUpperBoundaryCondition<BoundaryConditionsBase>>
        boundary_conditions_in_x,
    std::variant<std::unique_ptr<BoundaryConditionsBase>,
                 LowerUpperBoundaryCondition<BoundaryConditionsBase>>
        boundary_conditions_in_y,
    const Options::Context& context)
    : AlignedLattice(
          std::move(block_bounds), initial_refinement_levels,
          initial_number_of_grid_points, refined_refinement,
          refined_grid_points, blocks_to_exclude,
          std::holds_alternative<std::unique_ptr<BoundaryConditionsBase>>(
              boundary_conditions_in_x)
              ? std::get<std::unique_ptr<BoundaryConditionsBase>>(
                    boundary_conditions_in_x)
                    ->get_clone()
              : std::move(
                    std::get<
                        LowerUpperBoundaryCondition<BoundaryConditionsBase>>(
                        boundary_conditions_in_x)
                        .lower),

          std::holds_alternative<std::unique_ptr<BoundaryConditionsBase>>(
              boundary_conditions_in_x)
              ? std::get<std::unique_ptr<BoundaryConditionsBase>>(
                    boundary_conditions_in_x)
                    ->get_clone()
              : std::move(
                    std::get<
                        LowerUpperBoundaryCondition<BoundaryConditionsBase>>(
                        boundary_conditions_in_x)
                        .upper),
          std::holds_alternative<std::unique_ptr<BoundaryConditionsBase>>(
              boundary_conditions_in_y)
              ? std::get<std::unique_ptr<BoundaryConditionsBase>>(
                    boundary_conditions_in_y)
                    ->get_clone()
              : std::move(
                    std::get<
                        LowerUpperBoundaryCondition<BoundaryConditionsBase>>(
                        boundary_conditions_in_y)
                        .lower),
          std::holds_alternative<std::unique_ptr<BoundaryConditionsBase>>(
              boundary_conditions_in_y)
              ? std::get<std::unique_ptr<BoundaryConditionsBase>>(
                    boundary_conditions_in_y)
                    ->get_clone()
              : std::move(
                    std::get<
                        LowerUpperBoundaryCondition<BoundaryConditionsBase>>(
                        boundary_conditions_in_y)
                        .upper),
          context) {}

}  // namespace creators
}  // namespace domain
