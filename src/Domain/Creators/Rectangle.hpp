// Distributed under the MIT License.
// See LICENSE.txt for details.

/// \file
/// Defines class Rectangle.

#pragma once

#include <array>
#include <cstddef>
#include <memory>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "Domain/BoundaryConditions/BoundaryCondition.hpp"
#include "Domain/BoundaryConditions/GetBoundaryConditionsBase.hpp"
#include "Domain/Creators/DomainCreator.hpp"  // IWYU pragma: keep
#include "Domain/Creators/TimeDependence/TimeDependence.hpp"
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
}  // namespace CoordinateMaps

template <typename SourceFrame, typename TargetFrame, typename... Maps>
class CoordinateMap;
}  // namespace domain
/// \endcond

namespace domain {
namespace creators {
/// Create a 2D Domain consisting of a single Block.
class Rectangle : public DomainCreator<2> {
 public:
  using maps_list = tmpl::list<domain::CoordinateMap<
      Frame::BlockLogical, Frame::Inertial,
      CoordinateMaps::ProductOf2Maps<CoordinateMaps::Affine,
                                     CoordinateMaps::Affine>>>;

  struct LowerBound {
    using type = std::array<double, 2>;
    static constexpr Options::String help = {
        "Sequence of [x,y] for lower bounds."};
  };

  struct UpperBound {
    using type = std::array<double, 2>;
    static constexpr Options::String help = {
        "Sequence of [x,y] for upper bounds."};
  };
  struct IsPeriodicIn {
    using type = std::array<bool, 2>;
    static constexpr Options::String help = {
        "Sequence for [x,y], true if periodic."};
  };

  struct InitialRefinement {
    using type = std::array<size_t, 2>;
    static constexpr Options::String help = {
        "Initial refinement level in [x,y]."};
  };

  struct InitialGridPoints {
    using type = std::array<size_t, 2>;
    static constexpr Options::String help = {
        "Initial number of grid points in [x,y]."};
  };
  struct TimeDependence {
    using type =
        std::unique_ptr<domain::creators::time_dependence::TimeDependence<2>>;
    static constexpr Options::String help = {
        "The time dependence of the moving mesh domain."};
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
        : lower(std::move(lower_bc)), upper(std::move(upper_bc)){};
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
      tmpl::list<LowerBound, UpperBound, InitialRefinement, InitialGridPoints>;

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
          options_periodic>,
      tmpl::list<TimeDependence>>;

  static constexpr Options::String help{"Creates a 2D rectangle."};

  Rectangle(
      std::array<double, 2> lower_xy, std::array<double, 2> upper_xy,
      std::array<size_t, 2> initial_refinement_level_xy,
      std::array<size_t, 2> initial_number_of_grid_points_in_xy,
      std::array<bool, 2> is_periodic_in_xy,
      std::unique_ptr<domain::creators::time_dependence::TimeDependence<2>>
          time_dependence = nullptr);

  Rectangle(
      std::array<double, 2> lower_xy, std::array<double, 2> upper_xy,
      std::array<size_t, 2> initial_refinement_level_xy,
      std::array<size_t, 2> initial_number_of_grid_points_in_xy,
      std::unique_ptr<domain::BoundaryConditions::BoundaryCondition>
          boundary_condition_in_lower_x,
      std::unique_ptr<domain::BoundaryConditions::BoundaryCondition>
          boundary_condition_in_upper_x,
      std::unique_ptr<domain::BoundaryConditions::BoundaryCondition>
          boundary_condition_in_lower_y,
      std::unique_ptr<domain::BoundaryConditions::BoundaryCondition>
          boundary_condition_in_upper_y,
      std::unique_ptr<domain::creators::time_dependence::TimeDependence<2>>
          time_dependence,
      const Options::Context& context = {});

  template <typename BoundaryConditionsBase>
  Rectangle(
      std::array<double, 2> lower_xy, std::array<double, 2> upper_xy,
      std::array<size_t, 2> initial_refinement_level_xy,
      std::array<size_t, 2> initial_number_of_grid_points_in_xy,
      std::variant<std::unique_ptr<BoundaryConditionsBase>,
                   LowerUpperBoundaryCondition<BoundaryConditionsBase>>
          boundary_conditions_in_x,
      std::variant<std::unique_ptr<BoundaryConditionsBase>,
                   LowerUpperBoundaryCondition<BoundaryConditionsBase>>
          boundary_conditions_in_y,
      std::unique_ptr<domain::creators::time_dependence::TimeDependence<2>>
          time_dependence,
      const Options::Context& context = {});

  Rectangle() = default;
  Rectangle(const Rectangle&) = delete;
  Rectangle(Rectangle&&) = default;
  Rectangle& operator=(const Rectangle&) = delete;
  Rectangle& operator=(Rectangle&&) = default;
  ~Rectangle() override = default;

  Domain<2> create_domain() const override;

  std::vector<DirectionMap<
      2, std::unique_ptr<domain::BoundaryConditions::BoundaryCondition>>>
  external_boundary_conditions() const override;

  std::vector<std::array<size_t, 2>> initial_extents() const override;

  std::vector<std::array<size_t, 2>> initial_refinement_levels() const override;

  auto functions_of_time(const std::unordered_map<std::string, double>&
                             initial_expiration_times = {}) const
      -> std::unordered_map<
          std::string,
          std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>> override;

  std::vector<std::string> block_names() const override { return block_names_; }

 private:
  typename LowerBound::type lower_xy_{};
  typename UpperBound::type upper_xy_{};
  typename IsPeriodicIn::type is_periodic_in_xy_{};
  typename InitialRefinement::type initial_refinement_level_xy_{};
  typename InitialGridPoints::type initial_number_of_grid_points_in_xy_{};
  std::unique_ptr<domain::creators::time_dependence::TimeDependence<2>>
      time_dependence_{nullptr};
  std::unique_ptr<domain::BoundaryConditions::BoundaryCondition>
      boundary_condition_in_lower_x_;
  std::unique_ptr<domain::BoundaryConditions::BoundaryCondition>
      boundary_condition_in_upper_x_;
  std::unique_ptr<domain::BoundaryConditions::BoundaryCondition>
      boundary_condition_in_lower_y_;
  std::unique_ptr<domain::BoundaryConditions::BoundaryCondition>
      boundary_condition_in_upper_y_;
  inline static const std::vector<std::string> block_names_{"Rectangle"};
};

template <typename BoundaryConditionsBase>
Rectangle::Rectangle(
    std::array<double, 2> lower_xy, std::array<double, 2> upper_xy,
    std::array<size_t, 2> initial_refinement_level_xy,
    std::array<size_t, 2> initial_number_of_grid_points_in_xy,
    std::variant<std::unique_ptr<BoundaryConditionsBase>,
                 LowerUpperBoundaryCondition<BoundaryConditionsBase>>
        boundary_conditions_in_x,
    std::variant<std::unique_ptr<BoundaryConditionsBase>,
                 LowerUpperBoundaryCondition<BoundaryConditionsBase>>
        boundary_conditions_in_y,
    std::unique_ptr<domain::creators::time_dependence::TimeDependence<2>>
        time_dependence,
    const Options::Context& context)
    : Rectangle(
          lower_xy, upper_xy, initial_refinement_level_xy,
          initial_number_of_grid_points_in_xy,
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
          std::move(time_dependence), context) {}
}  // namespace creators
}  // namespace domain
