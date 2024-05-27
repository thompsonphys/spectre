// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <string>
#include <vector>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/TypeAliases.hpp"
#include "Domain/Tags.hpp"
#include "Elliptic/BoundaryConditions/BoundaryCondition.hpp"
#include "Elliptic/BoundaryConditions/BoundaryConditionType.hpp"
#include "Options/String.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/TMPL.hpp"

namespace ScalarSelfForce::BoundaryConditions {

class Angular : public elliptic::BoundaryConditions::BoundaryCondition<2> {
 private:
  using Base = elliptic::BoundaryConditions::BoundaryCondition<2>;

 public:
  struct MModeNumber {
    static constexpr Options::String help =
        "Mode number 'm' of the scalar field";
    using type = int;
  };

  static constexpr Options::String help = "Angular boundary condition";
  using options = tmpl::list<MModeNumber>;

  Angular() = default;
  Angular(const Angular&) = default;
  Angular& operator=(const Angular&) = default;
  Angular(Angular&&) = default;
  Angular& operator=(Angular&&) = default;
  ~Angular() = default;

  explicit Angular(int m_mode_number);

  /// \cond
  explicit Angular(CkMigrateMessage* m) : Base(m) {}
  using PUP::able::register_constructor;
  WRAPPED_PUPable_decl_template(Angular);
  /// \endcond

  std::unique_ptr<domain::BoundaryConditions::BoundaryCondition> get_clone()
      const override {
    return std::make_unique<Angular>(*this);
  }

  std::vector<elliptic::BoundaryConditionType> boundary_condition_types()
      const override {
    return {m_mode_number_ == 0 ? elliptic::BoundaryConditionType::Neumann
                                : elliptic::BoundaryConditionType::Dirichlet};
  }

  using argument_tags = tmpl::list<>;
  using volume_tags = tmpl::list<>;

  void apply(gsl::not_null<tnsr::I<DataVector, 2>*> field,
             gsl::not_null<tnsr::I<DataVector, 2>*> n_dot_field_gradient) const;

  using argument_tags_linearized = tmpl::list<>;
  using volume_tags_linearized = tmpl::list<>;

  void apply_linearized(gsl::not_null<tnsr::I<DataVector, 2>*> field_correction,
                        gsl::not_null<tnsr::I<DataVector, 2>*>
                            n_dot_field_gradient_correction) const;

 private:
  friend bool operator==(const Angular& lhs, const Angular& rhs);

  int m_mode_number_{};
};

bool operator!=(const Angular& lhs, const Angular& rhs);

}  // namespace ScalarSelfForce::BoundaryConditions
