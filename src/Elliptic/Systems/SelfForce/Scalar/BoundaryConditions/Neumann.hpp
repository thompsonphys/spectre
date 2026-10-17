// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <memory>
#include <vector>

#include "DataStructures/ComplexDataVector.hpp"
#include "DataStructures/Tensor/TypeAliases.hpp"
#include "Elliptic/BoundaryConditions/BoundaryCondition.hpp"
#include "Elliptic/BoundaryConditions/BoundaryConditionType.hpp"
#include "Options/String.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/TMPL.hpp"

namespace ScalarSelfForce::BoundaryConditions {

/*!
 * \brief Imposes a zero normal derivative (homogeneous Neumann condition)
 * \f$n_i F^i = 0\f$, where \f$F^i = \alpha^i \partial_i \Phi\f$ is the flux of
 * the scalar self-force system.
 *
 * \details Used at the equator of a domain in \f$\cos\theta \in [0, 1]\f$ for
 * an equatorial orbit. Every \f$(m, n)\f$-mode of the field is symmetric under
 * \f$\theta \to \pi - \theta\f$, so its derivative in \f$\cos\theta\f$
 * vanishes at the equator, where \f$\alpha^{\cos\theta}\f$ is nonzero.
 */
class Neumann : public elliptic::BoundaryConditions::BoundaryCondition<2> {
 private:
  using Base = elliptic::BoundaryConditions::BoundaryCondition<2>;

 public:
  static constexpr Options::String help =
      "Zero normal derivative (homogeneous Neumann condition), for example at "
      "the equator of a domain in cos(theta) with equatorial symmetry.";
  using options = tmpl::list<>;

  Neumann() = default;
  Neumann(const Neumann&) = default;
  Neumann& operator=(const Neumann&) = default;
  Neumann(Neumann&&) = default;
  Neumann& operator=(Neumann&&) = default;
  ~Neumann() override = default;

  /// \cond
  explicit Neumann(CkMigrateMessage* m);
  using PUP::able::register_constructor;
  WRAPPED_PUPable_decl_template(Neumann);
  /// \endcond

  std::unique_ptr<domain::BoundaryConditions::BoundaryCondition> get_clone()
      const override;

  std::vector<elliptic::BoundaryConditionType> boundary_condition_types()
      const override {
    return {elliptic::BoundaryConditionType::Neumann};
  }

  using argument_tags = tmpl::list<>;
  using volume_tags = tmpl::list<>;

  void apply(gsl::not_null<Scalar<ComplexDataVector>*> field,
             gsl::not_null<Scalar<ComplexDataVector>*> n_dot_field_gradient,
             const tnsr::i<ComplexDataVector, 2>& deriv_field) const;

  using argument_tags_linearized = tmpl::list<>;
  using volume_tags_linearized = tmpl::list<>;

  void apply_linearized(
      gsl::not_null<Scalar<ComplexDataVector>*> field_correction,
      gsl::not_null<Scalar<ComplexDataVector>*> n_dot_field_gradient_correction,
      const tnsr::i<ComplexDataVector, 2>& deriv_field_correction) const;

 private:
  friend bool operator==(const Neumann& lhs, const Neumann& rhs);
};

bool operator!=(const Neumann& lhs, const Neumann& rhs);

}  // namespace ScalarSelfForce::BoundaryConditions
