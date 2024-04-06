// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>
#include <cstddef>
#include <limits>
#include <pup.h>
#include <vector>

#include "DataStructures/DataBox/Prefixes.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Elliptic/Systems/SelfForce/Scalar/Tags.hpp"
#include "NumericalAlgorithms/LinearOperators/PartialDerivatives.hpp"
#include "NumericalAlgorithms/Spectral/Mesh.hpp"
#include "Options/String.hpp"
#include "PointwiseFunctions/InitialDataUtilities/Background.hpp"
#include "PointwiseFunctions/InitialDataUtilities/InitialGuess.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/TMPL.hpp"
#include "Utilities/TaggedTuple.hpp"

namespace ScalarSelfForce::AnalyticData {

class CircularOrbit : public elliptic::analytic_data::Background,
                      public elliptic::analytic_data::InitialGuess {
 public:
  struct ScalarCharge {
    static constexpr Options::String help = "Value of the scalar charge 'q'";
    using type = double;
  };
  struct BlackHoleMass {
    static constexpr Options::String help =
        "Kerr mass parameter 'M' of the black hole";
    using type = double;
  };
  struct BlackHoleSpin {
    static constexpr Options::String help =
        "Kerr dimensionless spin parameter 'chi' of the black hole";
    using type = double;
  };
  struct OrbitalRadius {
    static constexpr Options::String help =
        "Radius 'r_0' of the circular orbit";
    using type = double;
  };
  struct MModeNumber {
    static constexpr Options::String help =
        "Mode number 'm' of the scalar field";
    using type = int;
  };
  using options = tmpl::list<ScalarCharge, BlackHoleMass, BlackHoleSpin,
                             OrbitalRadius, MModeNumber>;
  static constexpr Options::String help =
      "Quasicircular orbit of a scalar point charge in Kerr spacetime";

  CircularOrbit() = default;
  CircularOrbit(const CircularOrbit&) = default;
  CircularOrbit& operator=(const CircularOrbit&) = default;
  CircularOrbit(CircularOrbit&&) = default;
  CircularOrbit& operator=(CircularOrbit&&) = default;
  ~CircularOrbit() = default;

  CircularOrbit(double scalar_charge, double black_hole_mass,
                double black_hole_spin, double orbital_radius,
                int m_mode_number);

  explicit CircularOrbit(CkMigrateMessage* m)
      : elliptic::analytic_data::Background(m),
        elliptic::analytic_data::InitialGuess(m) {}
  using PUP::able::register_constructor;
  WRAPPED_PUPable_decl_template(CircularOrbit);

  tnsr::I<double, 2> puncture_position() const;

  // Background
  tuples::TaggedTuple<Tags::Alpha, Tags::Beta, Tags::Gamma> variables(
      const tnsr::I<DataVector, 2>& x,
      tmpl::list<Tags::Alpha, Tags::Beta, Tags::Gamma> /*meta*/) const;

  // Initial guess
  tuples::TaggedTuple<Tags::MMode> variables(
      const tnsr::I<DataVector, 2>& x, tmpl::list<Tags::MMode> /*meta*/) const;

  // Fixed sources
  tuples::TaggedTuple<
      ::Tags::FixedSource<Tags::MMode>, Tags::SingularField,
      ::Tags::deriv<Tags::SingularField, tmpl::size_t<2>, Frame::Inertial>,
      Tags::BoyerLindquistRadius>
  variables(const tnsr::I<DataVector, 2>& x,
            tmpl::list<::Tags::FixedSource<Tags::MMode>, Tags::SingularField,
                       ::Tags::deriv<Tags::SingularField, tmpl::size_t<2>,
                                     Frame::Inertial>,
                       Tags::BoyerLindquistRadius> /*meta*/) const;

  template <typename... RequestedTags>
  tuples::TaggedTuple<RequestedTags...> variables(
      const tnsr::I<DataVector, 2>& x, const Mesh<2>& /*mesh*/,
      const InverseJacobian<DataVector, 2, Frame::ElementLogical,
                            Frame::Inertial>& /*inv_jacobian*/,
      tmpl::list<RequestedTags...> /*meta*/) const {
    return variables(x, tmpl::list<RequestedTags...>{});
  }

  // NOLINTNEXTLINE
  void pup(PUP::er& p) override;

 private:
  friend bool operator==(const CircularOrbit& lhs, const CircularOrbit& rhs);

  double scalar_charge_{std::numeric_limits<double>::signaling_NaN()};
  double black_hole_mass_{std::numeric_limits<double>::signaling_NaN()};
  double black_hole_spin_{std::numeric_limits<double>::signaling_NaN()};
  double orbital_radius_{std::numeric_limits<double>::signaling_NaN()};
  int m_mode_number_{};
};

bool operator!=(const CircularOrbit& lhs, const CircularOrbit& rhs);

}  // namespace ScalarSelfForce::AnalyticData
