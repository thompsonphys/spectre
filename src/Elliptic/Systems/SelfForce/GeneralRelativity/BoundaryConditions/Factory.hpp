// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include "Elliptic/Systems/SelfForce/GeneralRelativity/BoundaryConditions/Angular.hpp"
#include "Elliptic/Systems/SelfForce/GeneralRelativity/BoundaryConditions/Sommerfeld.hpp"
#include "Utilities/TMPL.hpp"

namespace GrSelfForce::BoundaryConditions {

using standard_boundary_conditions = tmpl::list<Angular, Sommerfeld>;

}
