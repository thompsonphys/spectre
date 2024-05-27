// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <string>

#include "DataStructures/DataBox/PrefixHelpers.hpp"
#include "DataStructures/DataBox/Prefixes.hpp"
#include "DataStructures/Variables.hpp"
#include "Domain/Tags/Faces.hpp"
#include "Elliptic/Actions/InitializeAnalyticSolution.hpp"
#include "Elliptic/Actions/InitializeFields.hpp"
#include "Elliptic/Actions/InitializeFixedSources.hpp"
#include "Elliptic/Actions/RunEventsAndTriggers.hpp"
#include "Elliptic/Amr/Actions.hpp"
#include "Elliptic/BoundaryConditions/Tags/BoundaryFields.hpp"
#include "Elliptic/DiscontinuousGalerkin/Actions/ApplyOperator.hpp"
#include "Elliptic/DiscontinuousGalerkin/Actions/InitializeDomain.hpp"
#include "Elliptic/DiscontinuousGalerkin/SubdomainOperator/InitializeSubdomain.hpp"
#include "Elliptic/DiscontinuousGalerkin/SubdomainOperator/SubdomainOperator.hpp"
#include "Elliptic/DiscontinuousGalerkin/Tags.hpp"
#include "Elliptic/Protocols/FirstOrderSystem.hpp"
#include "Elliptic/SubdomainPreconditioners/MinusLaplacian.hpp"
#include "Elliptic/Systems/GetSourcesComputer.hpp"
#include "Elliptic/Tags.hpp"
#include "IO/Observer/Helpers.hpp"
#include "NumericalAlgorithms/Convergence/Tags.hpp"
#include "Options/String.hpp"
#include "Parallel/PhaseControl/ExecutePhaseChange.hpp"
#include "ParallelAlgorithms/Actions/AddComputeTags.hpp"
#include "ParallelAlgorithms/Actions/InitializeItems.hpp"
#include "ParallelAlgorithms/Actions/RandomizeVariables.hpp"
#include "ParallelAlgorithms/Amr/Actions/Component.hpp"
#include "ParallelAlgorithms/Amr/Actions/Initialize.hpp"
#include "ParallelAlgorithms/Amr/Projectors/DefaultInitialize.hpp"
#include "ParallelAlgorithms/Amr/Projectors/Variables.hpp"
#include "ParallelAlgorithms/Amr/Tags.hpp"
#include "ParallelAlgorithms/LinearSolver/Actions/BuildMatrix.hpp"
#include "ParallelAlgorithms/LinearSolver/Actions/MakeIdentityIfSkipped.hpp"
#include "ParallelAlgorithms/LinearSolver/Gmres/Gmres.hpp"
#include "ParallelAlgorithms/LinearSolver/Multigrid/Actions/RestrictFields.hpp"
#include "ParallelAlgorithms/LinearSolver/Multigrid/Multigrid.hpp"
#include "ParallelAlgorithms/LinearSolver/Schwarz/Actions/CommunicateOverlapFields.hpp"
#include "ParallelAlgorithms/LinearSolver/Schwarz/Actions/ResetSubdomainSolver.hpp"
#include "ParallelAlgorithms/LinearSolver/Schwarz/Schwarz.hpp"
#include "ParallelAlgorithms/LinearSolver/Tags.hpp"
#include "ParallelAlgorithms/NonlinearSolver/NewtonRaphson/NewtonRaphson.hpp"
#include "ParallelAlgorithms/NonlinearSolver/Tags.hpp"
#include "PointwiseFunctions/InitialDataUtilities/AnalyticSolution.hpp"
#include "PointwiseFunctions/InitialDataUtilities/Background.hpp"
#include "PointwiseFunctions/InitialDataUtilities/InitialGuess.hpp"
#include "Utilities/ProtocolHelpers.hpp"
#include "Utilities/TMPL.hpp"

/// Items related to composing nonlinear elliptic solver executables
namespace elliptic {

/// Option tags for nonlinear elliptic solver executables
namespace OptionTags {

struct NonlinearSolverGroup {
  static std::string name() { return "NonlinearSolver"; }
  static constexpr Options::String help = "The iterative nonlinear solver";
};
struct NewtonRaphsonGroup {
  static std::string name() { return "NewtonRaphson"; }
  static constexpr Options::String help =
      "Options for the Newton-Raphson nonlinear solver";
  using group = NonlinearSolverGroup;
};
struct LinearSolverGroup {
  static std::string name() { return "LinearSolver"; }
  static constexpr Options::String help =
      "The iterative Krylov-subspace linear solver";
};
struct GmresGroup {
  static std::string name() { return "Gmres"; }
  static constexpr Options::String help = "Options for the GMRES linear solver";
  using group = LinearSolverGroup;
};
struct SchwarzSmootherGroup {
  static std::string name() { return "SchwarzSmoother"; }
  static constexpr Options::String help = "Options for the Schwarz smoother";
  using group = LinearSolverGroup;
};
struct MultigridGroup {
  static std::string name() { return "Multigrid"; }
  static constexpr Options::String help = "Options for the multigrid";
  using group = LinearSolverGroup;
};

}  // namespace OptionTags

/*!
 * \brief A complete nonlinear elliptic solver stack. Use to compose an
 * executable.
 *
 * The elliptic solver stack is described in detail in \cite Vu2021coj.
 *
 * Uses `Metavariables::volume_dim` and `Metavariables::system`. Also uses
 * `Metavariables` to instantiate parallel components.
 */
template <typename Metavariables>
struct Solver {
  static constexpr size_t volume_dim = Metavariables::volume_dim;
  using system = typename Metavariables::system;
  static_assert(
      tt::assert_conforms_to_v<system, elliptic::protocols::FirstOrderSystem>);
  static constexpr bool is_linear = false;
  //   std::is_same_v<elliptic::get_sources_computer<system, true>,
  //                  typename system::sources_computer>;

  using background_tag =
      elliptic::Tags::Background<elliptic::analytic_data::Background>;
  using initial_guess_tag =
      elliptic::Tags::InitialGuess<elliptic::analytic_data::InitialGuess>;

  /// These are the fields we solve for
  using fields_tag = ::Tags::Variables<typename system::primal_fields>;
  /// These are the fluxes corresponding to the fields, i.e. essentially their
  /// first derivatives. These are background fields for the linearized sources.
  using fluxes_tag = ::Tags::Variables<typename system::primal_fluxes>;
  /// These are the fixed sources, i.e. the RHS of the equations
  using fixed_sources_tag = db::add_tag_prefix<::Tags::FixedSource, fields_tag>;
  using operator_applied_to_fields_tag = tmpl::conditional_t<
      is_linear,
      // This is the linear operator applied to the fields. We only use it to
      // apply the operator to the initial guess, so an optimization would be to
      // re-use the `operator_applied_to_vars_tag` below. This optimization
      // needs a few minor changes to the parallel linear solver algorithm.
      db::add_tag_prefix<LinearSolver::Tags::OperatorAppliedTo, fields_tag>,
      db::add_tag_prefix<NonlinearSolver::Tags::OperatorAppliedTo, fields_tag>>;

  // For AMR iterations
  using amr_iteration_id =
      Convergence::Tags::IterationId<::amr::OptionTags::AmrGroup>;

  /// The nonlinear solver algorithm
  using nonlinear_solver = NonlinearSolver::newton_raphson::NewtonRaphson<
      Metavariables, fields_tag, OptionTags::NewtonRaphsonGroup,
      fixed_sources_tag, LinearSolver::multigrid::Tags::IsFinestGrid>;
  using nonlinear_solver_iteration_id =
      Convergence::Tags::IterationId<typename nonlinear_solver::options_group>;

  /// The linear solver algorithm. We use GMRES since the operator is not
  /// necessarily symmetric. Using CG here would be an optimization for
  /// symmetric problems.
  using linear_solver = LinearSolver::gmres::Gmres<
      Metavariables,
      tmpl::conditional_t<is_linear, fields_tag,
                          typename nonlinear_solver::linear_solver_fields_tag>,
      OptionTags::GmresGroup, true,
      tmpl::conditional_t<is_linear, fixed_sources_tag,
                          typename nonlinear_solver::linear_solver_source_tag>,
      LinearSolver::multigrid::Tags::IsFinestGrid>;
  using linear_solver_iteration_id =
      Convergence::Tags::IterationId<typename linear_solver::options_group>;

  /// Precondition each linear solver iteration with a multigrid V-cycle
  using multigrid = LinearSolver::multigrid::Multigrid<
      Metavariables, typename linear_solver::operand_tag,
      OptionTags::MultigridGroup, elliptic::dg::Tags::Massive,
      typename linear_solver::preconditioner_source_tag>;

  /// Smooth each multigrid level with a number of Schwarz smoothing steps
  using subdomain_operator =
      elliptic::dg::subdomain_operator::SubdomainOperator<
          system, OptionTags::SchwarzSmootherGroup>;
  using subdomain_preconditioners = tmpl::list<
      elliptic::subdomain_preconditioners::Registrars::MinusLaplacian<
          volume_dim, OptionTags::SchwarzSmootherGroup>>;
  using schwarz_smoother = LinearSolver::Schwarz::Schwarz<
      typename multigrid::smooth_fields_tag, OptionTags::SchwarzSmootherGroup,
      subdomain_operator, subdomain_preconditioners,
      typename multigrid::smooth_source_tag,
      LinearSolver::multigrid::Tags::MultigridLevel>;

  /// For the GMRES linear solver we need to apply the DG operator to its
  /// internal "operand" in every iteration of the algorithm.
  using vars_tag = typename linear_solver::operand_tag;
  using operator_applied_to_vars_tag =
      db::add_tag_prefix<LinearSolver::Tags::OperatorAppliedTo, vars_tag>;
  /// The correction fluxes can be stored in an arbitrary tag. We don't need to
  /// access them anywhere, they're just a memory buffer for the linearized
  /// operator.
  using fluxes_vars_tag =
      db::add_tag_prefix<NonlinearSolver::Tags::Correction, fluxes_tag>;

  /// Fields that may be observed to monitor the state of the solver
  using observe_fields = tmpl::conditional_t<
      is_linear,
      tmpl::append<
          typename fixed_sources_tag::tags_list,
          typename db::add_tag_prefix<LinearSolver::Tags::Operand,
                                      fields_tag>::tags_list,
          typename db::add_tag_prefix<LinearSolver::Tags::OperatorAppliedTo,
                                      fields_tag>::tags_list>,
      tmpl::append<
          typename fixed_sources_tag::tags_list,
          typename nonlinear_solver::linear_solver_fields_tag::tags_list,
          typename nonlinear_solver::linear_solver_source_tag::tags_list>>;

  /// Collect all reduction tags for observers
  using observed_reduction_data_tags =
      observers::collect_reduction_data_tags<tmpl::flatten<tmpl::list<
          tmpl::conditional_t<is_linear, tmpl::list<>, nonlinear_solver>,
          linear_solver, multigrid, schwarz_smoother>>>;

  template <bool Linearized>
  using dg_operator = elliptic::dg::Actions::DgOperator<
      system, Linearized,
      tmpl::conditional_t<Linearized, linear_solver_iteration_id,
                          nonlinear_solver_iteration_id>,
      tmpl::conditional_t<Linearized, vars_tag, fields_tag>,
      tmpl::conditional_t<Linearized, fluxes_vars_tag, fluxes_tag>,
      tmpl::conditional_t<Linearized, operator_applied_to_vars_tag,
                          operator_applied_to_fields_tag>>;

  using build_matrix = LinearSolver::Actions::BuildMatrix<
      linear_solver_iteration_id, vars_tag, operator_applied_to_vars_tag,
      domain::Tags::Coordinates<volume_dim, Frame::Inertial>,
      LinearSolver::multigrid::Tags::IsFinestGrid>;

  using build_matrix_actions = typename build_matrix::template actions<
      typename dg_operator<true>::apply_actions>;

  // For labeling the yaml option for RandomizeVariables
  struct RandomizeInitialGuess {};

  using init_analytic_solution_action =
      elliptic::Actions::InitializeOptionalAnalyticSolution<
          volume_dim, background_tag,
          tmpl::append<typename system::primal_fields,
                       typename system::primal_fluxes>,
          elliptic::analytic_data::AnalyticSolution>;

  using initialization_actions = tmpl::list<
      elliptic::dg::Actions::InitializeDomain<volume_dim>,
      tmpl::conditional_t<is_linear, tmpl::list<>,
                          typename nonlinear_solver::initialize_element>,
      typename linear_solver::initialize_element,
      typename multigrid::initialize_element,
      typename schwarz_smoother::initialize_element,
      Initialization::Actions::InitializeItems<
          ::amr::Initialization::Initialize<volume_dim>,
          elliptic::amr::Actions::Initialize>,
      elliptic::Actions::InitializeFields<system, initial_guess_tag>,
      ::Actions::RandomizeVariables<fields_tag, RandomizeInitialGuess>,
      elliptic::Actions::InitializeFixedSources<system, background_tag>,
      init_analytic_solution_action,
      elliptic::dg::Actions::initialize_operator<system, background_tag>,
      tmpl::conditional_t<
          is_linear, tmpl::list<>,
          ::Initialization::Actions::AddComputeTags<tmpl::list<
              // For linearized boundary conditions
              elliptic::Tags::BoundaryFieldsCompute<volume_dim, fields_tag>,
              elliptic::Tags::BoundaryFluxesCompute<volume_dim, fields_tag,
                                                    fluxes_tag>>>>>;

  using register_actions = tmpl::flatten<tmpl::list<
      tmpl::conditional_t<is_linear, typename build_matrix::register_actions,
                          typename nonlinear_solver::register_element>,
      typename multigrid::register_element,
      typename schwarz_smoother::register_element>>;

  template <typename Label>
  using smooth_actions = typename schwarz_smoother::template solve<
      typename dg_operator<true>::apply_actions, Label>;

  // These tags are communicated on subdomain overlaps to initialize the
  // subdomain geometry. AMR updates these tags, so we have to communicate them
  // after each AMR step.
  using subdomain_init_tags =
      tmpl::list<domain::Tags::Mesh<volume_dim>,
                 domain::Tags::Element<volume_dim>,
                 domain::Tags::NeighborMesh<volume_dim>>;

  /// This data needs to be communicated on subdomain overlap regions
  using communicated_overlap_tags = tmpl::flatten<tmpl::list<
      // For linearized sources
      fields_tag, fluxes_tag,
      // For linearized boundary conditions
      domain::make_faces_tags<volume_dim, typename system::primal_fields>,
      domain::make_faces_tags<
          volume_dim, db::wrap_tags_in<::Tags::NormalDotFlux,
                                       typename system::primal_fields>>>>;

  using init_subdomain_action =
      elliptic::dg::subdomain_operator::Actions::InitializeSubdomain<
          system, background_tag, typename schwarz_smoother::options_group,
          false>;

  template <typename StepActions>
  using linear_solve_actions = typename linear_solver::template solve<
      tmpl::list<
          // Multigrid preconditioning
          typename multigrid::template solve<
              typename dg_operator<true>::apply_actions,
              // Schwarz smoothing on each multigrid level
              smooth_actions<LinearSolver::multigrid::VcycleDownLabel>,
              smooth_actions<LinearSolver::multigrid::VcycleUpLabel>>,
          // Support disabling the preconditioner
          ::LinearSolver::Actions::make_identity_if_skipped<
              multigrid, typename dg_operator<true>::apply_actions>>,
      StepActions>;

  template <typename StepActions>
  using nonlinear_solve_actions = typename nonlinear_solver::template solve<
      typename dg_operator<false>::apply_actions,
      tmpl::list<
          // Transfer data down the multigrid hierachy
          LinearSolver::multigrid::Actions::ReceiveFieldsFromFinerGrid<
              volume_dim, tmpl::list<fields_tag, fluxes_tag>,
              typename multigrid::options_group>,
          LinearSolver::multigrid::Actions::SendFieldsToCoarserGrid<
              tmpl::list<fields_tag, fluxes_tag>,
              typename multigrid::options_group, void>,
          // Communicate data on subdomain overlap regions
          LinearSolver::Schwarz::Actions::SendOverlapFields<
              communicated_overlap_tags,
              typename schwarz_smoother::options_group, false,
              nonlinear_solver_iteration_id>,
          LinearSolver::Schwarz::Actions::ReceiveOverlapFields<
              volume_dim, communicated_overlap_tags,
              typename schwarz_smoother::options_group, false,
              nonlinear_solver_iteration_id>,
          // Reset Schwarz subdomain solver
          LinearSolver::Schwarz::Actions::ResetSubdomainSolver<
              typename schwarz_smoother::options_group>,
          // Linear solve for correction
          linear_solve_actions<tmpl::list<>>>,
      StepActions>;

  template <typename StepActions>
  using solve_actions = tmpl::flatten<tmpl::list<
      // Observe initial state
      elliptic::Actions::RunEventsAndTriggers<amr_iteration_id>,
      // Stop AMR iterations if complete
      elliptic::amr::Actions::StopAmr,
      // Run AMR or other phase changes if requested
      PhaseControl::Actions::ExecutePhaseChange, StepActions,
      // Increment AMR iteration ID
      elliptic::amr::Actions::IncrementIterationId,
      // Communicate subdomain geometry and initialize subdomain to account for
      // domain changes
      LinearSolver::Schwarz::Actions::SendOverlapFields<
          subdomain_init_tags, typename schwarz_smoother::options_group, false,
          amr_iteration_id>,
      LinearSolver::Schwarz::Actions::ReceiveOverlapFields<
          volume_dim, subdomain_init_tags,
          typename schwarz_smoother::options_group, false, amr_iteration_id>,
      init_subdomain_action,
      // Linear or nonlinear solve
      tmpl::conditional_t<
          is_linear,
          // Linear solve
          tmpl::list<
              // Apply the DG operator to the initial guess
              typename elliptic::dg::Actions::DgOperator<
                  system, true, linear_solver_iteration_id, fields_tag,
                  fluxes_vars_tag, operator_applied_to_fields_tag, vars_tag,
                  fluxes_vars_tag>::apply_actions,
              // Modify fixed sources with boundary conditions
              elliptic::dg::Actions::
                  ImposeInhomogeneousBoundaryConditionsOnSource<
                      system, fixed_sources_tag>,
              // Krylov solve
              linear_solve_actions<tmpl::list<>>>,
          // Nonlinear solve
          nonlinear_solve_actions<tmpl::list<>>>>>;

  template <typename Tag>
  using overlaps_tag =
      LinearSolver::Schwarz::Tags::Overlaps<Tag, volume_dim,
                                            OptionTags::SchwarzSmootherGroup>;

  using amr_projectors = tmpl::flatten<tmpl::list<
      elliptic::dg::ProjectGeometry<volume_dim>,
      tmpl::conditional_t<is_linear, tmpl::list<>,
                          typename nonlinear_solver::amr_projectors>,
      typename linear_solver::amr_projectors,
      typename multigrid::amr_projectors,
      typename schwarz_smoother::amr_projectors,
      ::amr::projectors::DefaultInitialize<tmpl::append<
          tmpl::list<domain::Tags::InitialExtents<volume_dim>,
                     domain::Tags::InitialRefinementLevels<volume_dim>>,
          // Tags communicated on subdomain overlaps. No need to project
          // these during AMR because they will be communicated.
          db::wrap_tags_in<
              overlaps_tag,
              tmpl::append<subdomain_init_tags,
                           tmpl::conditional_t<is_linear, tmpl::list<>,
                                               communicated_overlap_tags>>>,
          // Tags initialized on subdomains. No need to project these during
          // AMR because they will get re-initialized after communication.
          typename init_subdomain_action::simple_tags>>,
      ::amr::projectors::ProjectVariables<volume_dim, fields_tag>,
      elliptic::Actions::InitializeFixedSources<system, background_tag>,
      init_analytic_solution_action,
      elliptic::dg::Actions::amr_projectors<system, background_tag>,
      typename dg_operator<true>::amr_projectors,
      tmpl::conditional_t<is_linear, typename build_matrix::amr_projectors,
                          typename dg_operator<false>::amr_projectors>,
      elliptic::amr::Actions::Initialize>>;

  using component_list = tmpl::flatten<
      tmpl::list<tmpl::conditional_t<is_linear, tmpl::list<>,
                                     typename nonlinear_solver::component_list>,
                 typename linear_solver::component_list,
                 typename multigrid::component_list,
                 typename schwarz_smoother::component_list,
                 ::amr::Component<Metavariables>>>;
};

}  // namespace elliptic
