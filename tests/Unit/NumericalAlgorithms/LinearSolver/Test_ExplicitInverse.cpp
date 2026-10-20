// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <algorithm>
#include <blaze/math/DynamicMatrix.h>
#include <blaze/math/DynamicVector.h>
#include <complex>
#include <cstddef>
#include <functional>
#include <utility>

#include "DataStructures/ApplyMatrices.hpp"
#include "DataStructures/ComplexDataVector.hpp"
#include "DataStructures/DataBox/Tag.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Matrix.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Domain/Structure/Direction.hpp"
#include "Domain/Structure/ElementId.hpp"
#include "Framework/TestHelpers.hpp"
#include "Helpers/NumericalAlgorithms/LinearSolver/TestHelpers.hpp"
#include "NumericalAlgorithms/LinearSolver/ExplicitInverse.hpp"
#include "ParallelAlgorithms/LinearSolver/Schwarz/ElementCenteredSubdomainData.hpp"
#include "ParallelAlgorithms/LinearSolver/Schwarz/OverlapHelpers.hpp"
#include "Utilities/MakeWithValue.hpp"
#include "Utilities/TMPL.hpp"

namespace helpers = TestHelpers::LinearSolver;

namespace {
template <typename DataType>
struct ScalarFieldTag : db::SimpleTag {
  using type = Scalar<DataType>;
};

template <typename DataType>
struct VectorFieldTag : db::SimpleTag {
  using type = tnsr::I<DataType, 2>;
};

// Applies a dense matrix to the subdomain data in its iteration order, times a
// `scale` that is passed as operator argument. Provides shift weights that
// differ between components, points and overlaps.
template <typename SubdomainData>
struct OperatorWithShiftWeights {
  blaze::DynamicMatrix<std::complex<double>> matrix;

  void operator()(const gsl::not_null<SubdomainData*> result,
                  const SubdomainData& operand, const double scale) const {
    blaze::DynamicVector<std::complex<double>> contiguous_operand(
        operand.size());
    std::copy(operand.begin(), operand.end(), contiguous_operand.begin());
    const blaze::DynamicVector<std::complex<double>> contiguous_result =
        scale * matrix * contiguous_operand;
    std::copy(contiguous_result.begin(), contiguous_result.end(),
              result->begin());
  }

  void shift_weights(const gsl::not_null<SubdomainData*> weights,
                     const double scale) const {
    const auto set_weights = [&scale](const auto vars, const double offset) {
      auto& field = get<VectorFieldTag<ComplexDataVector>>(*vars);
      for (size_t d = 0; d < 2; ++d) {
        for (size_t i = 0; i < field.get(d).size(); ++i) {
          field.get(d)[i] = scale * (offset + 0.1 * static_cast<double>(i) +
                                     0.01 * static_cast<double>(d));
        }
      }
    };
    set_weights(make_not_null(&weights->element_data), 1.);
    for (auto& [overlap_id, overlap_weights] : weights->overlap_data) {
      set_weights(make_not_null(&overlap_weights),
                  overlap_id.direction().side() == Side::Lower ? 2. : 3.);
    }
  }
};
}  // namespace

namespace LinearSolver::Serial {

SPECTRE_TEST_CASE("Unit.LinearSolver.Serial.ExplicitInverse",
                  "[Unit][NumericalAlgorithms][LinearSolver]") {
  {
    INFO("Solve a simple matrix");
    const blaze::DynamicMatrix<double> matrix{{4., 1.}, {3., 1.}};
    const helpers::ApplyMatrix<double> linear_operator{matrix};
    const blaze::DynamicVector<double> source{1., 2.};
    const blaze::DynamicVector<double> expected_solution{-1., 5.};
    blaze::DynamicVector<double> solution(2);
    const ExplicitInverse<double> solver{"Matrix"};
    const auto has_converged =
        solver.solve(make_not_null(&solution), linear_operator, source);
    REQUIRE(has_converged);
    CHECK_ITERABLE_APPROX(solver.matrix_representation(), blaze::inv(matrix));
    CHECK_ITERABLE_APPROX(solution, expected_solution);
    std::ifstream matrix_file("Matrix.txt");
    std::string matrix_csv((std::istreambuf_iterator<char>(matrix_file)),
                           std::istreambuf_iterator<char>());
    CHECK(matrix_csv == "4 1\n3 1\n");
    {
      INFO("Resetting");
      ExplicitInverse<double> resetting_solver{};
      resetting_solver.solve(make_not_null(&solution), linear_operator, source);
      // Solving a different operator after resetting should work
      resetting_solver.reset();
      const blaze::DynamicMatrix<double> matrix2{{4., 1.}, {1., 3.}};
      const helpers::ApplyMatrix<double> linear_operator2{matrix2};
      const blaze::DynamicVector<double> expected_solution2{0.0909090909090909,
                                                            0.6363636363636364};
      resetting_solver.solve(make_not_null(&solution), linear_operator2,
                             source);
      CHECK_ITERABLE_APPROX(resetting_solver.matrix_representation(),
                            blaze::inv(matrix2));
      CHECK_ITERABLE_APPROX(solution, expected_solution2);
      // Without resetting, the solver should keep applying the cached
      // inverse even when solving a different operator
      solver.solve(make_not_null(&solution), linear_operator2, source);
      // Still the inverse of the operator we solved first
      CHECK_ITERABLE_APPROX(solver.matrix_representation(), blaze::inv(matrix));
      CHECK_ITERABLE_APPROX(solution, expected_solution);
    }
  }
  {
    INFO("Solve a complex matrix");
    const blaze::DynamicMatrix<std::complex<double>> matrix{
        {std::complex<double>(1., 2.), std::complex<double>(2., -1.)},
        {std::complex<double>(3., 4.), std::complex<double>(4., 1.)}};
    const helpers::ApplyMatrix<std::complex<double>> linear_operator{matrix};
    const blaze::DynamicVector<std::complex<double>> source{
        std::complex<double>(1., 1.), std::complex<double>(2., -3.)};
    const blaze::DynamicVector<std::complex<double>> expected_solution{
        std::complex<double>(0.45, -1.4), std::complex<double>(-1.2, 0.15)};
    blaze::DynamicVector<std::complex<double>> solution(2);
    const ExplicitInverse<std::complex<double>> solver{"Matrix"};
    const auto has_converged =
        solver.solve(make_not_null(&solution), linear_operator, source);
    REQUIRE(has_converged);
    CHECK_ITERABLE_APPROX(solver.matrix_representation(), blaze::inv(matrix));
    CHECK_ITERABLE_APPROX(solution, expected_solution);
    std::ifstream matrix_file("Matrix.txt");
    std::string matrix_csv((std::istreambuf_iterator<char>(matrix_file)),
                           std::istreambuf_iterator<char>());
    CHECK(matrix_csv == "(1,2) (2,-1)\n(3,4) (4,1)\n");
    {
      INFO("With complex shift");
      const ExplicitInverse<std::complex<double>> solver_shifted{
          std::nullopt, /* complex_shift */ 0.01};
      const auto has_converged_shifted = solver_shifted.solve(
          make_not_null(&solution), linear_operator, source);
      REQUIRE(has_converged_shifted);
      CHECK_ITERABLE_APPROX(
          solver_shifted.matrix_representation(),
          blaze::inv(matrix -
                     std::complex<double>(0., 0.01) *
                         blaze::IdentityMatrix<std::complex<double>>(2)));
    }
  }
  {
    INFO("Solve a heterogeneous data structure");
    using SubdomainData = ::LinearSolver::Schwarz::ElementCenteredSubdomainData<
        1, tmpl::list<ScalarFieldTag<DataVector>>>;

    const Matrix matrix_element{{4., 1., 1.}, {1., 1., 3.}, {0., 2., 0.}};
    const Matrix matrix_overlap{{4., 1.}, {3., 1.}};
    const ::LinearSolver::Schwarz::OverlapId<1> overlap_id{
        Direction<1>::lower_xi(), ElementId<1>{0}};
    const std::array<std::reference_wrapper<const Matrix>, 1> matrices_element{
        matrix_element};
    const std::array<std::reference_wrapper<const Matrix>, 1> matrices_overlap{
        matrix_overlap};
    const auto linear_operator = [&matrices_element, &matrices_overlap,
                                  &overlap_id](
                                     const gsl::not_null<SubdomainData*> result,
                                     const SubdomainData& operand) {
      apply_matrices(make_not_null(&result->element_data), matrices_element,
                     operand.element_data, Index<1>{3});
      apply_matrices(make_not_null(&result->overlap_data.at(overlap_id)),
                     matrices_overlap, operand.overlap_data.at(overlap_id),
                     Index<1>{2});
    };

    SubdomainData source{3};
    get(get<ScalarFieldTag<DataVector>>(source.element_data)) =
        DataVector{1., 2., 1.};
    source.overlap_data.emplace(overlap_id,
                                typename SubdomainData::OverlapData{2});
    get(get<ScalarFieldTag<DataVector>>(source.overlap_data.at(overlap_id))) =
        DataVector{1., 2.};
    auto expected_solution = make_with_value<SubdomainData>(source, 0.);
    get(get<ScalarFieldTag<DataVector>>(expected_solution.element_data)) =
        DataVector{0., 0.5, 0.5};
    get(get<ScalarFieldTag<DataVector>>(
        expected_solution.overlap_data.at(overlap_id))) = DataVector{-1., 5.};

    const ExplicitInverse<double> solver{};
    auto solution = make_with_value<SubdomainData>(source, 0.);
    solver.solve(make_not_null(&solution), linear_operator, source);
    CHECK(solver.size() == 5);
    Matrix expected_matrix(5, 5, 0.);
    blaze::submatrix(expected_matrix, 0, 0, 3, 3) = matrix_element;
    blaze::submatrix(expected_matrix, 3, 3, 2, 2) = matrix_overlap;
    blaze::invert(expected_matrix);
    CHECK_ITERABLE_APPROX(solver.matrix_representation(), expected_matrix);
    CHECK_VARIABLES_APPROX(solution.element_data,
                           expected_solution.element_data);
    CHECK_VARIABLES_APPROX(solution.overlap_data.at(overlap_id),
                           expected_solution.overlap_data.at(overlap_id));
  }
  {
    INFO("Complex shift weighted by the operator's shift weights");
    using SubdomainData = ::LinearSolver::Schwarz::ElementCenteredSubdomainData<
        1, tmpl::list<VectorFieldTag<ComplexDataVector>>>;
    const ::LinearSolver::Schwarz::OverlapId<1> lower_overlap_id{
        Direction<1>::lower_xi(), ElementId<1>{0}};
    const ::LinearSolver::Schwarz::OverlapId<1> upper_overlap_id{
        Direction<1>::upper_xi(), ElementId<1>{1}};
    SubdomainData source{3};
    // Insert the overlaps out of their sorted order
    source.overlap_data.emplace(upper_overlap_id,
                                typename SubdomainData::OverlapData{1});
    source.overlap_data.emplace(lower_overlap_id,
                                typename SubdomainData::OverlapData{2});
    const size_t size = source.size();
    REQUIRE(size == 12);
    {
      double value = 1.;
      for (auto& entry : source) {
        entry = std::complex<double>(value, -0.5 * value);
        value += 1.;
      }
    }
    blaze::DynamicMatrix<std::complex<double>> matrix(size, size);
    for (size_t i = 0; i < size; ++i) {
      for (size_t j = 0; j < size; ++j) {
        matrix(i, j) = std::complex<double>(
            1. / (1. + static_cast<double>(i + j)),
            i == j ? 2. : 0.1 * static_cast<double>(i) - 0.2);
      }
    }
    const double scale = 1.5;
    const OperatorWithShiftWeights<SubdomainData> linear_operator{matrix};
    const double complex_shift = 0.3;
    auto solution = make_with_value<SubdomainData>(source, 0.);
    {
      INFO("The weights multiply the shift");
      const ExplicitInverse<std::complex<double>> solver{std::nullopt,
                                                         complex_shift};
      solver.solve(make_not_null(&solution), linear_operator, source,
                   std::forward_as_tuple(scale));
      // The diagonal entry i belongs to the i-th entry of the subdomain data
      auto weights = make_with_value<SubdomainData>(source, 0.);
      linear_operator.shift_weights(make_not_null(&weights), scale);
      blaze::DynamicMatrix<std::complex<double>> expected_matrix =
          scale * matrix;
      size_t i = 0;
      for (const auto& weight : weights) {
        expected_matrix(i, i) -=
            std::complex<double>(0., complex_shift) * weight;
        ++i;
      }
      CHECK_ITERABLE_APPROX(solver.matrix_representation(),
                            blaze::inv(expected_matrix));
    }
    {
      INFO("No shift");
      const ExplicitInverse<std::complex<double>> solver{std::nullopt, 0.};
      solver.solve(make_not_null(&solution), linear_operator, source,
                   std::forward_as_tuple(scale));
      CHECK_ITERABLE_APPROX(
          solver.matrix_representation(),
          blaze::inv(
              blaze::DynamicMatrix<std::complex<double>>(scale * matrix)));
    }
  }
}

}  // namespace LinearSolver::Serial
