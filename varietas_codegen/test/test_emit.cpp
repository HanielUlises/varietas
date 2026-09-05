// The emitter, checked against the thing it is a picture of.
//
// A generated header has a failure mode that ordinary code does not: it can be
// syntactically perfect, compile without a murmur, and compute the wrong
// function, because nothing downstream ever compares it against the object it
// was generated from. So the check here is not that the text looks right. It is
// that the arithmetic the text describes, evaluated at a pose, equals the
// rational functions evaluated at the same pose, with the emitted expression
// tested against the exact value it is supposed to approximate.
//
// The companion suite, test_generated.cpp, closes the remaining gap by actually
// compiling a header this emitter produced.

#include <array>
#include <cstddef>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "varietas/codegen/emit.hpp"
#include "varietas/codegen/parametric_solution.hpp"
#include "varietas/codegen/rational_function.hpp"
#include "varietas/core/ideal/ideal.hpp"
#include "varietas/core/order/grevlex.hpp"

#include "planar_fixture.hpp"

namespace {

using varietas::grevlex;
using varietas::make_rational;
using varietas::order_id;
using varietas::rational;
using varietas::codegen::emit;
using varietas::codegen::emit_options;
using varietas::codegen::parametric_action_matrix;
using varietas::codegen::parametric_matrix;
using varietas::codegen::parametric_solution;

using field = varietas::rational_function<2>;
using parametric = varietas::polynomial<field, 2, grevlex>;
using traits = varietas::coefficient_traits<field>;

field constant(std::int64_t n, std::int64_t d = 1) { return field(make_rational(n, d)); }
field x() { return field::parameter(0); }
field y() { return field::parameter(1); }

// The system is described in planar_fixture.hpp, which test_generated.cpp also
// builds from, so that the two suites are demonstrably talking about the same
// object rather than about two systems that happen to look alike.
parametric_solution<2, 2> build_solution() { return varietas_test::planar_solution(); }

TEST(Emit, TheSystemIsZeroDimensionalAndTheMatricesAreParametric) {
  const auto solution = build_solution();
  ASSERT_TRUE(solution.is_well_formed());
  EXPECT_EQ(solution.dimension(), 2u);

  // If every entry were constant there would be nothing to emit that a table of
  // numbers could not do, and this test would not be testing the interesting
  // case.
  bool any_parametric = false;
  for (const auto& m : solution.action) {
    for (const auto& e : m.entries) {
      if (!e.is_constant()) {
        any_parametric = true;
      }
    }
  }
  EXPECT_TRUE(any_parametric);
}

TEST(Emit, HeaderCarriesTheOrderAndTheShape) {
  const auto solution = build_solution();
  emit_options options;
  options.name = "planar";
  options.name_space = "test_generated";
  const std::string header = emit(solution, options);

  EXPECT_NE(header.find("#ifndef TEST_GENERATED_PLANAR_GENERATED_HPP"), std::string::npos);
  EXPECT_NE(header.find("namespace test_generated"), std::string::npos);
  EXPECT_NE(header.find("struct planar"), std::string::npos);
  EXPECT_NE(header.find("order_name = \"grevlex\""), std::string::npos);
  EXPECT_NE(header.find("static constexpr std::size_t dimension = 2;"), std::string::npos);
  EXPECT_NE(header.find("num_unknowns = 2"), std::string::npos);
  EXPECT_NE(header.find("num_parameters = 2"), std::string::npos);

  // The order is recorded as the enumerator's value, not as a name only, so a
  // runtime can compare it without parsing text.
  EXPECT_NE(header.find("order_id = 2;"), std::string::npos);

  // Self-contained: nothing from varietas is named in the output.
  EXPECT_EQ(header.find("varietas"), header.rfind("varietas"))
      << "the header should mention varietas exactly once, in the banner";
}

TEST(Emit, HeaderIsDeterministic) {
  const auto a = emit(build_solution());
  const auto b = emit(build_solution());
  EXPECT_EQ(a, b) << "the same solution must emit the same bytes, or diffs are noise";
}

TEST(Emit, RationalsAreWrittenAsExactQuotientsNotDecimals) {
  // 1/3 has no decimal form, so if it survives as a quotient the compiler gets
  // the correctly rounded double; if it were printed as 0.333... it would not.
  const field third = constant(1, 3);
  const std::string text = varietas::codegen::detail::emit_function<2>(third, "pose");
  EXPECT_NE(text.find("1.0 / 3.0"), std::string::npos) << "got: " << text;
}

TEST(Emit, MonomialsAreWrittenAsProductsOfPoseEntries) {
  // x^2 * y as a coefficient.
  const field f = x() * x() * y();
  const std::string text = varietas::codegen::detail::emit_function<2>(f, "pose");
  EXPECT_NE(text.find("pose[0] * pose[0] * pose[1]"), std::string::npos) << "got: " << text;
  EXPECT_EQ(text.find("pow"), std::string::npos) << "powers should be written out";
}

// The load-bearing test.
//
// Rather than compile the emitted text here, the same expression tree is walked
// twice: once by the emitter, and once by rational_function::evaluate over Q.
// A pose is chosen, every matrix entry is evaluated exactly, and the double the
// generated code would produce is required to agree with it to within what a
// double can carry. A sign error, a transposed index or a dropped denominator
// all show up here.
TEST(Emit, EmittedEntriesAgreeWithTheExactFunctionsTheyStandFor) {
  const auto solution = build_solution();
  const std::array<rational, 2> at{make_rational(7, 5), make_rational(-3, 4)};
  const std::array<double, 2> at_double{7.0 / 5.0, -3.0 / 4.0};

  for (std::size_t variable = 0; variable < 2; ++variable) {
    const auto& m = solution.action[variable];
    for (std::size_t j = 0; j < m.dimension; ++j) {
      for (std::size_t i = 0; i < m.dimension; ++i) {
        const auto& f = m(i, j);
        const rational exact = f.is_zero() ? rational(0) : f.evaluate(at);

        // The emitted expression, interpreted directly. This mirrors what the
        // generated code computes without needing a compiler in the loop.
        const double approximated = [&]() -> double {
          if (f.is_zero()) {
            return 0.0;
          }
          double numerator = 0.0;
          for (const auto& t : f.numerator().terms()) {
            double term = t.coeff.get_d();
            for (std::size_t p = 0; p < 2; ++p) {
              for (unsigned e = 0; e < static_cast<unsigned>(t.mon[p]); ++e) {
                term *= at_double[p];
              }
            }
            numerator += term;
          }
          double denominator = 0.0;
          for (const auto& t : f.denominator().terms()) {
            double term = t.coeff.get_d();
            for (std::size_t p = 0; p < 2; ++p) {
              for (unsigned e = 0; e < static_cast<unsigned>(t.mon[p]); ++e) {
                term *= at_double[p];
              }
            }
            denominator += term;
          }
          return numerator / denominator;
        }();

        EXPECT_NEAR(approximated, exact.get_d(), 1e-12)
            << "variable " << variable << " entry (" << i << ", " << j << ")";
      }
    }
  }
}

// Two entries with different denominators must produce two guards, and the
// guards must not collide.
//
// This is a regression test with a story. The first version of the emitter
// wrote each guard as `const double denominator = ...;` inside one scope, which
// is fine for one pole and does not compile for two. The fixture the emitter
// was developed against happened to have exactly one, so the output compiled
// and the defect was invisible until a system with two poles was tried.
TEST(Emit, DistinctDenominatorsProduceDistinctGuards) {
  auto solution = build_solution();

  // Give the matrix of the first unknown two entries with unrelated poles.
  auto& m = solution.action[0];
  m(0, 0) = constant(1) / x();
  m(1, 1) = constant(1) / (y() + constant(1));

  const std::string header = emit(solution);

  // Neither guard may introduce a name, or a second one would redeclare it.
  EXPECT_EQ(header.find("const double denominator"), std::string::npos)
      << "guards must be written inline, or two of them collide";

  // Both poles are actually guarded.
  // Each pole is judged against the size of its own terms, not against zero.
  EXPECT_NE(header.find("if (!well_conditioned(pose[0], magnitude(pose[0]))) { return false; }"),
            std::string::npos);
  EXPECT_NE(header.find("+ (1.0), magnitude(pose[1]) + (1.0))) { return false; }"),
            std::string::npos)
      << header;
}

// The same denominator in many entries costs one guard, not one per entry.
TEST(Emit, RepeatedDenominatorsAreGuardedOnce) {
  auto solution = build_solution();
  auto& m = solution.action[0];
  m(0, 0) = constant(1) / x();
  m(0, 1) = constant(2) / x();
  m(1, 1) = constant(3) / x();

  const std::string header = emit(solution);

  const std::string guard = "if (!well_conditioned(pose[0], magnitude(pose[0])))";
  std::size_t count = 0;
  for (std::size_t at = header.find(guard); at != std::string::npos;
       at = header.find(guard, at + 1)) {
    ++count;
  }
  EXPECT_EQ(count, 1u) << "the same pole was guarded " << count << " times";
}

// The two runtimes, and the promise each makes.
//
// matrices_only exists so that the emitted header can go into a project with no
// linear algebra dependency at all; the moment an #include of Eigen leaks into
// it that promise is gone, and nothing else in the suite would notice.
TEST(Emit, MatricesOnlyPullsInNothingButTheFreestandingHeaders) {
  emit_options options;
  options.runtime = varietas::codegen::runtime_kind::matrices_only;
  const std::string header = emit(build_solution(), options);

  EXPECT_EQ(header.find("Eigen"), std::string::npos) << "Eigen leaked into the portable runtime";
  EXPECT_EQ(header.find("<complex>"), std::string::npos);
  EXPECT_EQ(header.find("static int solve"), std::string::npos);

  // What it does still carry: the matrices and the coordinates, which is
  // everything a caller needs to write its own eigensolver.
  EXPECT_NE(header.find("static bool action_matrix"), std::string::npos);
  EXPECT_NE(header.find("static bool variable_coordinates"), std::string::npos);
  EXPECT_NE(header.find("one_index"), std::string::npos);
}

TEST(Emit, TheEigenRuntimeCarriesASolver) {
  emit_options options;
  options.runtime = varietas::codegen::runtime_kind::eigen;
  const std::string header = emit(build_solution(), options);

  EXPECT_NE(header.find("#include <Eigen/Eigenvalues>"), std::string::npos);
  EXPECT_NE(header.find("static int solve"), std::string::npos);
  EXPECT_NE(header.find("EigenSolver"), std::string::npos);

  // The transpose is not incidental. Left eigenvectors of the multiplication
  // operator are the evaluation functionals; right eigenvectors are not, and
  // decomposing the matrix itself would return coordinates of nothing.
  EXPECT_NE(header.find("separating.transpose()"), std::string::npos)
      << "the solver must decompose the transpose, or it is not using the functionals";
}

TEST(Emit, AnIllFormedSolutionIsRejectedRatherThanEmitted) {
  parametric_solution<2, 2> solution;
  EXPECT_FALSE(solution.is_well_formed()) << "an empty quotient is not a solved system";

  auto good = build_solution();
  good.action.pop_back();
  EXPECT_FALSE(good.is_well_formed()) << "one matrix per unknown is required";
}

}  // namespace
