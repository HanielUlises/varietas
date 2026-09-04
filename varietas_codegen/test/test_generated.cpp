// Generated code, compiled and run.
//
// Everything else about the emitter is checked by comparing strings or by
// walking the same expression tree twice. This suite is the one that cannot be
// fooled that way: the header under test was produced by the emitter during
// this build, was handed to the compiler as ordinary source, and is called
// here. If the emitted text does not parse, the build fails; if it parses and
// computes the wrong thing, these expectations fail.
//
// The system is the one in planar_fixture.hpp, whose action matrices are known
// in closed form:
//
//   M_u = [ 0    xy ]      M_v = [ 0  xy^2 ]
//         [ 1/y  0  ]            [ 1  0    ]
//
// so the generated numbers can be checked against arithmetic done here, rather
// than against another run of the same machinery.

#include <array>
#include <cmath>
#include <cstddef>

#include <gtest/gtest.h>

#include "varietas/codegen/rational.hpp"
#include "varietas/core/ideal/ideal.hpp"
#include "varietas/core/order/grevlex.hpp"
#include "varietas/core/order/order_id.hpp"
#include "varietas/core/solve/spectral.hpp"

#include "planar.hpp"  // generated into the build tree by generate_planar

namespace {

using solver = varietas_generated::planar;

// Column-major, as the emitter documents and writes.
double entry(const std::array<double, 4>& m, std::size_t row, std::size_t column) {
  return m[column * solver::dimension + row];
}

TEST(Generated, ShapeAndOrderSurvivedTheRoundTrip) {
  EXPECT_EQ(solver::num_unknowns, 2u);
  EXPECT_EQ(solver::num_parameters, 2u);
  EXPECT_EQ(solver::dimension, 2u);

  // The recorded order must be the one the basis was actually computed under.
  // This is the check the README promises: a mismatch here is what stops a
  // header generated under one order from being used under another.
  EXPECT_EQ(solver::order_id,
            static_cast<std::uint8_t>(varietas::order_id::grevlex));
  EXPECT_STREQ(solver::order_name, "grevlex");
}

TEST(Generated, ActionMatricesMatchTheClosedForm) {
  const double x = 1.4;
  const double y = -0.75;
  const double pose[2] = {x, y};

  std::array<double, 4> mu{};
  ASSERT_TRUE(solver::action_matrix(0, pose, mu.data()));
  EXPECT_DOUBLE_EQ(entry(mu, 0, 0), 0.0);
  EXPECT_DOUBLE_EQ(entry(mu, 1, 0), 1.0 / y);
  EXPECT_DOUBLE_EQ(entry(mu, 0, 1), x * y);
  EXPECT_DOUBLE_EQ(entry(mu, 1, 1), 0.0);

  std::array<double, 4> mv{};
  ASSERT_TRUE(solver::action_matrix(1, pose, mv.data()));
  EXPECT_DOUBLE_EQ(entry(mv, 0, 0), 0.0);
  EXPECT_DOUBLE_EQ(entry(mv, 1, 0), 1.0);
  EXPECT_DOUBLE_EQ(entry(mv, 0, 1), x * y * y);
  EXPECT_DOUBLE_EQ(entry(mv, 1, 1), 0.0);
}

// The pole the parametric basis has, and the guard that reports it.
//
// Eliminating u divides by y, so at y = 0 the basis this header was generated
// from does not describe the system. The generated code has to say so: the
// alternative is a matrix full of infinities that an eigensolver would happily
// consume and answer wrongly from.
TEST(Generated, APoseOnThePoleIsRefusedRatherThanAnswered) {
  std::array<double, 4> m{};
  m.fill(-1.0);
  const double on_the_pole[2] = {1.0, 0.0};

  EXPECT_FALSE(solver::action_matrix(0, on_the_pole, m.data()))
      << "y = 0 is where u = v/y fails to be defined";
  for (double e : m) {
    EXPECT_DOUBLE_EQ(e, -1.0) << "out was written to despite the refusal";
  }

  // The v matrix has no denominator, so the same pose is perfectly fine there.
  EXPECT_TRUE(solver::action_matrix(1, on_the_pole, m.data()));
}

// The reason the matrices are worth generating: their eigenvalues are the
// coordinates of the solutions. For this system u^2 = x, so the spectrum of the
// u matrix must be the two square roots of x, and the trace and determinant are
// enough to say so without an eigensolver.
TEST(Generated, TheSpectrumOfTheUMatrixIsTheSolutionSet) {
  const double x = 2.25;
  const double pose[2] = {x, 0.5};

  std::array<double, 4> mu{};
  ASSERT_TRUE(solver::action_matrix(0, pose, mu.data()));

  const double trace = entry(mu, 0, 0) + entry(mu, 1, 1);
  const double determinant =
      entry(mu, 0, 0) * entry(mu, 1, 1) - entry(mu, 0, 1) * entry(mu, 1, 0);

  // Eigenvalues +sqrt(x) and -sqrt(x): they sum to zero and multiply to -x.
  // Note that the y in the (0,1) entry and the 1/y in the (1,0) entry cancel
  // in the determinant, which is how a basis with a pole at y = 0 still has a
  // spectrum that does not depend on y at all.
  EXPECT_NEAR(trace, 0.0, 1e-15);
  EXPECT_NEAR(determinant, -x, 1e-15);

  // Which is to say the characteristic polynomial is lambda^2 - x, the
  // generator we started from.
  const double lambda = std::sqrt(x);
  EXPECT_NEAR(lambda * lambda - x, 0.0, 1e-15);
}

TEST(Generated, PoseDependenceIsRealAndNotBakedIn) {
  std::array<double, 4> a{};
  std::array<double, 4> b{};
  const double first[2] = {1.0, 1.0};
  const double second[2] = {9.0, 1.0};

  ASSERT_TRUE(solver::action_matrix(0, first, a.data()));
  ASSERT_TRUE(solver::action_matrix(0, second, b.data()));
  EXPECT_NE(entry(a, 0, 1), entry(b, 0, 1))
      << "the matrix did not change with the pose, so nothing parametric was emitted";
}

// The coordinates that turn an eigenvector into a point.
//
// u is not a standard monomial here — it was eliminated — so its value at a
// solution is the normal form v/y evaluated there. That is what these
// coordinates are, and getting them wrong would produce a solver that returns
// the right number of points with the wrong u in each.
TEST(Generated, VariableCoordinatesAreTheNormalForms) {
  const double y = -0.4;
  const double pose[2] = {1.0, y};

  std::array<double, 2> cu{};
  ASSERT_TRUE(solver::variable_coordinates(0, pose, cu.data()));
  EXPECT_DOUBLE_EQ(cu[0], 0.0);        // no constant part
  EXPECT_DOUBLE_EQ(cu[1], 1.0 / y);    // u = v / y

  std::array<double, 2> cv{};
  ASSERT_TRUE(solver::variable_coordinates(1, pose, cv.data()));
  EXPECT_DOUBLE_EQ(cv[0], 0.0);
  EXPECT_DOUBLE_EQ(cv[1], 1.0);        // v = v
}

// The whole point, end to end: generated code, an eigendecomposition, and the
// solutions of the system it was generated from.
//
// The check is a residual, not a comparison against expected numbers. The
// solver is allowed to return the two branches in either order, and what has to
// hold of each is that it satisfies the equations we started from:
// u^2 = x and v = y u.
TEST(Generated, SolveReturnsPointsThatSatisfyTheOriginalEquations) {
  const double x = 2.25;
  const double y = 0.5;
  const double pose[2] = {x, y};

  std::array<double, 8> out{};
  solver::status state{};
  const int count = solver::solve(pose, out.data(), 4, &state);

  ASSERT_GE(count, 0) << "solve failed";
  EXPECT_EQ(state, solver::status::ok);
  EXPECT_EQ(count, 2) << "u^2 = x has two real roots at x = " << x;

  bool saw_positive = false;
  bool saw_negative = false;
  for (int k = 0; k < count; ++k) {
    const double u = out[static_cast<std::size_t>(k) * 2 + 0];
    const double v = out[static_cast<std::size_t>(k) * 2 + 1];

    EXPECT_NEAR(u * u - x, 0.0, 1e-9) << "solution " << k << " violates u^2 - x";
    EXPECT_NEAR(v - y * u, 0.0, 1e-9) << "solution " << k << " violates v - y u";

    if (u > 0.0) saw_positive = true;
    if (u < 0.0) saw_negative = true;
  }
  EXPECT_TRUE(saw_positive && saw_negative) << "both branches should be returned, not one twice";
}

// The complex case. x < 0 puts both roots off the real line, and a manipulator
// cannot be commanded to either, so the count is zero rather than two garbage
// points with the imaginary part quietly dropped.
TEST(Generated, ComplexSolutionsAreNotReportedAsReal) {
  const double pose[2] = {-1.0, 0.5};
  std::array<double, 8> out{};
  solver::status state{};

  const int count = solver::solve(pose, out.data(), 4, &state);
  ASSERT_GE(count, 0);
  EXPECT_EQ(count, 0) << "u^2 = -1 has no real solution";
}

TEST(Generated, SolveRefusesAPoseOnThePole) {
  const double on_the_pole[2] = {4.0, 0.0};
  std::array<double, 8> out{};
  solver::status state{};

  EXPECT_EQ(solver::solve(on_the_pole, out.data(), 4, &state), -1);
  EXPECT_EQ(state, solver::status::bad_pose);
}

TEST(Generated, SolveHonoursTheCapacityItIsGiven) {
  const double pose[2] = {2.25, 0.5};
  std::array<double, 8> out{};
  out.fill(-99.0);

  const int count = solver::solve(pose, out.data(), 1);
  EXPECT_EQ(count, 1) << "capacity 1 must not be exceeded";
  EXPECT_DOUBLE_EQ(out[2], -99.0) << "solve wrote past the capacity it was given";
}

TEST(Generated, SolveIsDeterministic) {
  const double pose[2] = {3.0, -1.25};
  std::array<double, 8> a{};
  std::array<double, 8> b{};
  const int first = solver::solve(pose, a.data(), 4);
  const int second = solver::solve(pose, b.data(), 4);
  ASSERT_EQ(first, second);
  for (int k = 0; k < first * 2; ++k) {
    EXPECT_DOUBLE_EQ(a[static_cast<std::size_t>(k)], b[static_cast<std::size_t>(k)])
        << "the same pose returned different points, so the separating form is not fixed";
  }
}

// Generated code against the library that generated it.
//
// Every other test here checks the generated solver against arithmetic written
// out by hand, which catches a wrong answer but not a shared misunderstanding.
// This one runs varietas_core's own spectral solver on the same system with the
// pose substituted before the computation, and requires the two to agree. The
// offline path — pose adjoined to the coefficient field, one basis for every
// pose, emitted, compiled — and the online path — a basis per pose, solved in
// place — are genuinely different computations, and their agreeing is the
// claim that makes emitting worth doing at all.
TEST(Generated, AgreesWithTheLibrarySolverOnTheSameSystem) {
  using varietas::grevlex;
  using varietas::make_rational;
  using varietas::rational;
  using numeric = varietas::polynomial<rational, 2, grevlex>;

  const double x = 2.25;
  const double y = 0.5;

  const auto u = numeric::variable(0);
  const auto v = numeric::variable(1);
  const varietas::ideal<rational, 2, grevlex> ideal({
      u * u - numeric::constant(make_rational(9, 4)),  // x = 2.25
      v - numeric::constant(make_rational(1, 2)) * u,  // y = 0.5
  });

  const auto reference = varietas::solve_zero_dimensional(ideal.basis());
  ASSERT_TRUE(reference.ok()) << varietas::to_string(reference.status);
  const auto expected = reference.real_points();

  const double pose[2] = {x, y};
  std::array<double, 8> out{};
  const int count = solver::solve(pose, out.data(), 4);
  ASSERT_GE(count, 0);
  ASSERT_EQ(static_cast<std::size_t>(count), expected.size())
      << "the two paths disagree on how many real solutions there are";

  // Neither path promises an order, so each expected point must be matched.
  for (const auto& want : expected) {
    bool matched = false;
    for (int k = 0; k < count && !matched; ++k) {
      const double du = out[static_cast<std::size_t>(k) * 2 + 0] - want[0];
      const double dv = out[static_cast<std::size_t>(k) * 2 + 1] - want[1];
      matched = std::abs(du) < 1e-9 && std::abs(dv) < 1e-9;
    }
    EXPECT_TRUE(matched) << "the library found (" << want[0] << ", " << want[1]
                         << ") and the generated solver did not";
  }
}

TEST(Generated, AnUnknownOutOfRangeIsRefused) {
  std::array<double, 4> m{};
  const double pose[2] = {1.0, 1.0};
  EXPECT_FALSE(solver::action_matrix(2, pose, m.data()));
  EXPECT_FALSE(solver::action_matrix(99, pose, m.data()));
}

}  // namespace
