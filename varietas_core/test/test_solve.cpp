#include <algorithm>
#include <array>
#include <cmath>
#include <vector>

#include <gtest/gtest.h>

#include "varietas/core/ideal/buchberger.hpp"
#include "varietas/core/ideal/ideal.hpp"
#include "varietas/core/order/grevlex.hpp"
#include "varietas/core/polynomial.hpp"
#include "varietas/core/quotient/action_matrix.hpp"
#include "varietas/core/quotient/quotient_basis.hpp"
#include "varietas/core/solve/spectral.hpp"

namespace {

using varietas::grevlex;

using poly2 = varietas::polynomial<double, 2, grevlex>;
using mon2 = varietas::monomial<2>;

poly2 x() { return poly2::variable(0); }
poly2 y() { return poly2::variable(1); }
poly2 constant(double c) { return poly2::constant(c); }

// The unit circle intersected with the line y = x: two points, at plus and
// minus the square root of one half in both coordinates.
std::vector<poly2> circle_and_line() {
  return {x() * x() + y() * y() - constant(1.0), x() - y()};
}

TEST(Quotient, StandardMonomialsOfAZeroDimensionalIdeal) {
  const auto basis = varietas::groebner_basis(circle_and_line());
  const auto quotient = varietas::standard_monomials(basis);

  ASSERT_TRUE(quotient.is_zero_dimensional);
  // Bezout: two curves of degree two and one meet in two points.
  EXPECT_EQ(quotient.dimension(), 2u);
  EXPECT_TRUE(quotient.contains(mon2::one()));
  EXPECT_EQ(quotient.index_of(mon2::one()), 0u);

  // No standard monomial is divisible by a leading monomial of the basis.
  for (const auto& m : quotient.monomials) {
    for (const auto& g : basis) {
      EXPECT_FALSE(mon2::divides(g.leading_monomial(), m));
    }
  }
}

TEST(Quotient, PositiveDimensionalIdealIsRejected) {
  // A single conic in the plane is a curve, so the quotient is infinite
  // dimensional and no action matrix exists.
  const auto basis = varietas::groebner_basis(
      std::vector<poly2>{x() * x() + y() * y() - constant(1.0)});
  const auto quotient = varietas::standard_monomials(basis);

  EXPECT_FALSE(quotient.is_zero_dimensional);
  EXPECT_TRUE(quotient.monomials.empty());
}

TEST(Quotient, UnitIdealHasTheZeroRingAsQuotient) {
  const auto basis = varietas::groebner_basis(
      std::vector<poly2>{x() - constant(1.0), x() - constant(2.0)});
  const auto quotient = varietas::standard_monomials(basis);

  EXPECT_TRUE(quotient.is_zero_dimensional);
  EXPECT_EQ(quotient.dimension(), 0u);
}

// The multiplication operators commute, and each satisfies the minimal
// polynomial of its variable on the quotient. Both are checked here because
// they are the properties the generated header depends on.
TEST(ActionMatrix, OperatorsCommuteAndRepresentMultiplication) {
  const auto basis = varietas::groebner_basis(circle_and_line());
  const auto quotient = varietas::standard_monomials(basis);
  ASSERT_TRUE(quotient.is_zero_dimensional);

  const auto mx = varietas::variable_action_matrix<double, 2, grevlex>(0, basis, quotient);
  const auto my = varietas::variable_action_matrix<double, 2, grevlex>(1, basis, quotient);

  EXPECT_LT((mx * my - my * mx).cwiseAbs().maxCoeff(), 1e-12);

  // x and y agree on the variety, so their operators agree on the quotient.
  EXPECT_LT((mx - my).cwiseAbs().maxCoeff(), 1e-12);

  // Multiplication by x twice equals multiplication by x^2, which reduces to
  // one half on this ideal.
  const auto identity = varietas::action_matrix_type::Identity(mx.rows(), mx.cols());
  EXPECT_LT((mx * mx - 0.5 * identity).cwiseAbs().maxCoeff(), 1e-12);
}

TEST(Spectral, SolvesTheIntersectionOfACircleAndALine) {
  const auto generators = circle_and_line();
  const auto basis = varietas::groebner_basis(generators);

  const auto solutions = varietas::solve_zero_dimensional(basis);
  ASSERT_TRUE(solutions.ok()) << varietas::to_string(solutions.status);
  ASSERT_EQ(solutions.points.size(), 2u);

  for (const auto& p : solutions.points) {
    EXPECT_LT(varietas::residual(generators, p), 1e-9);
  }

  auto real = solutions.real_points();
  ASSERT_EQ(real.size(), 2u);
  std::sort(real.begin(), real.end(),
            [](const std::array<double, 2>& a, const std::array<double, 2>& b) {
              return a[0] < b[0];
            });

  const double r = std::sqrt(0.5);
  EXPECT_NEAR(real[0][0], -r, 1e-9);
  EXPECT_NEAR(real[0][1], -r, 1e-9);
  EXPECT_NEAR(real[1][0], r, 1e-9);
  EXPECT_NEAR(real[1][1], r, 1e-9);
}

TEST(Spectral, ReportsComplexSolutionsSeparatelyFromRealOnes) {
  // The circle and the line y = x + 2 do not meet over the reals, but they
  // meet in two conjugate complex points. Reporting them is what distinguishes
  // an unreachable pose from an inconsistent system.
  const std::vector<poly2> generators{x() * x() + y() * y() - constant(1.0),
                                      y() - x() - constant(2.0)};
  const auto basis = varietas::groebner_basis(generators);

  const auto solutions = varietas::solve_zero_dimensional(basis);
  ASSERT_TRUE(solutions.ok()) << varietas::to_string(solutions.status);
  EXPECT_EQ(solutions.points.size(), 2u);
  EXPECT_TRUE(solutions.real_points().empty());

  for (const auto& p : solutions.points) {
    EXPECT_LT(varietas::residual(generators, p), 1e-9);
  }
}

TEST(Spectral, SolvesAThreeVariableSystemWithFourPoints) {
  using poly3 = varietas::polynomial<double, 3, grevlex>;
  const auto u = poly3::variable(0);
  const auto v = poly3::variable(1);
  const auto w = poly3::variable(2);

  // Two circles in the first two variables and a plane pinning the third.
  const std::vector<poly3> generators{
      u * u + v * v - poly3::constant(1.0),
      u * u - v * v - poly3::constant(0.5),
      w - u - poly3::constant(1.0),
  };
  const auto basis = varietas::groebner_basis(generators);

  const auto solutions = varietas::solve_zero_dimensional(basis);
  ASSERT_TRUE(solutions.ok()) << varietas::to_string(solutions.status);
  EXPECT_EQ(solutions.points.size(), 4u);

  for (const auto& p : solutions.points) {
    EXPECT_LT(varietas::residual(generators, p), 1e-8);
  }
  EXPECT_EQ(solutions.real_points().size(), 4u);
}

TEST(Spectral, ReportsTheViolatedHypothesisInsteadOfGuessing) {
  const auto curve = varietas::groebner_basis(
      std::vector<poly2>{x() * x() + y() * y() - constant(1.0)});
  EXPECT_EQ(varietas::solve_zero_dimensional(curve).status,
            varietas::solve_status::positive_dimensional);

  const auto inconsistent = varietas::groebner_basis(
      std::vector<poly2>{x() - constant(1.0), x() - constant(2.0)});
  EXPECT_EQ(varietas::solve_zero_dimensional(inconsistent).status,
            varietas::solve_status::empty_variety);
}

TEST(Spectral, IdealFacadeExposesTheSamePipeline) {
  const varietas::ideal<double, 2, grevlex> i(circle_and_line());

  ASSERT_TRUE(i.is_zero_dimensional());
  EXPECT_EQ(i.quotient().dimension(), 2u);
  EXPECT_FALSE(i.is_unit());
  EXPECT_TRUE(i.contains(x() * x() * 2.0 - constant(1.0)));

  const auto solutions = varietas::solve_zero_dimensional(i.basis());
  EXPECT_EQ(solutions.real_points().size(), 2u);

  // The statistics of the offline run are available for inspection.
  EXPECT_GT(i.statistics().pairs_generated, 0u);
}

// Inverse kinematics of a planar two-revolute arm with unit links, posed as a
// polynomial system in the cosines and sines of the joint angles subject to the
// Pythagorean relations. This is the smallest instance of the problem the
// library exists to solve, and it exercises the whole pipeline: ideal, Gröbner
// basis, finiteness verdict, action matrix, spectral solve.
//
// Variables are ordered (c1, s1, c2, s2). The forward map is
//
//     px = c1 + (c1 c2 - s1 s2),   py = s1 + (s1 c2 + c1 s2).
//
// A reachable interior pose admits exactly two configurations, elbow up and
// elbow down, and the quotient algebra has dimension two accordingly.
TEST(Spectral, PlanarTwoLinkInverseKinematics) {
  using poly4 = varietas::polynomial<double, 4, grevlex>;

  const auto c1 = poly4::variable(0);
  const auto s1 = poly4::variable(1);
  const auto c2 = poly4::variable(2);
  const auto s2 = poly4::variable(3);

  const double px = 1.0;
  const double py = 1.0;

  const std::vector<poly4> generators{
      c1 * c1 + s1 * s1 - poly4::constant(1.0),
      c2 * c2 + s2 * s2 - poly4::constant(1.0),
      c1 + (c1 * c2 - s1 * s2) - poly4::constant(px),
      s1 + (s1 * c2 + c1 * s2) - poly4::constant(py),
  };

  const varietas::ideal<double, 4, grevlex> arm(generators);
  ASSERT_TRUE(arm.is_zero_dimensional());
  EXPECT_EQ(arm.quotient().dimension(), 2u);

  const auto solutions = varietas::solve_zero_dimensional(arm.basis());
  ASSERT_TRUE(solutions.ok()) << varietas::to_string(solutions.status);

  auto real = solutions.real_points(1e-7);
  ASSERT_EQ(real.size(), 2u);

  for (const auto& q : real) {
    // The recovered values are genuine sines and cosines.
    EXPECT_NEAR(q[0] * q[0] + q[1] * q[1], 1.0, 1e-7);
    EXPECT_NEAR(q[2] * q[2] + q[3] * q[3], 1.0, 1e-7);
    // And they reproduce the commanded pose.
    EXPECT_NEAR(q[0] + (q[0] * q[2] - q[1] * q[3]), px, 1e-7);
    EXPECT_NEAR(q[1] + (q[1] * q[2] + q[0] * q[3]), py, 1e-7);
  }

  // The two configurations are the reflections of one another in the second
  // joint: same elbow cosine, opposite elbow sine.
  EXPECT_NEAR(real[0][2], real[1][2], 1e-7);
  EXPECT_NEAR(real[0][3], -real[1][3], 1e-7);
}

// A pose beyond the reach of the arm has no solution over the reals, yet the
// system remains consistent over the complex numbers: the ideal is proper and
// the solver returns complex configurations rather than claiming failure.
TEST(Spectral, UnreachablePoseIsComplexRatherThanInconsistent) {
  using poly4 = varietas::polynomial<double, 4, grevlex>;

  const auto c1 = poly4::variable(0);
  const auto s1 = poly4::variable(1);
  const auto c2 = poly4::variable(2);
  const auto s2 = poly4::variable(3);

  const std::vector<poly4> generators{
      c1 * c1 + s1 * s1 - poly4::constant(1.0),
      c2 * c2 + s2 * s2 - poly4::constant(1.0),
      c1 + (c1 * c2 - s1 * s2) - poly4::constant(3.0),
      s1 + (s1 * c2 + c1 * s2),
  };

  const varietas::ideal<double, 4, grevlex> arm(generators);
  EXPECT_FALSE(arm.is_unit());
  ASSERT_TRUE(arm.is_zero_dimensional());

  const auto solutions = varietas::solve_zero_dimensional(arm.basis());
  ASSERT_TRUE(solutions.ok()) << varietas::to_string(solutions.status);
  EXPECT_FALSE(solutions.points.empty());
  EXPECT_TRUE(solutions.real_points(1e-7).empty());
}

}  // namespace
