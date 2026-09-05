// A generated inverse kinematics solver, compiled and run against the arm it
// came from.
//
// The header included below was produced during this build by a program that
// knew a chain and nothing else — no pose, no target, no numbers beyond the
// link lengths. What it has to satisfy is the only thing an inverse kinematics
// solver has to satisfy: driving the joints to a returned configuration puts
// the tool where it was asked to go. That is checked here by pushing every
// solution back through the numerical forward kinematics in
// varietas_kinematics, which is an independent computation — trigonometry on
// doubles, with no ideal and no matrix anywhere in it.
//
// The planar 2R arm of unit links reaches a point exactly when it lies in the
// annulus 0 <= r <= 2, with two elbow configurations inside it, one on each
// boundary circle, and none outside. The solver is required to agree.

#include <array>
#include <cmath>
#include <cstddef>
#include <vector>

#include <gtest/gtest.h>

#include "varietas/kinematics/evaluate.hpp"
#include "varietas/kinematics/rationalize.hpp"

#include "arms.hpp"
#include "planar_2r_ik.hpp"  // generated into the build tree by generate_planar_2r

namespace {

using solver = varietas_generated::planar_2r_ik;

varietas::chain<double> numeric_arm() {
  return varietas::chain_cast<double>(varietas_test::planar_two_link());
}

// The tool position the arm reaches when its half-angle variables are t.
std::array<double, 2> tool_position(const std::array<double, 2>& t) {
  const std::vector<double> angles{varietas::angle_from_variable(t[0]),
                                   varietas::angle_from_variable(t[1])};
  const auto frame = varietas::forward_kinematics(numeric_arm(), angles);
  return {frame.translation()[0], frame.translation()[1]};
}

// The locus the parametric basis does not describe.
//
// Eliminating under grevlex divided by x^2 + y^2 + 2x, so the one basis holds
// everywhere except the circle (x + 1)^2 + y^2 = 1. Nothing about the arm is
// special there — it reaches those points perfectly well, with two elbow
// configurations like any other interior point — the circle is an artefact of
// how the basis was computed, and it is the generated guard's business to
// refuse rather than answer there.
double pole_value(double x, double y) { return x * x + y * y + 2.0 * x; }

struct solve_result {
  int count = 0;
  std::vector<std::array<double, 2>> configurations;
  solver::status state{};
};

solve_result solve_at(double x, double y) {
  const double pose[2] = {x, y};
  std::array<double, 8> out{};
  solve_result result;
  result.count = solver::solve(pose, out.data(), 4, &result.state);
  for (int k = 0; k < result.count; ++k) {
    result.configurations.push_back(
        {out[static_cast<std::size_t>(k) * 2 + 0], out[static_cast<std::size_t>(k) * 2 + 1]});
  }
  return result;
}

// --- the header describes the arm it was generated from --------------------

TEST(GeneratedIk, ShapeAndOrderSurvivedTheRoundTrip) {
  EXPECT_EQ(solver::num_unknowns, 2u);
  EXPECT_EQ(solver::num_parameters, 2u);
  // Two elbow configurations, which is the count the offline solve certified.
  EXPECT_EQ(solver::dimension, 2u);
  EXPECT_STREQ(solver::order_name, "grevlex");
}

// --- the solutions are solutions -------------------------------------------

TEST(GeneratedIk, EveryReturnedConfigurationReachesTheRequestedPoint) {
  // A point well inside the annulus and off both axes, so neither branch is a
  // special configuration.
  const double x = 1.1;
  const double y = 0.4;

  const auto result = solve_at(x, y);
  ASSERT_GE(result.count, 0) << "solve failed";
  EXPECT_EQ(result.state, solver::status::ok);
  EXPECT_EQ(result.count, 2) << "an interior point has two elbow configurations";

  for (const auto& t : result.configurations) {
    const auto reached = tool_position(t);
    EXPECT_NEAR(reached[0], x, 1e-9);
    EXPECT_NEAR(reached[1], y, 1e-9);
  }
}

TEST(GeneratedIk, TheTwoBranchesAreDistinct) {
  const auto result = solve_at(1.1, 0.4);
  ASSERT_EQ(result.count, 2);

  const auto& a = result.configurations[0];
  const auto& b = result.configurations[1];
  // Elbow up and elbow down differ in the second joint; returning one branch
  // twice would satisfy the residual check above and still be wrong.
  EXPECT_GT(std::abs(a[1] - b[1]), 1e-6)
      << "both branches should be returned, not one of them twice";
}

TEST(GeneratedIk, AgreesWithTheForwardMapAcrossTheAnnulus) {
  // A sweep rather than a point: the parametric basis claims to hold away from
  // its pole, and one pose is a weak test of that claim.
  int poses_with_two_branches = 0;
  for (int i = -8; i <= 8; ++i) {
    for (int j = -8; j <= 8; ++j) {
      const double x = 0.2 * i;
      const double y = 0.2 * j;
      const double r = std::hypot(x, y);
      if (r < 0.15 || r > 1.9) {
        continue;  // stay off the boundary circles, where the branches merge
      }

      const auto result = solve_at(x, y);
      if (result.count < 0) {
        // The guard refused this pose. That is its business, and the tests
        // above pin down exactly which poses it refuses; here it is enough that
        // a refusal is a refusal and not a wrong answer.
        EXPECT_EQ(result.state, solver::status::bad_pose);
        EXPECT_LT(std::abs(pole_value(x, y)), 1e-9)
            << "only poses on the pole should be refused, but (" << x << ", " << y
            << ") is not one";
        continue;
      }
      ASSERT_EQ(result.count, 2) << "interior point (" << x << ", " << y << ")";
      ++poses_with_two_branches;

      for (const auto& t : result.configurations) {
        const auto reached = tool_position(t);
        EXPECT_NEAR(reached[0], x, 1e-8) << "at (" << x << ", " << y << ")";
        EXPECT_NEAR(reached[1], y, 1e-8) << "at (" << x << ", " << y << ")";
      }
    }
  }
  EXPECT_GT(poses_with_two_branches, 100) << "the sweep should not have been nearly empty";
}

// --- and refuses what the arm cannot do ------------------------------------

// --- the pole, and how narrowly the guard catches it -----------------------

TEST(GeneratedIk, RefusesAPoseExactlyOnThePole) {
  // (-1, -1) lies on the circle exactly, and in binary too: 1 + 1 - 2 is zero
  // with no rounding anywhere, so the emitted guard sees a denominator of
  // exactly 0.0 and refuses.
  ASSERT_EQ(pole_value(-1.0, -1.0), 0.0);

  const auto result = solve_at(-1.0, -1.0);
  EXPECT_EQ(result.count, -1);
  EXPECT_EQ(result.state, solver::status::bad_pose);
}

// A pose that is on the pole mathematically and not in floating point.
//
// (-1.6, 0.8) satisfies x^2 + y^2 + 2x = 0 exactly over the rationals, but
// evaluated in doubles the same expression comes to about 4e-16 rather than to
// zero. A guard comparing against 0.0 would not fire, and the action matrices
// would then be formed by dividing by that 4e-16 — which is what this header
// used to do, returning two configurations of which one did not reach the
// target and neither the count nor the status said so.
//
// The guard now compares the denominator against the size of the terms that
// produced it, so total cancellation is recognised whether or not it landed
// exactly on zero, and the pose is refused.
TEST(GeneratedIk, RefusesAPoseWithinRoundingOfThePole) {
  const double x = 0.2 * -8;  // -1.6000000000000001
  const double y = 0.2 * 4;   //  0.8
  EXPECT_NE(pole_value(x, y), 0.0) << "the premise of this test is that rounding hides the pole";
  EXPECT_LT(std::abs(pole_value(x, y)), 1e-12) << "and that it is nevertheless on it";

  const auto result = solve_at(x, y);
  EXPECT_EQ(result.count, -1) << "a pose within rounding of the pole must be refused";
  EXPECT_EQ(result.state, solver::status::bad_pose);
}

// And a pose that is merely near the pole, not on it, is still answered.
//
// A relative guard has to refuse total cancellation without refusing every pose
// in the neighbourhood, or it would carve a hole out of the workspace instead
// of a circle.
TEST(GeneratedIk, APoseNearButNotOnThePoleIsStillAnswered) {
  // Off the circle by about a hundredth, which is enormous next to 1e-12.
  const double x = -1.59;
  const double y = 0.8;
  ASSERT_GT(std::abs(pole_value(x, y)), 1e-3);

  const auto result = solve_at(x, y);
  ASSERT_GE(result.count, 0) << "the guard must not refuse a pose merely close to the pole";
  ASSERT_EQ(result.count, 2);
  for (const auto& t : result.configurations) {
    const auto reached = tool_position(t);
    EXPECT_NEAR(reached[0], x, 1e-7);
    EXPECT_NEAR(reached[1], y, 1e-7);
  }
}

TEST(GeneratedIk, PointsOutsideTheAnnulusHaveNoRealSolution) {
  // r > 2 is out of reach; the two configurations are complex there, and a
  // manipulator cannot be commanded to either.
  const auto result = solve_at(2.5, 0.0);
  ASSERT_GE(result.count, 0);
  EXPECT_EQ(result.count, 0) << "a point beyond the arm's reach has no real configuration";
}

TEST(GeneratedIk, IsDeterministic) {
  const auto a = solve_at(1.1, 0.4);
  const auto b = solve_at(1.1, 0.4);
  ASSERT_EQ(a.count, b.count);
  for (std::size_t k = 0; k < a.configurations.size(); ++k) {
    EXPECT_DOUBLE_EQ(a.configurations[k][0], b.configurations[k][0]);
    EXPECT_DOUBLE_EQ(a.configurations[k][1], b.configurations[k][1]);
  }
}

}  // namespace
