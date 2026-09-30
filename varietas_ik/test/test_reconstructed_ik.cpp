// The parametric solve recovered from fixed poses.
//
// Two claims, checked separately. Where the symbolic solve over Q(p) is cheap
// enough to run, the reconstruction has to agree with it entry for entry, as
// rational functions and not merely as values. And where it is not, on a
// three-joint arm the decoupling cannot touch, the reconstruction has to be
// right exactly: at a configuration with rational half-angle tangents, the
// vector of standard monomials must be an eigenvector of every action matrix
// transposed, with the configuration's own coordinate as eigenvalue, over Q.

#include <array>
#include <cstddef>

#include <gtest/gtest.h>

#include "varietas/ik/decoupled_ik.hpp"
#include "varietas/ik/parametric_ik.hpp"
#include "varietas/ik/reconstructed_ik.hpp"
#include "varietas/kinematics/rationalize.hpp"

#include "arms.hpp"

namespace {

using varietas::make_rational;
using varietas::rational;

template <std::size_t N, std::size_t P>
void expect_same_solution(const varietas::codegen::parametric_solution<N, P>& a,
                          const varietas::codegen::parametric_solution<N, P>& b) {
  ASSERT_EQ(a.quotient.monomials, b.quotient.monomials);
  EXPECT_EQ(a.one_index, b.one_index);
  for (std::size_t i = 0; i < N; ++i) {
    ASSERT_EQ(a.action[i].entries.size(), b.action[i].entries.size());
    for (std::size_t e = 0; e < a.action[i].entries.size(); ++e) {
      EXPECT_EQ(a.action[i].entries[e], b.action[i].entries[e]) << "matrix " << i << " entry " << e;
    }
    for (std::size_t k = 0; k < a.dimension(); ++k) {
      EXPECT_EQ(a.variable_coordinates[i][k], b.variable_coordinates[i][k]);
    }
  }
}

TEST(ReconstructedIk, AgreesWithTheSymbolicSolveOnThePlanarArm) {
  const auto arm = varietas_test::planar_two_link();
  varietas::ik::reconstruction_report report;
  const auto reconstructed = varietas::ik::reconstructed_position_ik<2, 2>(arm, {0, 1}, &report);
  const auto symbolic = varietas::ik::parametric_position_ik<2, 2>(arm, {0, 1});
  ASSERT_TRUE(reconstructed.ok());
  ASSERT_TRUE(symbolic.ok());
  EXPECT_EQ(reconstructed.branches, 2u);
  EXPECT_TRUE(report.checks_passed);
  EXPECT_GE(report.exact_checks, 1u);
  expect_same_solution(reconstructed.solution, symbolic.solution);
}

// The reduced problem the decoupling poses for the anthropomorphic arm: two
// joints against a radius and a height.
TEST(ReconstructedIk, AgreesWithTheSymbolicSolveOnTheReducedArm) {
  const auto arm = varietas_test::anthropomorphic_three_link();
  const auto decoupled = varietas::ik::decoupled_position_ik<3>(arm);
  ASSERT_TRUE(decoupled.ok());
  // The reduced chain is the arm without its base joint, posed in the plane
  // the decoupling found; rebuilt here from the same pieces.
  varietas::chain<rational> reduced("reduced");
  for (std::size_t j = 1; j < arm.joints().size(); ++j) {
    reduced.add_joint(arm.joints()[j]);
  }
  reduced.set_tool(arm.tool());
  const auto plane = decoupled.frame.plane();
  const auto reconstructed = varietas::ik::reconstructed_position_ik<2, 2>(reduced, plane);
  const auto symbolic = varietas::ik::parametric_position_ik<2, 2>(reduced, plane);
  ASSERT_TRUE(reconstructed.ok());
  ASSERT_TRUE(symbolic.ok());
  expect_same_solution(reconstructed.solution, symbolic.solution);
}

// The arm whose base is displaced off its own axis, which the decoupling
// refuses, solved in full over Q(x, y, z) by reconstruction.
TEST(ReconstructedIk, SolvesAnArmTheDecouplingRefusesAndIsExactAtRationalConfigurations) {
  const auto arm = varietas_test::anthropomorphic_off_axis();
  ASSERT_FALSE(varietas::ik::decoupled_position_ik<3>(arm).ok());

  varietas::ik::reconstruction_report report;
  const auto result = varietas::ik::reconstructed_position_ik<3, 3>(arm, {0, 1, 2}, &report);
  ASSERT_TRUE(result.ok());
  EXPECT_EQ(result.branches, 4u);
  EXPECT_TRUE(report.checks_passed);

  const auto& s = result.solution;
  const auto map = varietas::rational_forward_kinematics<3, varietas::grevlex>(arm.fold_fixed_joints());
  const std::array<std::array<rational, 3>, 2> configurations{{
      {make_rational(1, 2), make_rational(-1, 3), make_rational(2, 5)},
      {make_rational(-3, 7), make_rational(5, 4), make_rational(1, 9)},
  }};
  for (const auto& t : configurations) {
    const rational denominator = map.denominator().evaluate(t);
    std::array<rational, 3> pose;
    for (std::size_t k = 0; k < 3; ++k) {
      pose[k] = map.translation(k).evaluate(t) / denominator;
    }
    const std::size_t d = s.dimension();
    std::vector<rational> w(d);
    for (std::size_t k = 0; k < d; ++k) {
      w[k] = s.quotient.monomials[k].evaluate<rational>(t);
    }
    for (std::size_t i = 0; i < 3; ++i) {
      for (std::size_t j = 0; j < d; ++j) {
        rational lhs = 0;
        for (std::size_t k = 0; k < d; ++k) {
          if (!s.action[i](k, j).is_zero()) {
            lhs += s.action[i](k, j).evaluate(pose) * w[k];
          }
        }
        EXPECT_EQ(lhs, t[i] * w[j]) << "unknown " << i << ", column " << j;
      }
    }
  }

  // The equations as posed travel with the solution, for the Newton steps.
  EXPECT_EQ(s.residual_numerators.size(), 3u);
  EXPECT_FALSE(s.residual_denominator.is_zero());
}

TEST(ReconstructedIk, RefusesByCountingAsTheSymbolicSolveDoes) {
  const auto arm = varietas_test::planar_three_link();
  const auto result = varietas::ik::reconstructed_position_ik<3, 2>(arm, {0, 1});
  EXPECT_EQ(result.status, varietas::ik::parametric_ik_status::underdetermined);
}

}  // namespace
