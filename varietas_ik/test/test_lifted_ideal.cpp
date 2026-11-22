// The pose ideal split at the middle of the chain, and lifted from Z/p.
//
// Three claims. The split equations saturate to the same ideal as the
// equations of the whole chain, so the lifted path, which uses them, answers
// the same question. The lifted basis is Buchberger's wherever Buchberger
// finishes. And on six joints, where Buchberger over Q does not, the lifted
// basis is right where it can be checked exactly: a configuration with
// rational half-angle tangents is a common zero of every element, over Q, and
// the count of branches is the one such an arm has.

#include <array>
#include <cstddef>
#include <vector>

#include <gtest/gtest.h>

#include "varietas/codegen/prime_field.hpp"
#include "varietas/core/order/grevlex.hpp"
#include "varietas/core/quotient/quotient_basis.hpp"
#include "varietas/ik/cast_chain.hpp"
#include "varietas/ik/lifted_ideal.hpp"
#include "varietas/kinematics/rationalize.hpp"

#include "arms.hpp"

namespace {

using varietas::grevlex;
using varietas::make_rational;
using varietas::rational;
using varietas::residue;

template <std::size_t N, class Coeff>
void split_agrees_with_the_whole_chain(const varietas::chain<Coeff>& robot,
                                       const std::array<Coeff, N>& configuration) {
  const auto map = varietas::rational_forward_kinematics<N, grevlex>(robot);
  const auto target = map.evaluate(configuration);
  const auto whole = varietas::kinematic_ideal_generators<N, grevlex>(
      robot, varietas::pose_residuals(map, target));
  const auto split = varietas::kinematic_ideal_generators<N, grevlex>(
      robot, varietas::split_pose_residuals<N, grevlex>(robot, target));
  EXPECT_EQ(split, whole);
  EXPECT_FALSE(varietas::is_unit_ideal(split));
}

// The first `joints` joints of an arm, with no tool: a chain of fewer joints
// than a pose has coordinates, whose poses are reachable only when they come
// from a configuration.
varietas::chain<rational> first_joints(const varietas::chain<rational>& arm, std::size_t joints) {
  varietas::chain<rational> out("truncated");
  for (std::size_t j = 0; j < joints; ++j) {
    out.add_joint(arm.joints()[j]);
  }
  return out;
}

TEST(SplitPose, SaturatesToTheIdealOfTheWholeChainOverQ) {
  split_agrees_with_the_whole_chain<3, rational>(
      varietas_test::anthropomorphic_off_axis().fold_fixed_joints(),
      {make_rational(1, 3), make_rational(1, 2), make_rational(-2, 7)});
  split_agrees_with_the_whole_chain<4, rational>(
      first_joints(varietas_test::industrial_six().fold_fixed_joints(), 4),
      {make_rational(1, 2), make_rational(-1, 3), make_rational(2, 5), make_rational(1, 4)});
}

TEST(SplitPose, SaturatesToTheIdealOfTheWholeChainOnSixJoints) {
  residue::select(2147483629);
  const auto robot =
      varietas::ik::cast_chain<residue>(varietas_test::skewed_wrist_six().fold_fixed_joints());
  const auto r = [](long n, long d) { return residue(n) / residue(d); };
  split_agrees_with_the_whole_chain<6, residue>(
      robot, {r(1, 2), r(-1, 3), r(2, 5), r(1, 4), r(-2, 3), r(3, 7)});
}

TEST(LiftedPose, IsBuchbergersBasisWhereBuchbergerFinishes) {
  const auto robot = first_joints(varietas_test::industrial_six().fold_fixed_joints(), 4);
  const std::array<rational, 4> configuration{make_rational(1, 2), make_rational(-1, 3),
                                              make_rational(2, 5), make_rational(1, 4)};
  const auto target =
      varietas::rational_forward_kinematics<4, grevlex>(robot).evaluate(configuration);
  varietas::modular::statistics s;
  const auto lifted = varietas::ik::lifted_pose_ideal_generators<4, grevlex>(robot, target, &s);
  EXPECT_TRUE(s.lifted);
  EXPECT_EQ(lifted, (varietas::pose_ideal_generators<4, grevlex>(robot, target)));
}

TEST(LiftedPose, IsBuchbergersBasisForAPosition) {
  const auto robot = varietas_test::anthropomorphic_off_axis().fold_fixed_joints();
  const varietas::vector3<rational> target(make_rational(3, 4), make_rational(-1, 2),
                                           make_rational(1, 3));
  EXPECT_EQ((varietas::ik::lifted_position_ideal_generators<3, grevlex>(robot, target)),
            (varietas::position_ideal_generators<3, grevlex>(robot, target)));
}

// A six-joint arm whose wrist axes meet but are skewed, at the pose of a
// rational configuration. Buchberger over Q has not been seen to finish this.
TEST(LiftedPose, SolvesASixJointArmAtAFullPoseExactly) {
  const auto robot = varietas_test::skewed_wrist_six().fold_fixed_joints();
  const std::array<rational, 6> configuration{make_rational(1, 2),  make_rational(-1, 3),
                                              make_rational(2, 5),  make_rational(1, 4),
                                              make_rational(-2, 3), make_rational(3, 7)};
  const auto target =
      varietas::rational_forward_kinematics<6, grevlex>(robot).evaluate(configuration);

  varietas::modular::statistics s;
  const auto basis = varietas::ik::lifted_pose_ideal_generators<6, grevlex>(robot, target, &s);
  ASSERT_TRUE(s.lifted);

  const auto quotient = varietas::standard_monomials(basis);
  ASSERT_TRUE(quotient.is_zero_dimensional);
  EXPECT_EQ(quotient.dimension(), 8u);  // a spherical wrist's eight branches
  for (const auto& g : basis) {
    EXPECT_EQ(g.evaluate(configuration), rational(0));
  }
}

}  // namespace
