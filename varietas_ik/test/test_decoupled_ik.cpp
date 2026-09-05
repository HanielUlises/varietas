// Solving an arm by not adjoining its first joint.
//
// The claim being tested is a cost claim as much as a correctness one: the
// three-joint arm that produced nothing over Q(x, y, z) in fifteen minutes
// reduces to a two-joint problem that solves in well under a second, and the
// branch count of the reduced problem, doubled, is the branch count of the arm.
// The geometry is checked here; that the reduced solution generates code which
// actually reaches the requested point is checked in test_generated_decoupled.

#include <cstddef>

#include <gtest/gtest.h>

#include "varietas/ik/decoupled_ik.hpp"

#include "arms.hpp"

namespace {

using varietas::ik::decoupled_position_ik;
using varietas::ik::decoupling_status;

TEST(decoupled_ik, the_anthropomorphic_arm_reduces) {
  const auto result = decoupled_position_ik<3>(varietas_test::anthropomorphic_three_link());

  ASSERT_TRUE(result.ok()) << varietas::ik::to_string(result.status);

  // Turning about z sweeps x into y and leaves z alone, so the reduced problem
  // is posed in x and z and it is y that the sweep generates.
  EXPECT_EQ(result.frame.axis, 2u);
  EXPECT_EQ(result.frame.radial, 0u);
  EXPECT_EQ(result.frame.swept, 1u);
  EXPECT_FALSE(result.frame.reversed);

  // Two joints against two coordinates, which is the size that solves quickly.
  EXPECT_EQ(result.reduced.dimension(), 2u);
  EXPECT_TRUE(result.reduced.is_well_formed());

  // Elbow up and elbow down, each reachable facing forwards or reversed: the
  // four branches an anthropomorphic arm is known to have, and the same count
  // that solving the whole system over Q at a fixed pose reports.
  EXPECT_EQ(result.branches, 4u);
  EXPECT_EQ(result.first_joint_name, "q1");
}

TEST(decoupled_ik, a_base_offset_along_its_own_axis_is_still_decoupled) {
  // A translation along z commutes with a rotation about z, so the pedestal
  // changes the height the reduced problem is posed at and nothing else.
  const auto result = decoupled_position_ik<3>(varietas_test::anthropomorphic_on_a_pedestal());

  ASSERT_TRUE(result.ok()) << varietas::ik::to_string(result.status);
  EXPECT_EQ(result.branches, 4u);
  EXPECT_EQ(result.frame.axis, 2u);
}

TEST(decoupled_ik, the_reduced_problem_carries_the_remaining_joints) {
  const auto result = decoupled_position_ik<3>(varietas_test::anthropomorphic_three_link());
  ASSERT_TRUE(result.ok());

  ASSERT_EQ(result.reduced.unknown_names.size(), 2u);
  EXPECT_EQ(result.reduced.unknown_names[0], "t_q2");
  EXPECT_EQ(result.reduced.unknown_names[1], "t_q3");

  // The parameters are the radius and the height, in that order.
  ASSERT_EQ(result.reduced.parameter_names.size(), 2u);
  EXPECT_EQ(result.reduced.parameter_names[0], "x");
  EXPECT_EQ(result.reduced.parameter_names[1], "z");
}

// --- what does not decouple ------------------------------------------------

TEST(decoupled_ik, a_planar_arm_does_not_decouple) {
  // The base of a planar arm turns in the plane the arm already works in, so
  // removing it leaves an arm that still moves in y. There is no sweep.
  const auto result = decoupled_position_ik<2>(varietas_test::planar_two_link());

  EXPECT_EQ(result.status, decoupling_status::does_not_reduce)
      << varietas::ik::to_string(result.status);
  // The geometry is refused on its own terms, before the reduced problem is
  // posed at all, so there is no reduced status to report.
  EXPECT_EQ(result.reduced_status, varietas::ik::parametric_ik_status::ok);
}

TEST(decoupled_ik, a_base_off_its_own_axis_does_not_decouple) {
  // The placement no longer commutes with the rotation: turning the base
  // carries the axis around with it, and no plane stays fixed.
  const auto result = decoupled_position_ik<3>(varietas_test::anthropomorphic_off_axis());

  EXPECT_EQ(result.status, decoupling_status::first_joint_is_displaced)
      << varietas::ik::to_string(result.status);
}

TEST(decoupled_ik, refuses_a_chain_with_the_wrong_number_of_joints) {
  const auto result = decoupled_position_ik<3>(varietas_test::planar_two_link());
  EXPECT_EQ(result.status, decoupling_status::wrong_degrees_of_freedom);
}

}  // namespace
