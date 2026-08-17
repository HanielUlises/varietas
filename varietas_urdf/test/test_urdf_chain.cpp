#include <cmath>
#include <random>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include <kdl/chain.hpp>
#include <kdl/chainfksolverpos_recursive.hpp>
#include <kdl/frames.hpp>
#include <kdl/tree.hpp>
#include <kdl_parser/kdl_parser.hpp>
#include <urdf/model.h>

#include "varietas/core/order/grevlex.hpp"
#include "varietas/kinematics/evaluate.hpp"
#include "varietas/kinematics/rationalize.hpp"
#include "varietas/urdf/urdf_chain.hpp"

namespace {

using varietas::chain;
using varietas::chain_cast;
using varietas::forward_kinematics;
using varietas::rational;
using varietas::rigid_transform;
using varietas::urdf_import::chain_from_model;
using varietas::urdf_import::import_report;
using varietas::urdf_import::import_status;

constexpr const char* kRoot = "lbr_iiwa_link_0";
constexpr const char* kTip = "lbr_iiwa_link_7";

std::string fixture_path() { return std::string(VARIETAS_URDF_TEST_DATA) + "/lbr_iiwa14.urdf"; }

urdf::Model load_model() {
  urdf::Model model;
  EXPECT_TRUE(model.initFile(fixture_path())) << fixture_path();
  return model;
}

// --- the rational recovery -------------------------------------------------

TEST(rational_approximation, recovers_the_decimal_the_author_wrote) {
  using varietas::urdf_import::rationalize;

  // 0.1575 is not a binary fraction, so reading it as a double and asking for
  // the exact rational gives a number with a fifty digit numerator; asking for
  // the best rational of bounded denominator gives back the decimal itself.
  EXPECT_EQ(rationalize(0.1575), varietas::make_rational(63, 400));
  EXPECT_EQ(rationalize(0.2045), varietas::make_rational(409, 2000));
  EXPECT_EQ(rationalize(0.081), varietas::make_rational(81, 1000));
  EXPECT_EQ(rationalize(-0.5), varietas::make_rational(-1, 2));
  EXPECT_EQ(rationalize(0.0), varietas::make_rational(0));

  // An irrational is approximated, not recovered, and the bound is respected.
  const rational third = rationalize(1.0 / 3.0);
  EXPECT_EQ(third, varietas::make_rational(1, 3));
  const rational pi = rationalize(M_PI, 100);
  EXPECT_EQ(pi, varietas::make_rational(311, 99));
  EXPECT_LE(pi.get_den(), 100);
}

TEST(rational_approximation, right_angles_are_recovered_exactly_from_truncated_decimals) {
  using varietas::urdf_import::snap_rotation;

  // The quaternion of a quarter turn about x, as the iiwa file's truncated
  // 1.57079632679 produces it: irrational entries, and yet projectively
  // integral, so the exact rotation is recovered with no deviation of its own.
  const double half = 1.57079632679 / 2.0;
  const auto snap = snap_rotation(std::sin(half), 0.0, 0.0, std::cos(half));

  EXPECT_EQ(snap.w, rational(1));
  EXPECT_EQ(snap.x, rational(1));
  EXPECT_EQ(snap.y, rational(0));
  EXPECT_EQ(snap.z, rational(0));

  // What is left is the distance from the file's angle to the right angle it
  // meant, which is the file's own truncation and not an error of ours. The
  // recovery therefore does move the robot, and says so.
  EXPECT_GT(snap.deviation_radians, 0.0);
  EXPECT_LT(snap.deviation_radians, 1e-11);
  EXPECT_FALSE(snap.round_trips);

  // And what it moves onto is exactly a rotation, which the stated one was not.
  EXPECT_DOUBLE_EQ(orthogonality_defect(snap.exact()), 0.0);
}

TEST(rational_approximation, an_oblique_rotation_is_reported_rather_than_rounded) {
  using varietas::urdf_import::snap_rotation;

  // A rotation by one radian about z has no exact rational quaternion, so the
  // recovery leaves a deviation far above any sane tolerance instead of
  // quietly producing a robot that is not the one described.
  const auto snap = snap_rotation(0.0, 0.0, std::sin(0.5), std::cos(0.5), 1000);
  EXPECT_GT(snap.deviation_radians, 1e-9);
  // It is still an exact rotation, just not the right one.
  EXPECT_DOUBLE_EQ(orthogonality_defect(snap.exact()), 0.0);
}

// --- importing the model ---------------------------------------------------

TEST(urdf_chain, imports_the_iiwa_exactly) {
  urdf::Model model = load_model();

  chain<rational> robot;
  const import_report report = chain_from_model(model, kRoot, kTip, robot);

  ASSERT_TRUE(report.ok()) << to_string(report.status) << ": " << report.detail;
  EXPECT_EQ(robot.degrees_of_freedom(), 7u);
  EXPECT_EQ(robot.joints().size(), 7u);
  EXPECT_EQ(robot.joints().front().name, "lbr_iiwa_joint_1");
  EXPECT_EQ(robot.joints().back().name, "lbr_iiwa_joint_7");

  EXPECT_EQ(report.joints.size(), 7u);

  // Every rotation in the file is a multiple of a right angle and every length
  // a short decimal, so all of it is representable. The lengths come back
  // indistinguishable from the file; the rotations do not, and that is the
  // recovery working rather than failing — it has moved each axis onto the
  // right angle the file truncated, by about 4e-12 radians.
  EXPECT_GT(report.max_rotation_deviation, 0.0);
  EXPECT_LT(report.max_rotation_deviation, 1e-11);
  EXPECT_FALSE(report.indistinguishable_from_file());

  // The translation residues are the file's decimal-to-binary rounding, below
  // the resolution of a double, so no length has moved.
  EXPECT_LT(report.max_translation_deviation, 1e-16);

  // The one joint written as an exact rpy of zeros is recovered without moving
  // at all, which is the case the flag is there to distinguish.
  EXPECT_EQ(report.joints[0].rotation_deviation, 0.0);
  EXPECT_TRUE(report.joints[0].unmoved);
  EXPECT_FALSE(report.joints[1].unmoved);

  // The recovered chain is valid over the exact field, which is the property
  // that a rounded import would fail: every origin is exactly orthogonal.
  EXPECT_TRUE(robot.validate().ok());

  // And the recovered lengths are the decimals of the file, not doubles.
  EXPECT_EQ(robot.joints()[0].origin.translation()[2], varietas::make_rational(63, 400));
  EXPECT_EQ(robot.joints()[6].origin.translation()[1], varietas::make_rational(81, 1000));
}

TEST(urdf_chain, carries_the_joint_limits_across) {
  urdf::Model model = load_model();

  chain<rational> robot;
  ASSERT_TRUE(chain_from_model(model, kRoot, kTip, robot).ok());

  for (const auto& j : robot.joints()) {
    EXPECT_TRUE(j.has_limits) << j.name;
    EXPECT_LT(j.lower, j.upper);
  }
  EXPECT_NEAR(robot.joints()[0].lower, -2.96705972839, 1e-12);
  EXPECT_NEAR(robot.joints()[1].upper, 2.09439510239, 1e-12);
}

TEST(urdf_chain, refuses_a_model_whose_geometry_is_not_exactly_representable) {
  urdf::Model model = load_model();

  // A tolerance below the file's own truncation of pi makes the model
  // unrepresentable, and the importer says so and names the joint instead of
  // returning a chain that is exact about the wrong robot.
  varietas::urdf_import::import_options strict;
  strict.rotation_tolerance = 1e-15;

  chain<rational> robot;
  const import_report report = chain_from_model(model, kRoot, kTip, robot, strict);

  EXPECT_FALSE(report.ok());
  EXPECT_EQ(report.status, import_status::rotation_deviation_exceeded);
  EXPECT_EQ(report.detail, "lbr_iiwa_joint_2");
  EXPECT_EQ(robot.degrees_of_freedom(), 0u) << "the chain must be left untouched";
}

TEST(urdf_chain, names_a_tip_that_is_not_in_the_model) {
  urdf::Model model = load_model();

  chain<rational> robot;
  const import_report report = chain_from_model(model, kRoot, "no_such_link", robot);

  EXPECT_EQ(report.status, import_status::tip_link_not_found);
  EXPECT_EQ(report.detail, "no_such_link");
}

TEST(urdf_chain, finds_the_tip_of_a_serial_model) {
  urdf::Model model = load_model();
  EXPECT_EQ(varietas::urdf_import::sole_tip_link(model), kTip);
}

// --- agreement with the reference implementation ---------------------------

// The point of the exercise: the chain recovered exactly from the file has the
// same forward kinematics as the one KDL builds from that file by reading the
// decimals directly. If the recovery had changed the robot, this is where it
// would show, and it is checked over the whole configuration space rather than
// at one pose.
TEST(urdf_chain, forward_kinematics_agrees_with_kdl) {
  urdf::Model model = load_model();

  chain<rational> exact_robot;
  ASSERT_TRUE(chain_from_model(model, kRoot, kTip, exact_robot).ok());
  const chain<double> robot = chain_cast<double>(exact_robot);

  KDL::Tree tree;
  ASSERT_TRUE(kdl_parser::treeFromUrdfModel(model, tree));
  KDL::Chain kdl_chain;
  ASSERT_TRUE(tree.getChain(kRoot, kTip, kdl_chain));
  ASSERT_EQ(kdl_chain.getNrOfJoints(), robot.degrees_of_freedom());
  KDL::ChainFkSolverPos_recursive solver(kdl_chain);

  std::mt19937 generator(20260817);
  std::uniform_real_distribution<double> uniform(-2.0, 2.0);

  double worst_position = 0.0;
  double worst_orientation = 0.0;

  for (int trial = 0; trial < 200; ++trial) {
    std::vector<double> values(robot.degrees_of_freedom());
    KDL::JntArray kdl_values(static_cast<unsigned int>(values.size()));
    for (std::size_t i = 0; i < values.size(); ++i) {
      values[i] = uniform(generator);
      kdl_values(static_cast<unsigned int>(i)) = values[i];
    }

    KDL::Frame reference;
    ASSERT_GE(solver.JntToCart(kdl_values, reference), 0);
    const rigid_transform<double> ours = forward_kinematics(robot, values);

    for (std::size_t i = 0; i < 3; ++i) {
      worst_position =
          std::max(worst_position, std::fabs(ours.translation()[i] - reference.p[i]));
    }

    varietas::matrix3<double> theirs;
    for (std::size_t i = 0; i < 3; ++i) {
      for (std::size_t j = 0; j < 3; ++j) {
        theirs(i, j) = reference.M(static_cast<int>(i), static_cast<int>(j));
      }
    }
    worst_orientation = std::max(
        worst_orientation,
        varietas::urdf_import::rotation_distance(ours.rotation(), theirs));
  }

  // The residual is the file's truncated pi propagated through the arm, about a
  // nanometre over a metre of reach, and not an error of the recovery.
  EXPECT_LT(worst_position, 1e-9);
  EXPECT_LT(worst_orientation, 1e-9);
}

// The tangent half-angle substitution, on a real robot.
//
// The planar arm in varietas_kinematics checks the rationalisation against a
// map written a few lines away from it. This checks it against KDL, on seven
// joints whose axes are turned by right angles at every link — the composition
// the closure property has to survive. Substituting t = tan(q / 2) into the
// rational map must reproduce the pose the reference implementation computes
// from the angles directly.
TEST(urdf_chain, the_rational_map_agrees_with_kdl_on_seven_joints) {
  urdf::Model model = load_model();

  chain<rational> exact_robot;
  ASSERT_TRUE(chain_from_model(model, kRoot, kTip, exact_robot).ok());
  const chain<double> robot = chain_cast<double>(exact_robot);
  ASSERT_EQ(robot.degrees_of_freedom(), 7u);

  const auto map =
      varietas::rational_forward_kinematics<7, varietas::grevlex>(robot);

  // Seven revolute joints, so seven factors 1 + t_i^2 and a denominator of
  // degree fourteen, whatever the fixed geometry between them.
  for (std::size_t i = 0; i < 7; ++i) {
    EXPECT_EQ(map.denominator_exponents()[i], 1) << "joint " << i;
  }

  KDL::Tree tree;
  ASSERT_TRUE(kdl_parser::treeFromUrdfModel(model, tree));
  KDL::Chain kdl_chain;
  ASSERT_TRUE(tree.getChain(kRoot, kTip, kdl_chain));
  KDL::ChainFkSolverPos_recursive solver(kdl_chain);

  std::mt19937 generator(20260818);
  std::uniform_real_distribution<double> uniform(-2.0, 2.0);

  double worst_position = 0.0;
  double worst_orientation = 0.0;

  for (int trial = 0; trial < 50; ++trial) {
    std::array<double, 7> point{};
    KDL::JntArray kdl_values(7);
    for (std::size_t i = 0; i < 7; ++i) {
      const double angle = uniform(generator);
      kdl_values(static_cast<unsigned int>(i)) = angle;
      point[i] = varietas::variable_from_angle(angle);
    }

    KDL::Frame reference;
    ASSERT_GE(solver.JntToCart(kdl_values, reference), 0);
    const rigid_transform<double> ours = map.evaluate(point);

    for (std::size_t i = 0; i < 3; ++i) {
      worst_position =
          std::max(worst_position, std::fabs(ours.translation()[i] - reference.p[i]));
    }

    varietas::matrix3<double> theirs;
    for (std::size_t i = 0; i < 3; ++i) {
      for (std::size_t j = 0; j < 3; ++j) {
        theirs(i, j) = reference.M(static_cast<int>(i), static_cast<int>(j));
      }
    }
    worst_orientation = std::max(
        worst_orientation,
        varietas::urdf_import::rotation_distance(ours.rotation(), theirs));
  }

  EXPECT_LT(worst_position, 1e-9);
  EXPECT_LT(worst_orientation, 1e-9);
}

}  // namespace
