// Inverse kinematics at one pose, orientation included.
//
// This is the capability the emitter cannot have. A generated solver adjoins
// the target to the coefficient field, and orientation cannot be adjoined at
// any workable size: a general pose needs six parameters and two is the
// working limit. Solving one pose at a time gives up the generated code and
// gets the other nine equations in return.
//
// The targets here are built by evaluating the exact forward map at rational
// half-angle values, rather than by rounding sines and cosines. That is not
// tidiness: a target assembled from doubles has denominators near 2^52, and
// Buchberger over coefficients that size does not finish.

#include <array>
#include <cmath>
#include <cstddef>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include <urdf/model.h>

#include "varietas/codegen/rational.hpp"
#include "varietas/core/order/grevlex.hpp"
#include "varietas/core/quotient/quotient_basis.hpp"
#include "varietas/core/solve/spectral.hpp"
#include "varietas/kinematics/evaluate.hpp"
#include "varietas/kinematics/rationalize.hpp"
#include "varietas/urdf/urdf_chain.hpp"

namespace {

using varietas::grevlex;
using varietas::rational;

varietas::chain<rational> load(const std::string& file) {
  urdf::Model model;
  EXPECT_TRUE(model.initFile(std::string(VARIETAS_URDF_TEST_DATA) + "/" + file));

  varietas::chain<rational> robot;
  const auto report = varietas::urdf_import::chain_from_model(
      model, model.getRoot()->name, varietas::urdf_import::sole_tip_link(model), robot);
  EXPECT_TRUE(report.ok());
  return robot.fold_fixed_joints();
}

// The pose the arm actually reaches at the given half-angle values, exactly.
template <std::size_t N>
varietas::rigid_transform<rational> pose_at(const varietas::chain<rational>& robot,
                                            const std::array<rational, N>& t) {
  const auto map = varietas::rational_forward_kinematics<N, grevlex>(robot);
  const rational denominator = map.denominator().evaluate(t);

  varietas::matrix3<rational> rotation;
  for (std::size_t i = 0; i < 3; ++i) {
    for (std::size_t j = 0; j < 3; ++j) {
      rotation(i, j) = map.rotation(i, j).evaluate(t) / denominator;
    }
  }
  varietas::vector3<rational> translation(map.translation(0).evaluate(t) / denominator,
                                          map.translation(1).evaluate(t) / denominator,
                                          map.translation(2).evaluate(t) / denominator);
  return varietas::rigid_transform<rational>(rotation, translation);
}

TEST(PoseIk, AFullPoseIsSolvedAndTheConfigurationRecovered) {
  const auto robot = load("anthropomorphic_3r.urdf");
  ASSERT_EQ(robot.degrees_of_freedom(), 3u);

  // t = (0, 0, 1) is q = (0, 0, pi/2): the elbow bent square.
  const std::array<rational, 3> t{rational(0), rational(0), rational(1)};
  const auto target = pose_at<3>(robot, t);

  const auto basis = varietas::pose_ideal_generators<3, grevlex>(robot, target);
  const auto quotient = varietas::standard_monomials(basis);
  ASSERT_TRUE(quotient.is_zero_dimensional);

  // Three joints against the twelve equations of a full pose: the orientation
  // pins down what the position alone would leave as four postures.
  EXPECT_EQ(quotient.dimension(), 1u);

  const auto solutions = varietas::solve_zero_dimensional(basis);
  ASSERT_TRUE(solutions.ok());
  const auto real = solutions.real_points();
  ASSERT_EQ(real.size(), 1u);

  EXPECT_NEAR(varietas::angle_from_variable(real[0][0]), 0.0, 1e-9);
  EXPECT_NEAR(varietas::angle_from_variable(real[0][1]), 0.0, 1e-9);
  EXPECT_NEAR(varietas::angle_from_variable(real[0][2]), M_PI / 2.0, 1e-9);
}

TEST(PoseIk, AnOrientationTheArmCannotProduceIsRefusedAsUnreachable) {
  const auto robot = load("anthropomorphic_3r.urdf");

  // The arm's orientation is always a yaw followed by a pitch, so any roll at
  // all is outside its reach. Rolling the reachable pose by a right angle,
  // exactly, so that nothing is blamed on rounding, must leave no solution.
  const std::array<rational, 3> t{rational(0), rational(0), rational(1)};
  const auto reachable = pose_at<3>(robot, t);

  varietas::matrix3<rational> roll;  // a quarter turn about x, exactly
  roll(0, 0) = rational(1); roll(0, 1) = rational(0);  roll(0, 2) = rational(0);
  roll(1, 0) = rational(0); roll(1, 1) = rational(0);  roll(1, 2) = rational(-1);
  roll(2, 0) = rational(0); roll(2, 1) = rational(1);  roll(2, 2) = rational(0);

  const varietas::rigid_transform<rational> rolled(roll * reachable.rotation(),
                                                   reachable.translation());
  const auto basis = varietas::pose_ideal_generators<3, grevlex>(robot, rolled);
  const auto quotient = varietas::standard_monomials(basis);

  // The unit ideal: no configuration, rather than a configuration that misses.
  EXPECT_EQ(quotient.dimension(), 0u);
}

TEST(PoseIk, PositionAloneLeavesTheFourPosturesTheDecoupledSolverAlsoFinds) {
  const auto robot = load("anthropomorphic_3r.urdf");

  const std::array<rational, 3> t{rational(1, 4), rational(1, 3), rational(1, 2)};
  const auto target = pose_at<3>(robot, t);

  const auto basis =
      varietas::position_ideal_generators<3, grevlex>(robot, target.translation());
  const auto quotient = varietas::standard_monomials(basis);
  ASSERT_TRUE(quotient.is_zero_dimensional);

  // The same count varietas_ik's decoupled solve reports for this arm, reached
  // here without decoupling anything.
  EXPECT_EQ(quotient.dimension(), 4u);

  const auto solutions = varietas::solve_zero_dimensional(basis);
  ASSERT_TRUE(solutions.ok());
  const auto real = solutions.real_points();
  EXPECT_EQ(real.size(), 4u);

  // Every one of them puts the tool where it was asked to go, checked against
  // the numerical forward map rather than against the algebra that produced it.
  const auto numeric = varietas::chain_cast<double>(robot);
  for (const auto& point : real) {
    const std::vector<double> angles{varietas::angle_from_variable(point[0]),
                                     varietas::angle_from_variable(point[1]),
                                     varietas::angle_from_variable(point[2])};
    const auto frame = varietas::forward_kinematics(numeric, angles);
    for (std::size_t i = 0; i < 3; ++i) {
      EXPECT_NEAR(frame.translation()[i],
                  varietas::coefficient_traits<rational>::to_double(target.translation()[i]),
                  1e-9);
    }
  }
}

}  // namespace
