// Six joints and a full pose, through the wrist centre.
//
// The decomposition is checked exactly, on the chain over Q: where the wrist
// centre is, and that an arm whose wrist axes do not meet is refused. The
// generated solvers are checked the way a caller would use them: poses are made
// by sending random configurations through the forward map, and for every pose
// the configuration it came from has to be among those returned, and every
// configuration returned has to reproduce the pose.

#include <algorithm>
#include <cmath>
#include <random>
#include <vector>

#include <gtest/gtest.h>

#include "varietas/ik/runtime.hpp"
#include "varietas/ik/spherical_wrist.hpp"
#include "varietas/kinematics/evaluate.hpp"

#include "arms.hpp"
#include "industrial_ik.hpp"
#include "skewed_ik.hpp"

namespace {

using varietas::rational;

TEST(SphericalWrist, TheWristCentreIsFoundExactly) {
  const auto w = varietas::ik::decompose_spherical_wrist(varietas_test::industrial_six());
  ASSERT_TRUE(w.ok());
  // 25 mm + 420 mm out, 400 + 455 + 35 mm up.
  EXPECT_EQ(w.centre[0], rational(89, 200));
  EXPECT_EQ(w.centre[1], rational(0));
  EXPECT_EQ(w.centre[2], rational(89, 100));
  // 80 mm behind the flange, along the last axis.
  EXPECT_EQ(w.centre_in_tool[0], rational(-2, 25));
  EXPECT_EQ(w.arm.degrees_of_freedom(), 3u);
}

// The commonest wrist rolls, pitches and rolls, so its first and last axes are
// collinear at the zero configuration. That is its singularity, not a defect,
// and must not be mistaken for a wrist without a centre.
TEST(SphericalWrist, CollinearFirstAndLastWristAxesAreAccepted) {
  EXPECT_TRUE(varietas::ik::decompose_spherical_wrist(varietas_test::industrial_six()).ok());
}

TEST(SphericalWrist, AWristWhoseAxesDoNotMeetIsRefused) {
  const auto w = varietas::ik::decompose_spherical_wrist(varietas_test::broken_wrist_six());
  EXPECT_EQ(w.status, varietas::ik::wrist_status::wrist_axes_do_not_meet);
}

TEST(SphericalWrist, OnlySixRevoluteJointsAreDecomposed) {
  const auto w = varietas::ik::decompose_spherical_wrist(varietas_test::anthropomorphic_three_link());
  EXPECT_EQ(w.status, varietas::ik::wrist_status::wrong_joints);
}

double wrapped(double a) {
  a = std::fmod(a + M_PI, 2.0 * M_PI);
  return a < 0 ? a + M_PI : a - M_PI;
}

template <class Solver>
void every_pose_is_solved(const varietas::chain<rational>& exact) {
  const auto arm = varietas::chain_cast<double>(exact);
  std::mt19937 rng(20260930u);
  std::uniform_real_distribution<double> angle(-M_PI, M_PI);
  double worst_position = 0.0;
  double worst_rotation = 0.0;
  for (int trial = 0; trial < 2000; ++trial) {
    std::vector<double> q(6);
    for (double& v : q) {
      v = angle(rng);
    }
    const auto target = varietas::forward_kinematics(arm, q);
    double position[3];
    double rotation[9];
    for (int i = 0; i < 3; ++i) {
      position[i] = target.translation()[i];
      for (int j = 0; j < 3; ++j) {
        rotation[3 * i + j] = target.rotation()(i, j);
      }
    }
    double out[Solver::max_configurations * 6];
    const int found = Solver::solve(position, rotation, out,
                                    static_cast<int>(Solver::max_configurations));
    ASSERT_GE(found, 1) << "trial " << trial;
    ASSERT_LE(found, 8);

    bool original = false;
    for (int k = 0; k < found; ++k) {
      const std::vector<double> s(out + 6 * k, out + 6 * k + 6);
      const auto reached = varietas::forward_kinematics(arm, s);
      for (int i = 0; i < 3; ++i) {
        worst_position = std::max(worst_position,
                                  std::abs(reached.translation()[i] - position[i]));
        for (int j = 0; j < 3; ++j) {
          worst_rotation = std::max(worst_rotation,
                                    std::abs(reached.rotation()(i, j) - rotation[3 * i + j]));
        }
      }
      double distance = 0.0;
      for (int i = 0; i < 6; ++i) {
        distance = std::max(distance, std::abs(wrapped(s[i] - q[i])));
      }
      original = original || distance < 1e-6;
    }
    EXPECT_TRUE(original) << "the configuration the pose came from is missing, trial " << trial;
  }
  EXPECT_LT(worst_position, 1e-9);
  EXPECT_LT(worst_rotation, 1e-9);
}

TEST(SphericalWrist, TheIndustrialArmSolvesEveryPoseItCanStrike) {
  every_pose_is_solved<varietas_generated::industrial_ik>(varietas_test::industrial_six());
}

// A wrist whose consecutive axes are not perpendicular reaches only some
// orientations, so some arm configurations have no completion and the count
// varies; what may not vary is that the pose's own configuration is found.
TEST(SphericalWrist, ASkewedWristOnAnArmThatDoesNotDecoupleSolvesEveryPose) {
  every_pose_is_solved<varietas_generated::skewed_ik>(varietas_test::skewed_wrist_six());
}

// The ranges of the description are respected: every configuration
// solve_within_limits returns lies inside them, and the configuration a pose
// came from is among those returned whenever it was inside them itself.
TEST(SphericalWrist, SolutionsAreFittedIntoTheJointRanges) {
  using solver = varietas_generated::industrial_ik;
  const auto exact = varietas_test::industrial_six_limited();
  const auto arm = varietas::chain_cast<double>(exact);
  const auto limits = varietas::ik::joint_limits::of(exact);
  std::mt19937 rng(11u);
  std::uniform_real_distribution<double> angle(-M_PI, M_PI);
  int fewer = 0;
  for (int trial = 0; trial < 2000; ++trial) {
    std::vector<double> q(6);
    for (double& v : q) {
      v = angle(rng);
    }
    const bool inside = limits.fit(q.data());
    const auto target = varietas::forward_kinematics(arm, q);
    double position[3];
    double rotation[9];
    for (int i = 0; i < 3; ++i) {
      position[i] = target.translation()[i];
      for (int j = 0; j < 3; ++j) {
        rotation[3 * i + j] = target.rotation()(i, j);
      }
    }
    double all[solver::max_configurations * 6];
    double kept[solver::max_configurations * 6];
    const int found = solver::solve(position, rotation, all, 16);
    const int within = solver::solve_within_limits(position, rotation, kept, 16);
    ASSERT_GE(within, 0);
    ASSERT_LE(within, found);
    fewer += within < found;
    bool original = false;
    for (int k = 0; k < within; ++k) {
      for (int i = 0; i < 6; ++i) {
        EXPECT_GE(kept[6 * k + i], limits.lower[i] - 1e-9);
        EXPECT_LE(kept[6 * k + i], limits.upper[i] + 1e-9);
      }
      double distance = 0.0;
      for (int i = 0; i < 6; ++i) {
        distance = std::max(distance, std::abs(kept[6 * k + i] - q[i]));
      }
      original = original || distance < 1e-6;
    }
    if (inside) {
      EXPECT_TRUE(original) << "trial " << trial;
    }
  }
  // The ranges are narrower than a turn, so some configurations must go.
  EXPECT_GT(fewer, 0);
}

// The solver evaluated in process returns what the generated header returns,
// configuration for configuration, on both routes the arm can take.
template <class Solver>
void runtime_agrees(const varietas::chain<rational>& exact) {
  std::string why;
  const auto runtime = varietas::ik::build_wrist_solver(exact, &why);
  ASSERT_TRUE(runtime.has_value()) << why;
  const auto arm = varietas::chain_cast<double>(exact);
  std::mt19937 rng(5u);
  std::uniform_real_distribution<double> angle(-M_PI, M_PI);
  for (int trial = 0; trial < 500; ++trial) {
    std::vector<double> q(6);
    for (double& v : q) {
      v = angle(rng);
    }
    const auto target = varietas::forward_kinematics(arm, q);
    double position[3];
    double rotation[9];
    for (int i = 0; i < 3; ++i) {
      position[i] = target.translation()[i];
      for (int j = 0; j < 3; ++j) {
        rotation[3 * i + j] = target.rotation()(i, j);
      }
    }
    double a[16 * 6];
    double b[16 * 6];
    const int na = runtime->solve(position, rotation, a, 16);
    const int nb = Solver::solve(position, rotation, b, 16);
    ASSERT_EQ(na, nb) << "trial " << trial;
    for (int k = 0; k < na; ++k) {
      double best = 1e9;
      for (int l = 0; l < nb; ++l) {
        double d = 0.0;
        for (int i = 0; i < 6; ++i) {
          d = std::max(d, std::abs(std::remainder(a[6 * k + i] - b[6 * l + i], 2.0 * M_PI)));
        }
        best = std::min(best, d);
      }
      EXPECT_LT(best, 1e-8) << "trial " << trial;
    }
  }
}

TEST(SphericalWrist, TheRuntimeSolverMatchesTheGeneratedHeaders) {
  runtime_agrees<varietas_generated::industrial_ik>(varietas_test::industrial_six());
  runtime_agrees<varietas_generated::skewed_ik>(varietas_test::skewed_wrist_six());
}

}  // namespace
