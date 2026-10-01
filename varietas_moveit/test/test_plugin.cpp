// The plugin as MoveIt uses it: loaded by name, initialised for a group, and
// asked for poses that the arm's own forward kinematics produced.

#include <cmath>
#include <memory>
#include <random>
#include <vector>

#include <gtest/gtest.h>
#include <moveit/kinematics_base/kinematics_base.h>
#include <moveit/robot_state/robot_state.h>
#include <pluginlib/class_loader.hpp>
#include <rclcpp/rclcpp.hpp>

#include "robot.hpp"
#include "varietas_moveit/varietas_kinematics_plugin.hpp"

namespace {

class Plugin : public ::testing::Test {
 protected:
  void SetUp() override {
    model_ = industrial_robot_model();
    ASSERT_TRUE(model_) << "the industrial arm's URDF or SRDF did not load";
    node_ = std::make_shared<rclcpp::Node>("test_varietas_plugin");
    loader_ = std::make_unique<pluginlib::ClassLoader<kinematics::KinematicsBase>>(
        "moveit_core", "kinematics::KinematicsBase");
    solver_ = loader_->createSharedInstance("varietas_moveit/VarietasKinematicsPlugin");
    ASSERT_TRUE(solver_->initialize(node_, *model_, "manipulator", "base_link", {"flange"}, 0.1));
    group_ = model_->getJointModelGroup("manipulator");
  }

  // A configuration inside the joint ranges, and the pose it puts the flange at.
  std::vector<double> random_configuration(std::mt19937& rng) const {
    moveit::core::RobotState state(model_);
    random_numbers::RandomNumberGenerator generator(rng());
    state.setToRandomPositions(group_, generator);
    std::vector<double> q;
    state.copyJointGroupPositions(group_, q);
    return q;
  }

  geometry_msgs::msg::Pose pose_of(const std::vector<double>& q) const {
    std::vector<geometry_msgs::msg::Pose> poses;
    EXPECT_TRUE(solver_->getPositionFK({"flange"}, q, poses));
    return poses.front();
  }

  static double pose_error(const geometry_msgs::msg::Pose& a, const geometry_msgs::msg::Pose& b) {
    const double dp = std::sqrt(std::pow(a.position.x - b.position.x, 2) +
                                std::pow(a.position.y - b.position.y, 2) +
                                std::pow(a.position.z - b.position.z, 2));
    const double dot = std::abs(a.orientation.w * b.orientation.w + a.orientation.x * b.orientation.x +
                                a.orientation.y * b.orientation.y + a.orientation.z * b.orientation.z);
    return std::max(dp, 1.0 - std::min(1.0, dot));
  }

  bool within_limits(const std::vector<double>& q) const {
    const auto bounds = group_->getActiveJointModelsBounds();
    for (std::size_t i = 0; i < q.size(); ++i) {
      const auto& b = bounds[i]->front();
      if (b.position_bounded_ && (q[i] < b.min_position_ - 1e-9 || q[i] > b.max_position_ + 1e-9)) {
        return false;
      }
    }
    return true;
  }

  moveit::core::RobotModelPtr model_;
  rclcpp::Node::SharedPtr node_;
  std::unique_ptr<pluginlib::ClassLoader<kinematics::KinematicsBase>> loader_;
  kinematics::KinematicsBasePtr solver_;
  const moveit::core::JointModelGroup* group_ = nullptr;
};

TEST_F(Plugin, EveryReachablePoseIsSolvedFromAnySeed) {
  std::mt19937 rng(42u);
  int solved = 0;
  double worst = 0.0;
  const int trials = 2000;
  for (int trial = 0; trial < trials; ++trial) {
    const auto pose = pose_of(random_configuration(rng));
    const auto seed = random_configuration(rng);
    std::vector<double> solution;
    moveit_msgs::msg::MoveItErrorCodes code;
    if (solver_->getPositionIK(pose, seed, solution, code)) {
      ++solved;
      EXPECT_TRUE(within_limits(solution));
      worst = std::max(worst, pose_error(pose_of(solution), pose));
    }
  }
  EXPECT_EQ(solved, trials) << "a pose struck by a configuration inside the ranges went unsolved";
  EXPECT_LT(worst, 1e-9);
}

TEST_F(Plugin, AllConfigurationsAreReturnedAndTheOriginalIsAmongThem) {
  std::mt19937 rng(7u);
  for (int trial = 0; trial < 500; ++trial) {
    const auto original = random_configuration(rng);
    const auto pose = pose_of(original);
    std::vector<std::vector<double>> solutions;
    kinematics::KinematicsResult result;
    ASSERT_TRUE(solver_->getPositionIK({pose}, original, solutions, result,
                                       kinematics::KinematicsQueryOptions()));
    EXPECT_EQ(result.kinematic_error, kinematics::KinematicErrors::OK);
    ASSERT_LE(solutions.size(), 8u);
    bool found = false;
    for (const auto& s : solutions) {
      EXPECT_LT(pose_error(pose_of(s), pose), 1e-9);
      EXPECT_TRUE(within_limits(s));
      double d = 0.0;
      for (std::size_t i = 0; i < s.size(); ++i) {
        d = std::max(d, std::abs(s[i] - original[i]));
      }
      found = found || d < 1e-6;
    }
    EXPECT_TRUE(found) << "trial " << trial;
  }
}

// The seed decides which configuration comes first, not whether one is found:
// seeded at the configuration a pose came from, that configuration is returned.
TEST_F(Plugin, TheConfigurationNearestTheSeedComesFirst) {
  std::mt19937 rng(3u);
  for (int trial = 0; trial < 500; ++trial) {
    const auto original = random_configuration(rng);
    std::vector<double> solution;
    moveit_msgs::msg::MoveItErrorCodes code;
    ASSERT_TRUE(solver_->getPositionIK(pose_of(original), original, solution, code));
    for (std::size_t i = 0; i < original.size(); ++i) {
      EXPECT_NEAR(solution[i], original[i], 1e-6) << "trial " << trial << " joint " << i;
    }
  }
}

// A callback that rejects a configuration, as a collision check does, sends
// the search on to the next one rather than ending it.
TEST_F(Plugin, ARejectedConfigurationGivesWayToTheNext) {
  std::mt19937 rng(9u);
  const auto original = random_configuration(rng);
  const auto pose = pose_of(original);
  int calls = 0;
  const kinematics::KinematicsBase::IKCallbackFn reject_first =
      [&calls](const geometry_msgs::msg::Pose&, const std::vector<double>&,
               moveit_msgs::msg::MoveItErrorCodes& verdict) {
        verdict.val = calls++ == 0 ? moveit_msgs::msg::MoveItErrorCodes::GOAL_IN_COLLISION
                                   : moveit_msgs::msg::MoveItErrorCodes::SUCCESS;
      };
  std::vector<double> solution;
  moveit_msgs::msg::MoveItErrorCodes code;
  ASSERT_TRUE(solver_->searchPositionIK(pose, original, 0.1, solution, reject_first, code));
  EXPECT_EQ(calls, 2);
  double d = 0.0;
  for (std::size_t i = 0; i < solution.size(); ++i) {
    d = std::max(d, std::abs(solution[i] - original[i]));
  }
  EXPECT_GT(d, 1e-6) << "the rejected configuration was returned";
}

TEST_F(Plugin, AnUnreachablePoseIsReportedAtOnce) {
  geometry_msgs::msg::Pose far;
  far.position.x = 5.0;
  far.orientation.w = 1.0;
  std::vector<double> solution;
  moveit_msgs::msg::MoveItErrorCodes code;
  EXPECT_FALSE(solver_->getPositionIK(far, std::vector<double>(6, 0.0), solution, code));
  EXPECT_EQ(code.val, moveit_msgs::msg::MoveItErrorCodes::NO_IK_SOLUTION);
}

}  // namespace

int main(int argc, char** argv) {
  rclcpp::init(argc, argv);
  ::testing::InitGoogleTest(&argc, argv);
  const int result = RUN_ALL_TESTS();
  rclcpp::shutdown();
  return result;
}
