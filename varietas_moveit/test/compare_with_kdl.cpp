// varietas beside MoveIt's default kinematics solver, on the industrial arm.
//
// Both plugins are loaded by name and initialised for the same group, and
// asked for the same poses from the same seeds. A pose is made by sending a
// random configuration inside the joint ranges through the forward map, so
// every pose is reachable. Reported: how often each finds a configuration,
// what a call costs, and how far the answer lands from the pose.
//
//   ros2 run varietas_moveit compare_with_kdl     (built with the tests)

#include <chrono>
#include <cmath>
#include <cstdio>
#include <memory>
#include <random>
#include <vector>

#include <moveit/kinematics_base/kinematics_base.h>
#include <moveit/robot_state/robot_state.h>
#include <pluginlib/class_loader.hpp>
#include <rclcpp/rclcpp.hpp>

#include "robot.hpp"

namespace {

struct tally {
  int solved = 0;
  std::vector<double> microseconds;
  double worst = 0.0;
};

double median(std::vector<double> v) {
  std::sort(v.begin(), v.end());
  return v.empty() ? 0.0 : v[v.size() / 2];
}

double percentile(std::vector<double> v, double p) {
  std::sort(v.begin(), v.end());
  return v.empty() ? 0.0 : v[static_cast<std::size_t>(p * (v.size() - 1))];
}

}  // namespace

int main(int argc, char** argv) {
  rclcpp::init(argc, argv);
  auto model = industrial_robot_model();
  auto node = std::make_shared<rclcpp::Node>("compare_with_kdl");
  pluginlib::ClassLoader<kinematics::KinematicsBase> loader("moveit_core",
                                                           "kinematics::KinematicsBase");
  auto varietas = loader.createSharedInstance("varietas_moveit/VarietasKinematicsPlugin");
  auto kdl = loader.createSharedInstance("kdl_kinematics_plugin/KDLKinematicsPlugin");
  if (!varietas->initialize(node, *model, "manipulator", "base_link", {"flange"}, 0.1) ||
      !kdl->initialize(node, *model, "manipulator", "base_link", {"flange"}, 0.1)) {
    std::fprintf(stderr, "a plugin did not initialise\n");
    return 1;
  }
  const auto* group = model->getJointModelGroup("manipulator");

  std::mt19937 rng(20260930u);
  random_numbers::RandomNumberGenerator generator(1234u);
  moveit::core::RobotState state(model);
  const int trials = 2000;
  const double timeout = 0.05;
  tally v, k;
  for (int t = 0; t < trials; ++t) {
    state.setToRandomPositions(group, generator);
    std::vector<double> target_q;
    state.copyJointGroupPositions(group, target_q);
    std::vector<geometry_msgs::msg::Pose> poses;
    varietas->getPositionFK({"flange"}, target_q, poses);
    state.setToRandomPositions(group, generator);
    std::vector<double> seed;
    state.copyJointGroupPositions(group, seed);

    for (int which = 0; which < 2; ++which) {
      auto& solver = which == 0 ? varietas : kdl;
      tally& out = which == 0 ? v : k;
      std::vector<double> solution;
      moveit_msgs::msg::MoveItErrorCodes code;
      const auto start = std::chrono::steady_clock::now();
      const bool ok = solver->searchPositionIK(poses.front(), seed, timeout, solution, code);
      out.microseconds.push_back(
          std::chrono::duration<double, std::micro>(std::chrono::steady_clock::now() - start).count());
      if (ok) {
        ++out.solved;
        std::vector<geometry_msgs::msg::Pose> reached;
        varietas->getPositionFK({"flange"}, solution, reached);
        const auto& a = reached.front().position;
        const auto& b = poses.front().position;
        out.worst = std::max(out.worst, std::sqrt(std::pow(a.x - b.x, 2) + std::pow(a.y - b.y, 2) +
                                                  std::pow(a.z - b.z, 2)));
      }
    }
  }
  std::printf("%d reachable poses, random seeds, timeout %.0f ms for KDL\n\n", trials,
              timeout * 1000.0);
  std::printf("%-10s %10s %14s %14s %16s\n", "", "solved", "median", "99th centile",
              "worst position");
  for (int which = 0; which < 2; ++which) {
    const tally& t = which == 0 ? v : k;
    std::printf("%-10s %9.1f%% %11.1f us %11.1f us %14.1e m\n", which == 0 ? "varietas" : "KDL",
                100.0 * t.solved / trials, median(t.microseconds), percentile(t.microseconds, 0.99),
                t.worst);
  }
  rclcpp::shutdown();
  return 0;
}
