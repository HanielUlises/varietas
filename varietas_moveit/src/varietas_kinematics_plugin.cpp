#include "varietas_moveit/varietas_kinematics_plugin.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>

#include <Eigen/Geometry>
#include <moveit/robot_model/robot_model.h>
#include <moveit/robot_state/robot_state.h>
#include <pluginlib/class_list_macros.hpp>

#include "varietas/urdf/urdf_chain.hpp"

namespace varietas_moveit {

namespace {

// Not LOGGER: KinematicsBase has a static member of that name, which would
// shadow this one inside the member functions.
const rclcpp::Logger kLogger = rclcpp::get_logger("varietas_moveit");

constexpr double kTurn = 6.28318530717958647692;

std::string without_leading_slash(const std::string& frame) {
  return !frame.empty() && frame.front() == '/' ? frame.substr(1) : frame;
}

}  // namespace

bool VarietasKinematicsPlugin::initialize(const rclcpp::Node::SharedPtr& node,
                                          const moveit::core::RobotModel& robot_model,
                                          const std::string& group_name,
                                          const std::string& base_frame,
                                          const std::vector<std::string>& tip_frames,
                                          double search_discretization) {
  node_ = node;
  storeValues(robot_model, group_name, base_frame, tip_frames, search_discretization);

  group_ = robot_model.getJointModelGroup(group_name);
  if (group_ == nullptr) {
    RCLCPP_ERROR(kLogger, "group '%s' does not exist", group_name.c_str());
    return false;
  }
  if (tip_frames.size() != 1) {
    RCLCPP_ERROR(kLogger, "group '%s': one tip frame is supported, %zu were given",
                 group_name.c_str(), tip_frames.size());
    return false;
  }
  const auto& active = group_->getActiveJointModels();
  if (active.size() != 6) {
    RCLCPP_ERROR(kLogger, "group '%s' has %zu active joints; this plugin solves six-joint arms "
                         "whose last three axes meet",
                 group_name.c_str(), active.size());
    return false;
  }

  // The chain between the group's base and tip, from the URDF MoveIt holds,
  // recovered exactly.
  const auto urdf = robot_model.getURDF();
  if (!urdf) {
    RCLCPP_ERROR(kLogger, "the robot model carries no URDF");
    return false;
  }
  const std::string base = without_leading_slash(base_frame);
  const std::string tip = without_leading_slash(tip_frames.front());
  varietas::chain<varietas::rational> exact;
  const auto imported = varietas::urdf_import::chain_from_model(*urdf, base, tip, exact);
  if (!imported.ok()) {
    RCLCPP_ERROR(kLogger, "the chain from '%s' to '%s' could not be recovered exactly: %s (%s)",
                 base.c_str(), tip.c_str(), varietas::urdf_import::to_string(imported.status),
                 imported.detail.c_str());
    return false;
  }

  const auto start = std::chrono::steady_clock::now();
  std::string why;
  solver_ = varietas::ik::build_wrist_solver(exact, &why);
  if (!solver_) {
    RCLCPP_ERROR(kLogger, "group '%s' was refused: %s", group_name.c_str(), why.c_str());
    return false;
  }
  const double seconds =
      std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();

  // The chain's joints, in its order, matched to the group's by name, with the
  // ranges MoveIt knows, which may override the URDF's.
  const auto& group_names = group_->getActiveJointModelNames();
  const auto bounds = group_->getActiveJointModelsBounds();
  limits_ = varietas::ik::joint_limits{};
  group_index_.clear();
  // Held in a variable: iterating fold_fixed_joints().joints() directly would
  // walk the joints of a temporary chain destroyed before the loop begins.
  const varietas::chain<varietas::rational> folded = exact.fold_fixed_joints();
  for (const auto& j : folded.joints()) {
    if (!j.is_actuated()) {
      continue;
    }
    const auto it = std::find(group_names.begin(), group_names.end(), j.name);
    if (it == group_names.end()) {
      RCLCPP_ERROR(kLogger, "joint '%s' of the chain is not an active joint of group '%s'",
                   j.name.c_str(), group_name.c_str());
      solver_.reset();
      return false;
    }
    const auto index = static_cast<std::size_t>(it - group_names.begin());
    group_index_.push_back(index);
    const auto& variable = bounds[index]->front();
    limits_.limited.push_back(variable.position_bounded_);
    limits_.lower.push_back(variable.min_position_);
    limits_.upper.push_back(variable.max_position_);
  }

  joint_names_ = group_names;
  link_names_ = {tip};
  RCLCPP_INFO(kLogger,
              "group '%s': split at the wrist, arm %s in %.3f s; up to %zu configurations a pose",
              group_name.c_str(), solver_->arm().is_decoupled() ? "decoupled" : "reconstructed",
              seconds, solver_->max_configurations());
  return true;
}

std::vector<std::vector<double>> VarietasKinematicsPlugin::solve_all(
    const geometry_msgs::msg::Pose& pose, const std::vector<double>& seed) const {
  std::vector<std::vector<double>> out;
  if (!solver_) {
    return out;
  }
  const Eigen::Quaterniond q(pose.orientation.w, pose.orientation.x, pose.orientation.y,
                             pose.orientation.z);
  const Eigen::Matrix3d r = q.normalized().toRotationMatrix();
  const double position[3] = {pose.position.x, pose.position.y, pose.position.z};
  double rotation[9];
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      rotation[3 * i + j] = r(i, j);
    }
  }

  const std::size_t capacity = solver_->max_configurations();
  std::vector<double> raw(capacity * 6);
  const int found = solver_->solve(position, rotation, raw.data(), static_cast<int>(capacity));
  for (int k = 0; k < found; ++k) {
    std::vector<double> configuration(6, 0.0);
    bool reachable = true;
    for (std::size_t c = 0; c < 6 && reachable; ++c) {
      const std::size_t g = group_index_[c];
      const double angle = raw[static_cast<std::size_t>(k) * 6 + c];
      const double near = g < seed.size() ? seed[g] : 0.0;
      // Of the angle's representatives a whole turn apart, the one inside the
      // range nearest the seed, so that a joint that can reach an angle two
      // ways is not sent round a full turn to get there.
      double best = std::numeric_limits<double>::quiet_NaN();
      double best_distance = std::numeric_limits<double>::infinity();
      const double centre = angle + kTurn * std::round((near - angle) / kTurn);
      for (int turns = -2; turns <= 2; ++turns) {
        const double v = centre + kTurn * turns;
        if (limits_.limited[c] && (v < limits_.lower[c] - 1e-9 || v > limits_.upper[c] + 1e-9)) {
          continue;
        }
        if (std::abs(v - near) < best_distance) {
          best_distance = std::abs(v - near);
          best = v;
        }
      }
      if (std::isnan(best)) {
        reachable = false;
        break;
      }
      configuration[g] = best;
    }
    if (reachable) {
      out.push_back(std::move(configuration));
    }
  }

  const auto distance = [&](const std::vector<double>& a) {
    double d = 0.0;
    for (std::size_t i = 0; i < a.size() && i < seed.size(); ++i) {
      d += (a[i] - seed[i]) * (a[i] - seed[i]);
    }
    return d;
  };
  std::stable_sort(out.begin(), out.end(), [&](const std::vector<double>& a,
                                               const std::vector<double>& b) {
    return distance(a) < distance(b);
  });
  return out;
}

bool VarietasKinematicsPlugin::search(const geometry_msgs::msg::Pose& ik_pose,
                                      const std::vector<double>& ik_seed_state,
                                      const std::vector<double>* consistency_limits,
                                      std::vector<double>& solution,
                                      const IKCallbackFn* solution_callback,
                                      moveit_msgs::msg::MoveItErrorCodes& error_code) const {
  if (!solver_) {
    error_code.val = moveit_msgs::msg::MoveItErrorCodes::FAILURE;
    return false;
  }
  for (const auto& candidate : solve_all(ik_pose, ik_seed_state)) {
    if (consistency_limits != nullptr && !consistency_limits->empty()) {
      bool consistent = true;
      for (std::size_t i = 0; i < candidate.size() && i < consistency_limits->size(); ++i) {
        if (std::abs(candidate[i] - ik_seed_state[i]) > (*consistency_limits)[i]) {
          consistent = false;
          break;
        }
      }
      if (!consistent) {
        continue;
      }
    }
    if (solution_callback != nullptr && *solution_callback) {
      moveit_msgs::msg::MoveItErrorCodes verdict;
      (*solution_callback)(ik_pose, candidate, verdict);
      if (verdict.val != moveit_msgs::msg::MoveItErrorCodes::SUCCESS) {
        continue;  // rejected, say for a collision; the next configuration may do
      }
    }
    solution = candidate;
    error_code.val = moveit_msgs::msg::MoveItErrorCodes::SUCCESS;
    return true;
  }
  error_code.val = moveit_msgs::msg::MoveItErrorCodes::NO_IK_SOLUTION;
  return false;
}

bool VarietasKinematicsPlugin::getPositionIK(const geometry_msgs::msg::Pose& ik_pose,
                                             const std::vector<double>& ik_seed_state,
                                             std::vector<double>& solution,
                                             moveit_msgs::msg::MoveItErrorCodes& error_code,
                                             const kinematics::KinematicsQueryOptions&) const {
  return search(ik_pose, ik_seed_state, nullptr, solution, nullptr, error_code);
}

bool VarietasKinematicsPlugin::getPositionIK(const std::vector<geometry_msgs::msg::Pose>& ik_poses,
                                             const std::vector<double>& ik_seed_state,
                                             std::vector<std::vector<double>>& solutions,
                                             kinematics::KinematicsResult& result,
                                             const kinematics::KinematicsQueryOptions&) const {
  solutions.clear();
  if (ik_poses.empty()) {
    result.kinematic_error = kinematics::KinematicErrors::EMPTY_TIP_POSES;
    return false;
  }
  if (ik_poses.size() != 1) {
    result.kinematic_error = kinematics::KinematicErrors::MULTIPLE_TIPS_NOT_SUPPORTED;
    return false;
  }
  solutions = solve_all(ik_poses.front(), ik_seed_state);
  result.kinematic_error = solutions.empty() ? kinematics::KinematicErrors::NO_SOLUTION
                                             : kinematics::KinematicErrors::OK;
  // Every configuration there is was examined, so the search was complete.
  result.solution_percentage = 1.0;
  return !solutions.empty();
}

bool VarietasKinematicsPlugin::searchPositionIK(const geometry_msgs::msg::Pose& ik_pose,
                                                const std::vector<double>& ik_seed_state,
                                                double, std::vector<double>& solution,
                                                moveit_msgs::msg::MoveItErrorCodes& error_code,
                                                const kinematics::KinematicsQueryOptions&) const {
  return search(ik_pose, ik_seed_state, nullptr, solution, nullptr, error_code);
}

bool VarietasKinematicsPlugin::searchPositionIK(const geometry_msgs::msg::Pose& ik_pose,
                                                const std::vector<double>& ik_seed_state,
                                                double,
                                                const std::vector<double>& consistency_limits,
                                                std::vector<double>& solution,
                                                moveit_msgs::msg::MoveItErrorCodes& error_code,
                                                const kinematics::KinematicsQueryOptions&) const {
  return search(ik_pose, ik_seed_state, &consistency_limits, solution, nullptr, error_code);
}

bool VarietasKinematicsPlugin::searchPositionIK(const geometry_msgs::msg::Pose& ik_pose,
                                                const std::vector<double>& ik_seed_state,
                                                double, std::vector<double>& solution,
                                                const IKCallbackFn& solution_callback,
                                                moveit_msgs::msg::MoveItErrorCodes& error_code,
                                                const kinematics::KinematicsQueryOptions&) const {
  return search(ik_pose, ik_seed_state, nullptr, solution, &solution_callback, error_code);
}

bool VarietasKinematicsPlugin::searchPositionIK(const geometry_msgs::msg::Pose& ik_pose,
                                                const std::vector<double>& ik_seed_state,
                                                double,
                                                const std::vector<double>& consistency_limits,
                                                std::vector<double>& solution,
                                                const IKCallbackFn& solution_callback,
                                                moveit_msgs::msg::MoveItErrorCodes& error_code,
                                                const kinematics::KinematicsQueryOptions&) const {
  return search(ik_pose, ik_seed_state, &consistency_limits, solution, &solution_callback,
                error_code);
}

bool VarietasKinematicsPlugin::getPositionFK(const std::vector<std::string>& link_names,
                                             const std::vector<double>& joint_angles,
                                             std::vector<geometry_msgs::msg::Pose>& poses) const {
  if (!group_ || joint_angles.size() != joint_names_.size()) {
    return false;
  }
  moveit::core::RobotState state(robot_model_);
  state.setToDefaultValues();
  state.setJointGroupPositions(group_, joint_angles);
  state.update();
  const Eigen::Isometry3d base_inverse = state.getGlobalLinkTransform(base_frame_).inverse();
  poses.clear();
  for (const auto& name : link_names) {
    const Eigen::Isometry3d t = base_inverse * state.getGlobalLinkTransform(name);
    const Eigen::Quaterniond q(t.rotation());
    geometry_msgs::msg::Pose pose;
    pose.position.x = t.translation().x();
    pose.position.y = t.translation().y();
    pose.position.z = t.translation().z();
    pose.orientation.w = q.w();
    pose.orientation.x = q.x();
    pose.orientation.y = q.y();
    pose.orientation.z = q.z();
    poses.push_back(pose);
  }
  return true;
}

}  // namespace varietas_moveit

PLUGINLIB_EXPORT_CLASS(varietas_moveit::VarietasKinematicsPlugin, kinematics::KinematicsBase)
