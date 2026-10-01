#ifndef VARIETAS_MOVEIT_VARIETAS_KINEMATICS_PLUGIN_HPP
#define VARIETAS_MOVEIT_VARIETAS_KINEMATICS_PLUGIN_HPP

#include <optional>
#include <string>
#include <vector>

#include <moveit/kinematics_base/kinematics_base.h>
#include <moveit/robot_model/joint_model_group.h>

#include "varietas/ik/runtime.hpp"

namespace varietas_moveit {

// A MoveIt kinematics plugin for six-joint arms whose last three axes meet.
//
// Everything expensive happens in initialize(), once, when MoveIt loads the
// plugin for a group: the group's chain is read from the URDF MoveIt already
// holds and recovered exactly, split at its wrist centre, and the arm that
// places the centre is solved, by the decoupling in milliseconds where the
// base sweeps a plane and by reconstruction in seconds where it does not. A
// query is then the solver evaluated at the pose: every configuration that
// reaches it, fitted by whole turns into the joint ranges MoveIt knows, which
// may come from joint_limits.yaml rather than the URDF, ordered by distance
// from the seed.
//
// So a query never searches and never depends on its seed for whether it
// succeeds, only for which answer comes first. A pose with no configuration
// inside the ranges is reported as having none, at once, rather than after a
// timeout.
class VarietasKinematicsPlugin : public kinematics::KinematicsBase {
 public:
  bool initialize(const rclcpp::Node::SharedPtr& node, const moveit::core::RobotModel& robot_model,
                  const std::string& group_name, const std::string& base_frame,
                  const std::vector<std::string>& tip_frames, double search_discretization) override;

  bool getPositionIK(const geometry_msgs::msg::Pose& ik_pose, const std::vector<double>& ik_seed_state,
                     std::vector<double>& solution, moveit_msgs::msg::MoveItErrorCodes& error_code,
                     const kinematics::KinematicsQueryOptions& options =
                         kinematics::KinematicsQueryOptions()) const override;

  // Every configuration, which is what this plugin is for.
  bool getPositionIK(const std::vector<geometry_msgs::msg::Pose>& ik_poses,
                     const std::vector<double>& ik_seed_state,
                     std::vector<std::vector<double>>& solutions,
                     kinematics::KinematicsResult& result,
                     const kinematics::KinematicsQueryOptions& options) const override;

  bool searchPositionIK(const geometry_msgs::msg::Pose& ik_pose,
                        const std::vector<double>& ik_seed_state, double timeout,
                        std::vector<double>& solution, moveit_msgs::msg::MoveItErrorCodes& error_code,
                        const kinematics::KinematicsQueryOptions& options =
                            kinematics::KinematicsQueryOptions()) const override;

  bool searchPositionIK(const geometry_msgs::msg::Pose& ik_pose,
                        const std::vector<double>& ik_seed_state, double timeout,
                        const std::vector<double>& consistency_limits,
                        std::vector<double>& solution, moveit_msgs::msg::MoveItErrorCodes& error_code,
                        const kinematics::KinematicsQueryOptions& options =
                            kinematics::KinematicsQueryOptions()) const override;

  bool searchPositionIK(const geometry_msgs::msg::Pose& ik_pose,
                        const std::vector<double>& ik_seed_state, double timeout,
                        std::vector<double>& solution, const IKCallbackFn& solution_callback,
                        moveit_msgs::msg::MoveItErrorCodes& error_code,
                        const kinematics::KinematicsQueryOptions& options =
                            kinematics::KinematicsQueryOptions()) const override;

  bool searchPositionIK(const geometry_msgs::msg::Pose& ik_pose,
                        const std::vector<double>& ik_seed_state, double timeout,
                        const std::vector<double>& consistency_limits,
                        std::vector<double>& solution, const IKCallbackFn& solution_callback,
                        moveit_msgs::msg::MoveItErrorCodes& error_code,
                        const kinematics::KinematicsQueryOptions& options =
                            kinematics::KinematicsQueryOptions()) const override;

  bool getPositionFK(const std::vector<std::string>& link_names,
                     const std::vector<double>& joint_angles,
                     std::vector<geometry_msgs::msg::Pose>& poses) const override;

  const std::vector<std::string>& getJointNames() const override { return joint_names_; }
  const std::vector<std::string>& getLinkNames() const override { return link_names_; }

  // Every configuration that reaches `pose`, in the group's joint order,
  // fitted into the ranges and nearest `seed` first. Public so that it can be
  // called without the MoveIt plumbing around it.
  std::vector<std::vector<double>> solve_all(const geometry_msgs::msg::Pose& pose,
                                             const std::vector<double>& seed) const;

  // Whether the arm decoupled, for logs and tests.
  bool decoupled() const { return solver_ && solver_->arm().is_decoupled(); }

 private:
  bool search(const geometry_msgs::msg::Pose& ik_pose, const std::vector<double>& ik_seed_state,
              const std::vector<double>* consistency_limits, std::vector<double>& solution,
              const IKCallbackFn* solution_callback,
              moveit_msgs::msg::MoveItErrorCodes& error_code) const;

  const moveit::core::JointModelGroup* group_ = nullptr;
  std::optional<varietas::ik::runtime_wrist_solver> solver_;
  varietas::ik::joint_limits limits_;  // in chain order
  // group_index_[k] is the position in the group's joint vector of the chain's
  // k-th joint.
  std::vector<std::size_t> group_index_;
  std::vector<std::string> joint_names_;
  std::vector<std::string> link_names_;
};

}  // namespace varietas_moveit

#endif
