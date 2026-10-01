#ifndef VARIETAS_MOVEIT_TEST_ROBOT_HPP
#define VARIETAS_MOVEIT_TEST_ROBOT_HPP

#include <memory>
#include <string>

#include <ament_index_cpp/get_package_share_directory.hpp>
#include <moveit/robot_model/robot_model.h>
#include <srdfdom/model.h>
#include <urdf/model.h>

// The industrial arm from varietas_urdf's test data, with the one planning
// group a MoveIt configuration for it would declare: the chain from the base
// to the flange.
inline moveit::core::RobotModelPtr industrial_robot_model() {
  const std::string path = ament_index_cpp::get_package_share_directory("varietas_urdf") +
                           "/data/industrial_6r.urdf";
  auto urdf = std::make_shared<urdf::Model>();
  if (!urdf->initFile(path)) {
    return nullptr;
  }
  const std::string srdf_text =
      "<?xml version=\"1.0\"?><robot name=\"industrial_6r\">"
      "<group name=\"manipulator\"><chain base_link=\"base_link\" tip_link=\"flange\"/></group>"
      "</robot>";
  auto srdf = std::make_shared<srdf::Model>();
  if (!srdf->initString(*urdf, srdf_text)) {
    return nullptr;
  }
  return std::make_shared<moveit::core::RobotModel>(urdf, srdf);
}

#endif
