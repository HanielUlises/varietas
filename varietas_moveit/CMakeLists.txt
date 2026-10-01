cmake_minimum_required(VERSION 3.16)
project(varietas_moveit)

find_package(ament_cmake REQUIRED)
find_package(moveit_core REQUIRED)
find_package(pluginlib REQUIRED)
find_package(rclcpp REQUIRED)
find_package(urdf REQUIRED)
find_package(varietas_ik REQUIRED)
find_package(varietas_urdf REQUIRED)
find_package(Eigen3 REQUIRED)

add_library(varietas_moveit_plugin SHARED src/varietas_kinematics_plugin.cpp)
target_compile_features(varietas_moveit_plugin PUBLIC cxx_std_17)
target_include_directories(varietas_moveit_plugin PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>)
ament_target_dependencies(varietas_moveit_plugin moveit_core pluginlib rclcpp urdf)
target_link_libraries(varietas_moveit_plugin varietas_urdf::varietas_urdf Eigen3::Eigen)

pluginlib_export_plugin_description_file(moveit_core varietas_kinematics_plugin_description.xml)

install(TARGETS varietas_moveit_plugin
  EXPORT export_varietas_moveit
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  RUNTIME DESTINATION bin)
install(DIRECTORY include/ DESTINATION include)

if(BUILD_TESTING)
  find_package(ament_cmake_gtest REQUIRED)
  find_package(ament_index_cpp REQUIRED)

  # The plugin is loaded by name through pluginlib, as MoveIt loads it, so the
  # test also checks that it is exported where MoveIt looks.
  ament_add_gtest(test_varietas_plugin test/test_plugin.cpp TIMEOUT 300)
  target_link_libraries(test_varietas_plugin varietas_moveit_plugin)
  ament_target_dependencies(test_varietas_plugin moveit_core pluginlib rclcpp ament_index_cpp)

  # Not a test: a comparison with MoveIt's default solver, run by hand.
  add_executable(compare_with_kdl test/compare_with_kdl.cpp)
  target_link_libraries(compare_with_kdl varietas_moveit_plugin)
  ament_target_dependencies(compare_with_kdl moveit_core pluginlib rclcpp ament_index_cpp)
endif()

ament_export_dependencies(moveit_core pluginlib rclcpp)
ament_package()
