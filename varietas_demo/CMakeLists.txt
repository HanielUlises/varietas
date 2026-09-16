cmake_minimum_required(VERSION 3.16)
project(varietas_demo)

find_package(ament_cmake REQUIRED)
find_package(rclcpp REQUIRED)
find_package(sensor_msgs REQUIRED)
find_package(visualization_msgs REQUIRED)
find_package(urdf REQUIRED)
find_package(tf2_ros REQUIRED)
find_package(geometry_msgs REQUIRED)
find_package(varietas_urdf REQUIRED)

add_executable(sweep_node src/sweep_node.cpp)
target_compile_features(sweep_node PRIVATE cxx_std_17)
ament_target_dependencies(sweep_node rclcpp sensor_msgs visualization_msgs geometry_msgs tf2_ros urdf)
target_link_libraries(sweep_node varietas_urdf::varietas_urdf)

# The branch demonstration draws every configuration the generated solver
# returns, so it needs a generated solver. It is emitted here, during this
# build, from the URDF, by the path urdf_codegen --decouple takes: the
# demonstration compiles against varietas's own output rather than a copy
# checked in beside it, and an emitter that regresses stops the build.
add_executable(generate_branch_solver src/generate_branch_solver.cpp)
target_compile_features(generate_branch_solver PRIVATE cxx_std_17)
ament_target_dependencies(generate_branch_solver urdf)
target_link_libraries(generate_branch_solver varietas_urdf::varietas_urdf)

set(branch_arm ${CMAKE_CURRENT_SOURCE_DIR}/urdf/anthropomorphic_offset_3r.urdf)
set(generated_dir ${CMAKE_CURRENT_BINARY_DIR}/generated)
set(branch_header ${generated_dir}/branch_ik.hpp)
add_custom_command(
  OUTPUT ${branch_header}
  COMMAND ${CMAKE_COMMAND} -E make_directory ${generated_dir}
  COMMAND generate_branch_solver ${branch_arm} ${branch_header}
  DEPENDS generate_branch_solver ${branch_arm}
  COMMENT "Sweeping the base joint of the anthropomorphic arm and emitting ${branch_header}")
add_custom_target(branch_header_target DEPENDS ${branch_header})

add_executable(branch_node src/branch_node.cpp)
target_compile_features(branch_node PRIVATE cxx_std_17)
ament_target_dependencies(branch_node rclcpp sensor_msgs visualization_msgs geometry_msgs tf2_ros urdf)
target_link_libraries(branch_node varietas_urdf::varietas_urdf)
target_include_directories(branch_node PRIVATE ${generated_dir})
add_dependencies(branch_node branch_header_target)

install(TARGETS sweep_node branch_node DESTINATION lib/${PROJECT_NAME})
install(DIRECTORY launch rviz urdf DESTINATION share/${PROJECT_NAME})

ament_package()
