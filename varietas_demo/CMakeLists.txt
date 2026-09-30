cmake_minimum_required(VERSION 3.16)
project(varietas_demo)

find_package(ament_cmake REQUIRED)
find_package(rclcpp REQUIRED)
find_package(sensor_msgs REQUIRED)
find_package(visualization_msgs REQUIRED)
find_package(urdf REQUIRED)
find_package(tf2_ros REQUIRED)
find_package(geometry_msgs REQUIRED)
find_package(std_msgs REQUIRED)
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

# The solve-all demonstration draws an arm the decoupling refuses, so its
# solver is the full one over Q(x, y, z), reconstructed from fixed poses. That
# takes about half a minute and is paid here, once, rather than by anyone who
# runs the demonstration. find_package(Eigen3) is for the iteration the node
# sets beside the generated solver; the generated header needs it anyway.
find_package(Eigen3 REQUIRED)

add_executable(generate_general_solver src/generate_general_solver.cpp)
target_compile_features(generate_general_solver PRIVATE cxx_std_17)
ament_target_dependencies(generate_general_solver urdf)
target_link_libraries(generate_general_solver varietas_urdf::varietas_urdf)

set(general_arm ${CMAKE_CURRENT_SOURCE_DIR}/urdf/general_3r.urdf)
set(general_header ${generated_dir}/general_ik.hpp)
add_custom_command(
  OUTPUT ${general_header}
  COMMAND ${CMAKE_COMMAND} -E make_directory ${generated_dir}
  COMMAND generate_general_solver ${general_arm} ${general_header}
  DEPENDS generate_general_solver ${general_arm}
  COMMENT "Reconstructing the full solver of general_3r from fixed poses into ${general_header}")
add_custom_target(general_header_target DEPENDS ${general_header})

add_executable(solve_all_node src/solve_all_node.cpp)
target_compile_features(solve_all_node PRIVATE cxx_std_17)
ament_target_dependencies(solve_all_node rclcpp visualization_msgs std_msgs urdf)
target_link_libraries(solve_all_node varietas_urdf::varietas_urdf Eigen3::Eigen)
target_include_directories(solve_all_node PRIVATE ${generated_dir})
add_dependencies(solve_all_node general_header_target)

install(TARGETS sweep_node branch_node solve_all_node DESTINATION lib/${PROJECT_NAME})
install(DIRECTORY launch rviz urdf DESTINATION share/${PROJECT_NAME})

ament_package()
