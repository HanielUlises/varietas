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

install(TARGETS sweep_node DESTINATION lib/${PROJECT_NAME})
install(DIRECTORY launch rviz DESTINATION share/${PROJECT_NAME})

ament_package()
