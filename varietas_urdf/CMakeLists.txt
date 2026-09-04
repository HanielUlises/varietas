cmake_minimum_required(VERSION 3.16)
project(varietas_urdf)

find_package(ament_cmake REQUIRED)
find_package(urdf REQUIRED)
find_package(varietas_codegen REQUIRED)
find_package(varietas_ik REQUIRED)
find_package(varietas_kinematics REQUIRED)

add_library(${PROJECT_NAME} INTERFACE)
target_compile_features(${PROJECT_NAME} INTERFACE cxx_std_17)
target_include_directories(${PROJECT_NAME} INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>)
target_link_libraries(${PROJECT_NAME} INTERFACE
  varietas_kinematics::varietas_kinematics
  varietas_codegen::varietas_codegen
  varietas_ik::varietas_ik
  ${urdf_TARGETS})

add_executable(urdf_report src/urdf_report.cpp)
target_link_libraries(urdf_report ${PROJECT_NAME})
install(TARGETS urdf_report DESTINATION lib/${PROJECT_NAME})

# The pipeline as a command: a URDF in, a header out.
add_executable(urdf_codegen src/urdf_codegen.cpp)
target_link_libraries(urdf_codegen ${PROJECT_NAME})
install(TARGETS urdf_codegen DESTINATION lib/${PROJECT_NAME})

install(DIRECTORY include/ DESTINATION include)
install(TARGETS ${PROJECT_NAME} EXPORT ${PROJECT_NAME}_export)
install(DIRECTORY test/data DESTINATION share/${PROJECT_NAME})

ament_export_targets(${PROJECT_NAME}_export HAS_LIBRARY_TARGET)
ament_export_dependencies(varietas_kinematics varietas_codegen varietas_ik urdf)

if(BUILD_TESTING)
  find_package(ament_cmake_gtest REQUIRED)
  find_package(kdl_parser REQUIRED)
  find_package(orocos_kdl REQUIRED)

  # KDL is a test dependency only: it is the reference the recovered chain is
  # checked against, not something the pipeline runs on.
  ament_add_gtest(test_urdf_chain test/test_urdf_chain.cpp)
  target_compile_definitions(test_urdf_chain PRIVATE
    VARIETAS_URDF_TEST_DATA="${CMAKE_CURRENT_SOURCE_DIR}/test/data")
  target_link_libraries(test_urdf_chain ${PROJECT_NAME} ${kdl_parser_TARGETS} orocos-kdl)
endif()

ament_package()
