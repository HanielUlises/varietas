cmake_minimum_required(VERSION 3.16)
project(varietas_ik)

find_package(ament_cmake REQUIRED)
find_package(varietas_core REQUIRED)
find_package(varietas_codegen REQUIRED)
find_package(varietas_kinematics REQUIRED)

add_library(${PROJECT_NAME} INTERFACE)
target_compile_features(${PROJECT_NAME} INTERFACE cxx_std_17)
target_include_directories(${PROJECT_NAME} INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>)
target_link_libraries(${PROJECT_NAME} INTERFACE
  varietas_core::varietas_core
  varietas_codegen::varietas_codegen
  varietas_kinematics::varietas_kinematics)

install(DIRECTORY include/ DESTINATION include)
install(TARGETS ${PROJECT_NAME} EXPORT ${PROJECT_NAME}_export)

ament_export_targets(${PROJECT_NAME}_export HAS_LIBRARY_TARGET)
ament_export_dependencies(varietas_core varietas_codegen varietas_kinematics)

if(BUILD_TESTING)
  find_package(ament_cmake_gtest REQUIRED)

  ament_add_gtest(test_parametric_ik test/test_parametric_ik.cpp TIMEOUT 300)
  target_link_libraries(test_parametric_ik ${PROJECT_NAME})
  target_include_directories(test_parametric_ik PRIVATE test)

  # The same discipline varietas_codegen applies to its own emitter, applied to
  # the whole pipeline: a program that links the bridge solves a real arm over
  # Q(x, y) during the build and leaves a header behind, and the test below
  # includes that header as ordinary source. Text that does not parse is a build
  # failure, and the numbers are checked against the closed form of the arm.
  add_executable(generate_planar_2r test/generate_planar_2r.cpp)
  target_link_libraries(generate_planar_2r ${PROJECT_NAME})
  target_include_directories(generate_planar_2r PRIVATE test)

  set(generated_dir ${CMAKE_CURRENT_BINARY_DIR}/generated)
  set(generated_header ${generated_dir}/planar_2r_ik.hpp)
  add_custom_command(
    OUTPUT ${generated_header}
    COMMAND ${CMAKE_COMMAND} -E make_directory ${generated_dir}
    COMMAND generate_planar_2r ${generated_header}
    DEPENDS generate_planar_2r
    COMMENT "Solving the planar 2R arm over Q(x, y) and emitting ${generated_header}")
  add_custom_target(planar_2r_header DEPENDS ${generated_header})

  ament_add_gtest(test_generated_ik test/test_generated_ik.cpp TIMEOUT 300)
  target_link_libraries(test_generated_ik ${PROJECT_NAME})
  target_include_directories(test_generated_ik PRIVATE ${generated_dir} test)
  add_dependencies(test_generated_ik planar_2r_header)
endif()

ament_package()
