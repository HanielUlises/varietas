cmake_minimum_required(VERSION 3.16)
project(varietas_kinematics)

find_package(ament_cmake REQUIRED)
find_package(varietas_core REQUIRED)

add_library(${PROJECT_NAME} INTERFACE)
target_compile_features(${PROJECT_NAME} INTERFACE cxx_std_17)
target_include_directories(${PROJECT_NAME} INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>)
target_link_libraries(${PROJECT_NAME} INTERFACE varietas_core::varietas_core)

install(DIRECTORY include/ DESTINATION include)
install(TARGETS ${PROJECT_NAME} EXPORT ${PROJECT_NAME}_export)

ament_export_targets(${PROJECT_NAME}_export HAS_LIBRARY_TARGET)
ament_export_dependencies(varietas_core)

if(BUILD_TESTING)
  find_package(ament_cmake_gtest REQUIRED)

  # The library itself needs nothing but varietas_core; the exact field is a
  # test dependency only, so that the same chain code is exercised over double
  # and over the rationals the offline pipeline will actually use.
  find_package(varietas_codegen REQUIRED)

  ament_add_gtest(test_chain test/test_chain.cpp)
  target_link_libraries(test_chain ${PROJECT_NAME} varietas_codegen::varietas_codegen)

  ament_add_gtest(test_rationalize test/test_rationalize.cpp)
  target_link_libraries(test_rationalize ${PROJECT_NAME} varietas_codegen::varietas_codegen)

  ament_add_gtest(test_trigonometric test/test_trigonometric.cpp TIMEOUT 300)
  target_link_libraries(test_trigonometric ${PROJECT_NAME} varietas_codegen::varietas_codegen)

  ament_add_gtest(test_workspace test/test_workspace.cpp TIMEOUT 300)
  target_link_libraries(test_workspace ${PROJECT_NAME} varietas_codegen::varietas_codegen)
endif()

ament_package()
