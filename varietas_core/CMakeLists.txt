cmake_minimum_required(VERSION 3.16)
project(varietas_core)

find_package(ament_cmake REQUIRED)
find_package(Eigen3 REQUIRED)

add_library(${PROJECT_NAME} INTERFACE)
target_compile_features(${PROJECT_NAME} INTERFACE cxx_std_17)
target_include_directories(${PROJECT_NAME} INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>)
target_link_libraries(${PROJECT_NAME} INTERFACE Eigen3::Eigen)

install(DIRECTORY include/ DESTINATION include)
install(TARGETS ${PROJECT_NAME} EXPORT ${PROJECT_NAME}_export)

ament_export_targets(${PROJECT_NAME}_export HAS_LIBRARY_TARGET)
ament_export_dependencies(Eigen3)

if(BUILD_TESTING)
  find_package(ament_cmake_gtest REQUIRED)
  ament_add_gtest(test_algebra test/test_algebra.cpp)
  target_link_libraries(test_algebra ${PROJECT_NAME})

  ament_add_gtest(test_ideal test/test_ideal.cpp)
  target_link_libraries(test_ideal ${PROJECT_NAME})

  ament_add_gtest(test_solve test/test_solve.cpp)
  target_link_libraries(test_solve ${PROJECT_NAME})

  ament_add_gtest(test_saturation test/test_saturation.cpp)
  target_link_libraries(test_saturation ${PROJECT_NAME})
endif()

ament_package()
