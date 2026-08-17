cmake_minimum_required(VERSION 3.16)
project(varietas_codegen)

find_package(ament_cmake REQUIRED)
find_package(varietas_core REQUIRED)

# GMP ships no config file and no reliable .pc, so it is located by hand. Both
# the C library and the C++ wrapper are needed: mpq_class lives in gmpxx.
find_path(GMPXX_INCLUDE_DIR gmpxx.h)
find_library(GMP_LIBRARY NAMES gmp)
find_library(GMPXX_LIBRARY NAMES gmpxx)

if(NOT GMPXX_INCLUDE_DIR OR NOT GMP_LIBRARY OR NOT GMPXX_LIBRARY)
  message(FATAL_ERROR
    "GMP with its C++ interface is required by varietas_codegen. "
    "Install libgmp-dev (Debian/Ubuntu) or gmp-devel (Fedora).")
endif()

add_library(${PROJECT_NAME} INTERFACE)
target_compile_features(${PROJECT_NAME} INTERFACE cxx_std_17)
target_include_directories(${PROJECT_NAME} INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
  ${GMPXX_INCLUDE_DIR})
target_link_libraries(${PROJECT_NAME} INTERFACE
  varietas_core::varietas_core
  ${GMPXX_LIBRARY}
  ${GMP_LIBRARY})

install(DIRECTORY include/ DESTINATION include)
install(TARGETS ${PROJECT_NAME} EXPORT ${PROJECT_NAME}_export)

ament_export_targets(${PROJECT_NAME}_export HAS_LIBRARY_TARGET)
ament_export_dependencies(varietas_core)
ament_export_include_directories(${GMPXX_INCLUDE_DIR})
ament_export_libraries(${GMPXX_LIBRARY} ${GMP_LIBRARY})

if(BUILD_TESTING)
  find_package(ament_cmake_gtest REQUIRED)

  ament_add_gtest(test_rational test/test_rational.cpp)
  target_link_libraries(test_rational ${PROJECT_NAME})

  # The core algebra and ideal suites, run a second time over the exact field.
  ament_add_gtest(test_exact_algebra test/test_exact_algebra.cpp)
  target_link_libraries(test_exact_algebra ${PROJECT_NAME})

  ament_add_gtest(test_exact_groebner test/test_exact_groebner.cpp)
  target_link_libraries(test_exact_groebner ${PROJECT_NAME})
endif()

ament_package()
