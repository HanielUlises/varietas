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

  ament_add_gtest(test_rational_function test/test_rational_function.cpp)
  target_link_libraries(test_rational_function ${PROJECT_NAME})

  ament_add_gtest(test_emit test/test_emit.cpp)
  target_link_libraries(test_emit ${PROJECT_NAME})
  target_include_directories(test_emit PRIVATE test)

  # Generated code, compiled.
  #
  # The emitter's other tests compare text or walk the same expression tree
  # twice; neither can catch output that does not parse, or that parses and
  # means something else. So the header is generated during the build by a
  # program that links the emitter, and test_generated includes it as ordinary
  # source. A syntax error in the emitted text is a build failure here, which
  # is where it should be, and the numbers are checked against the closed form
  # of the same system.
  add_executable(generate_planar test/generate_planar.cpp)
  target_link_libraries(generate_planar ${PROJECT_NAME})
  target_include_directories(generate_planar PRIVATE test)

  set(generated_dir ${CMAKE_CURRENT_BINARY_DIR}/generated)
  set(generated_header ${generated_dir}/planar.hpp)
  add_custom_command(
    OUTPUT ${generated_header}
    COMMAND ${CMAKE_COMMAND} -E make_directory ${generated_dir}
    COMMAND generate_planar ${generated_header}
    DEPENDS generate_planar
    COMMENT "Emitting ${generated_header} from the parametric solution")
  add_custom_target(planar_header DEPENDS ${generated_header})

  ament_add_gtest(test_generated test/test_generated.cpp)
  target_link_libraries(test_generated ${PROJECT_NAME})
  target_include_directories(test_generated PRIVATE ${generated_dir})
  add_dependencies(test_generated planar_header)

  # The core algebra and ideal suites, run a second time over the exact field.
  ament_add_gtest(test_exact_algebra test/test_exact_algebra.cpp)
  target_link_libraries(test_exact_algebra ${PROJECT_NAME})

  ament_add_gtest(test_exact_groebner test/test_exact_groebner.cpp)
  target_link_libraries(test_exact_groebner ${PROJECT_NAME})

  ament_add_gtest(test_exact_saturation test/test_exact_saturation.cpp)
  target_link_libraries(test_exact_saturation ${PROJECT_NAME})
endif()

ament_package()
