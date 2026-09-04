// Writes the header that test_generated.cpp compiles.
//
// This is the offline half in miniature and in the shape it will really be
// used: a program that solves a parametric system once and leaves a header
// behind. Running it as part of the build is what makes test_generated a test
// of generated code rather than of a checked-in artefact that may no longer
// correspond to the emitter.

#include <cstdio>
#include <fstream>
#include <string>

#include "varietas/codegen/emit.hpp"

#include "planar_fixture.hpp"

int main(int argc, char** argv) {
  if (argc != 2) {
    std::fprintf(stderr, "usage: %s <output-header>\n", argv[0]);
    return 2;
  }

  varietas::codegen::emit_options options;
  options.name = "planar";
  options.name_space = "varietas_generated";
  options.source_note = "u^2 - x, v - y u, solved over Q(x, y) by varietas_codegen's test fixture.";
  options.runtime = varietas::codegen::runtime_kind::eigen;

  const std::string header = varietas::codegen::emit(varietas_test::planar_solution(), options);

  std::ofstream out(argv[1]);
  if (!out) {
    std::fprintf(stderr, "cannot write %s\n", argv[1]);
    return 1;
  }
  out << header;
  if (!out) {
    std::fprintf(stderr, "write to %s failed\n", argv[1]);
    return 1;
  }
  return 0;
}
