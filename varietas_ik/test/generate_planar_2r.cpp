// The pipeline, end to end, run at build time.
//
// A chain goes in and a header comes out, and nothing in between is told what
// pose it will be asked about. This is the program the README's missing link
// was: varietas_kinematics poses the problem, varietas_ik solves it over
// Q(x, y), varietas_codegen writes it down.

#include <cstdio>
#include <fstream>
#include <string>

#include "varietas/codegen/emit.hpp"
#include "varietas/ik/parametric_ik.hpp"

#include "arms.hpp"

int main(int argc, char** argv) {
  if (argc != 2) {
    std::fprintf(stderr, "usage: %s <output-header>\n", argv[0]);
    return 2;
  }

  // The arm is planar, so only x and y are adjoined; the bridge checks that the
  // z equation is vacuous before agreeing to drop it.
  const auto result =
      varietas::ik::parametric_position_ik<2, 2>(varietas_test::planar_two_link(), {0, 1});
  if (!result.ok()) {
    std::fprintf(stderr, "parametric inverse kinematics failed: %s\n",
                 varietas::ik::to_string(result.status));
    return 1;
  }

  varietas::codegen::emit_options options;
  options.name = "planar_2r_ik";
  options.name_space = "varietas_generated";
  options.source_note =
      "Position inverse kinematics of the planar 2R arm of unit links, solved once over "
      "Q(x, y) by varietas_ik. The unknowns are t = tan(q/2).";
  options.runtime = varietas::codegen::runtime_kind::eigen;

  const std::string header = varietas::codegen::emit(result.solution, options);

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
