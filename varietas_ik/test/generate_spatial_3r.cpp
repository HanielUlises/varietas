// A solver for the anthropomorphic arm, written at build time.
//
// The arm this comes from does not solve over Q(x, y, z) in any time worth
// waiting for. What is emitted here is the two-joint problem left after its
// base joint is swept out, which solves in a fraction of a second, together
// with the arctangent that puts the base joint back — so the header is a
// solver for the whole arm and not for most of it.

#include <cstdio>
#include <fstream>
#include <string>

#include "varietas/codegen/emit.hpp"
#include "varietas/ik/decoupled_ik.hpp"
#include "varietas/ik/emit_decoupled.hpp"

#include "arms.hpp"

int main(int argc, char** argv) {
  if (argc != 2) {
    std::fprintf(stderr, "usage: %s <output-header>\n", argv[0]);
    return 2;
  }

  const auto result =
      varietas::ik::decoupled_position_ik<3>(varietas_test::anthropomorphic_three_link());
  if (!result.ok()) {
    std::fprintf(stderr, "the arm did not decouple: %s\n",
                 varietas::ik::to_string(result.status));
    return 1;
  }

  varietas::codegen::emit_options options;
  options.name = "anthropomorphic_ik";
  options.name_space = "varietas_generated";
  options.source_note =
      "The anthropomorphic 3R arm. Its base yaw was swept out rather than adjoined: the "
      "struct below solves shoulder and elbow against a radius and a height over "
      "Q(radius, height), and the wrapper after it puts the base joint back.";
  options.runtime = varietas::codegen::runtime_kind::eigen;

  const std::string header = varietas::ik::emit_decoupled(result, options);

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
