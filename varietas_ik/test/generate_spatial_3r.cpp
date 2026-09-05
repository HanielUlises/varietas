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

namespace {

// One arm, solved and written out. Returns false with a reason on stderr.
bool emit_arm(const varietas::chain<varietas::rational>& robot, const std::string& name,
              const std::string& note, const char* path) {
  const auto result = varietas::ik::decoupled_position_ik<3>(robot);
  if (!result.ok()) {
    std::fprintf(stderr, "the arm %s did not decouple: %s\n", name.c_str(),
                 varietas::ik::to_string(result.status));
    return false;
  }

  varietas::codegen::emit_options options;
  options.name = name;
  options.name_space = "varietas_generated";
  options.source_note = note;
  options.runtime = varietas::codegen::runtime_kind::eigen;

  std::ofstream out(path);
  if (!out) {
    std::fprintf(stderr, "cannot write %s\n", path);
    return false;
  }
  out << varietas::ik::emit_decoupled(result, options);
  if (!out) {
    std::fprintf(stderr, "write to %s failed\n", path);
    return false;
  }
  return true;
}

}  // namespace

int main(int argc, char** argv) {
  if (argc != 3) {
    std::fprintf(stderr, "usage: %s <header> <reversed-header>\n", argv[0]);
    return 2;
  }

  if (!emit_arm(varietas_test::anthropomorphic_three_link(), "anthropomorphic_ik",
                "The anthropomorphic 3R arm. Its base yaw was swept out rather than adjoined: "
                "the struct below solves shoulder and elbow against a radius and a height over "
                "Q(radius, height), and the wrapper after it puts the base joint back.",
                argv[1])) {
    return 1;
  }

  // The same arm with its base axis reversed. The reduced problem is identical;
  // the sign the wrapper applies to the arctangent is not, and emitting both is
  // the only way that sign gets compiled and run.
  if (!emit_arm(varietas_test::anthropomorphic_reversed_base(), "anthropomorphic_reversed_ik",
                "The anthropomorphic 3R arm with its base turning about -z. Identical to the "
                "forward version except for the sign of the recovered base angle.",
                argv[2])) {
    return 1;
  }

  return 0;
}
