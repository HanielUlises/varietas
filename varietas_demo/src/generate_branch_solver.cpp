// The solver the branch demonstration draws, written at build time.
//
// This is the same path `urdf_codegen --decouple` takes, run here so that the
// demonstration compiles against a header varietas produced during this build
// rather than against a copy checked in beside it. If the emitter regresses,
// the demonstration stops building, which is the only guarantee worth having.
//
// Three joints, because that is what the decomposition is for: the base yaw is
// swept out and the shoulder and elbow are solved over Q(radius, height). An
// arm of another size is refused here rather than silently reinterpreted.

#include <cstdio>
#include <fstream>
#include <string>

#include <urdf/model.h>

#include "varietas/codegen/emit.hpp"
#include "varietas/ik/decoupled_ik.hpp"
#include "varietas/ik/emit_decoupled.hpp"
#include "varietas/urdf/urdf_chain.hpp"

int main(int argc, char** argv) {
  if (argc != 3) {
    std::fprintf(stderr, "usage: %s <file.urdf> <output.hpp>\n", argv[0]);
    return 2;
  }

  urdf::Model model;
  if (!model.initFile(argv[1])) {
    std::fprintf(stderr, "could not parse %s\n", argv[1]);
    return 1;
  }

  const std::string tip = varietas::urdf_import::sole_tip_link(model);
  if (tip.empty()) {
    std::fprintf(stderr, "%s is branched; the demonstration wants a serial arm\n", argv[1]);
    return 1;
  }

  varietas::chain<varietas::rational> exact;
  const auto report =
      varietas::urdf_import::chain_from_model(model, model.getRoot()->name, tip, exact);
  if (!report.ok()) {
    std::fprintf(stderr, "import refused: %s (%s)\n",
                 varietas::urdf_import::to_string(report.status), report.detail.c_str());
    return 1;
  }

  if (exact.degrees_of_freedom() != 3) {
    std::fprintf(stderr,
                 "%s has %zu actuated joints; the decoupled demonstration is written for "
                 "three, one to sweep out and two to solve for\n",
                 argv[1], exact.degrees_of_freedom());
    return 1;
  }

  const auto solution = varietas::ik::decoupled_position_ik<3>(exact);
  if (!solution.ok()) {
    std::fprintf(stderr, "the arm did not decouple: %s\n",
                 varietas::ik::to_string(solution.status));
    return 1;
  }

  varietas::codegen::emit_options options;
  options.name = "branch_ik";
  options.name_space = "varietas_generated";
  options.runtime = varietas::codegen::runtime_kind::eigen;
  options.source_note =
      "Emitted for the varietas_demo branch demonstration. The base joint was swept out "
      "rather than adjoined; the wrapper returns joint angles for the whole arm.";

  std::ofstream out(argv[2]);
  if (!out) {
    std::fprintf(stderr, "cannot write %s\n", argv[2]);
    return 1;
  }
  out << varietas::ik::emit_decoupled(solution, options);
  if (!out) {
    std::fprintf(stderr, "write to %s failed\n", argv[2]);
    return 1;
  }

  std::fprintf(stderr, "%s: %zu branches, wrote %s\n", exact.name().c_str(), solution.branches,
               argv[2]);
  return 0;
}
