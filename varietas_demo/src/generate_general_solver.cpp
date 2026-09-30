// The solver the solve-all demonstration draws, written at build time.
//
// The arm is general_3r, whose three axes are pairwise skew. The decoupling
// refuses it and the symbolic solve over Q(x, y, z) does not finish, so this
// runs the reconstruction from fixed poses, the path `urdf_codegen
// --reconstruct` takes, and emits the full three-parameter solver. It takes
// about half a minute, which the build pays once so that no caller does.
//
// What the reconstruction found out about itself is written into the header
// too, as a struct of constants after the solver: how long it took, how many
// primes and samples it used, how many exact checks it passed, and the reason
// the decoupling gave for refusing. The demonstration shows those facts rather
// than restating them, so they cannot drift from what actually happened.

#include <chrono>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>

#include <urdf/model.h>

#include "varietas/codegen/emit.hpp"
#include "varietas/ik/decoupled_ik.hpp"
#include "varietas/ik/reconstructed_ik.hpp"
#include "varietas/urdf/urdf_chain.hpp"

namespace {

// A C++ string literal holding `text`.
std::string quoted(const std::string& text) {
  std::string out = "\"";
  for (const char c : text) {
    if (c == '"' || c == '\\') {
      out += '\\';
    }
    out += c;
  }
  return out + "\"";
}

}  // namespace

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
  varietas::chain<varietas::rational> exact;
  const auto imported =
      varietas::urdf_import::chain_from_model(model, model.getRoot()->name, tip, exact);
  if (!imported.ok()) {
    std::fprintf(stderr, "import refused: %s (%s)\n",
                 varietas::urdf_import::to_string(imported.status), imported.detail.c_str());
    return 1;
  }
  if (exact.degrees_of_freedom() != 3) {
    std::fprintf(stderr, "%s has %zu actuated joints; the demonstration is written for three\n",
                 argv[1], exact.degrees_of_freedom());
    return 1;
  }

  // Asked first, so that the refusal is the library's own sentence.
  const auto decoupled = varietas::ik::decoupled_position_ik<3>(exact);
  const std::string decoupling =
      decoupled.ok() ? std::string("admitted") : varietas::ik::to_string(decoupled.status);

  varietas::ik::reconstruction_report report;
  const auto start = std::chrono::steady_clock::now();
  const auto solved =
      varietas::ik::reconstructed_position_ik<3, 3>(exact, {0, 1, 2}, &report);
  const double seconds =
      std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
  if (!solved.ok()) {
    std::fprintf(stderr, "reconstruction refused: %s\n",
                 varietas::ik::to_string(solved.status));
    return 1;
  }

  std::ostringstream provenance;
  provenance << "// What the reconstruction that produced general_ik found out about itself.\n"
             << "struct general_ik_provenance {\n"
             << "  static constexpr double seconds = " << seconds << ";\n"
             << "  static constexpr int primes = " << report.sampling.primes << ";\n"
             << "  static constexpr int samples = " << report.sampling.samples << ";\n"
             << "  static constexpr int exact_checks = " << report.exact_checks << ";\n"
             << "  static constexpr int branches = " << solved.branches << ";\n"
             << "  static constexpr const char* decoupling = " << quoted(decoupling) << ";\n"
             << "};";

  varietas::codegen::emit_options options;
  options.name = "general_ik";
  options.name_space = "varietas_generated";
  options.runtime = varietas::codegen::runtime_kind::eigen;
  options.source_note =
      "Emitted for the varietas_demo solve-all demonstration, by reconstruction from "
      "fixed-pose solves over prime fields, checked exactly at rational poses. The unknowns "
      "are t = tan(q/2).";
  options.epilogue = provenance.str();

  std::ofstream out(argv[2]);
  if (!out) {
    std::fprintf(stderr, "cannot write %s\n", argv[2]);
    return 1;
  }
  out << varietas::codegen::emit(solved.solution, options);
  if (!out) {
    std::fprintf(stderr, "write to %s failed\n", argv[2]);
    return 1;
  }
  std::fprintf(stderr, "%s: %zu branches in %.1f s over %zu primes, decoupling %s; wrote %s\n",
               exact.name().c_str(), solved.branches, seconds, report.sampling.primes,
               decoupling.c_str(), argv[2]);
  return 0;
}
