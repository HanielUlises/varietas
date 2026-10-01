// Six-joint solvers with a spherical wrist, written at build time.
//
// Two arms. The industrial one decouples at the arm, so its header carries the
// decoupled three-joint solver; the skewed one does not, so its arm is
// reconstructed from fixed poses over Q(x, y, z). Between them every path the
// wrist emitter can take is compiled and run by test_spherical_wrist.

#include <cstdio>
#include <fstream>
#include <string>

#include "varietas/ik/decoupled_ik.hpp"
#include "varietas/ik/reconstructed_ik.hpp"
#include "varietas/ik/spherical_wrist.hpp"

#include "arms.hpp"

namespace {

bool write(const varietas::chain<varietas::rational>& robot, const std::string& name,
           const char* path) {
  const auto wrist = varietas::ik::decompose_spherical_wrist(robot);
  if (!wrist.ok()) {
    std::fprintf(stderr, "%s: %s\n", name.c_str(), varietas::ik::to_string(wrist.status));
    return false;
  }
  varietas::codegen::emit_options options;
  options.name = name;
  options.runtime = varietas::codegen::runtime_kind::eigen;
  std::string text;
  const auto decoupled = varietas::ik::decoupled_position_ik<3>(wrist.arm);
  if (decoupled.ok()) {
    text = varietas::ik::emit_spherical_wrist(wrist, decoupled, options);
  } else {
    const auto reconstructed =
        varietas::ik::reconstructed_position_ik<3, 3>(wrist.arm, {0, 1, 2});
    if (!reconstructed.ok()) {
      std::fprintf(stderr, "%s: the arm was not solved\n", name.c_str());
      return false;
    }
    text = varietas::ik::emit_spherical_wrist(wrist, reconstructed.solution, options);
  }
  std::ofstream out(path);
  out << text;
  return static_cast<bool>(out);
}

}  // namespace

int main(int argc, char** argv) {
  if (argc != 3) {
    std::fprintf(stderr, "usage: %s <industrial header> <skewed header>\n", argv[0]);
    return 2;
  }
  return write(varietas_test::industrial_six(), "industrial_ik", argv[1]) &&
                 write(varietas_test::skewed_wrist_six(), "skewed_ik", argv[2])
             ? 0
             : 1;
}
