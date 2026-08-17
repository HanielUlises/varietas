// Reads a URDF and reports what it costs to make it exact, joint by joint.
//
// Usage: urdf_report <file.urdf> [tip_link] [root_link]
//
// The table is the audit a model has to pass before any Gröbner basis computed
// from it means anything: which placements were recovered without moving, which
// were moved onto the exact geometry the file's decimals were approximating,
// and by how far.

#include <cstdio>
#include <string>

#include <urdf/model.h>

#include "varietas/kinematics/evaluate.hpp"
#include "varietas/urdf/urdf_chain.hpp"

int main(int argc, char** argv) {
  if (argc < 2) {
    std::fprintf(stderr, "usage: %s <file.urdf> [tip_link] [root_link]\n", argv[0]);
    return 2;
  }

  urdf::Model model;
  if (!model.initFile(argv[1])) {
    std::fprintf(stderr, "could not parse %s\n", argv[1]);
    return 1;
  }

  const std::string tip =
      argc > 2 ? argv[2] : varietas::urdf_import::sole_tip_link(model);
  const std::string root = argc > 3 ? argv[3] : model.getRoot()->name;
  if (tip.empty()) {
    std::fprintf(stderr, "%s is branched; name the tip link explicitly\n", argv[1]);
    return 1;
  }

  varietas::chain<varietas::rational> robot;
  const auto report = varietas::urdf_import::chain_from_model(model, root, tip, robot);

  std::printf("model            %s\n", model.getName().c_str());
  std::printf("chain            %s -> %s\n", root.c_str(), tip.c_str());

  if (!report.ok()) {
    std::printf("refused          %s (%s)\n",
                varietas::urdf_import::to_string(report.status), report.detail.c_str());
    return 1;
  }

  std::printf("degrees of freedom %zu\n\n", robot.degrees_of_freedom());
  std::printf("%-24s %-10s %14s %14s\n", "joint", "type", "rotation", "translation");
  std::printf("%-24s %-10s %14s %14s\n", "", "", "moved (rad)", "moved (m)");
  std::printf("%s\n", std::string(66, '-').c_str());

  for (std::size_t i = 0; i < report.joints.size(); ++i) {
    const auto& recovery = report.joints[i];
    const char* type = "fixed";
    switch (robot.joints()[i].type) {
      case varietas::joint_type::revolute:
        type = "revolute";
        break;
      case varietas::joint_type::prismatic:
        type = "prismatic";
        break;
      case varietas::joint_type::fixed:
        break;
    }
    std::printf("%-24s %-10s %14.3e %14.3e%s\n", recovery.name.c_str(), type,
                recovery.rotation_deviation, recovery.translation_deviation,
                recovery.unmoved ? "   unmoved" : "");
  }

  std::printf("\nworst rotation     %.3e rad\n", report.max_rotation_deviation);
  std::printf("worst translation  %.3e m\n", report.max_translation_deviation);
  std::printf("exact chain valid  %s\n", robot.validate().ok() ? "yes" : "no");

  // The recovered origins are exactly orthogonal; the file's were not. This is
  // the number that decides whether the offline pipeline is entitled to run.
  double worst_defect = 0.0;
  for (const auto& j : robot.joints()) {
    worst_defect = std::max(worst_defect, orthogonality_defect(j.origin.rotation()));
  }
  std::printf("orthogonality      %.3e (exact requires 0)\n", worst_defect);
  return 0;
}
