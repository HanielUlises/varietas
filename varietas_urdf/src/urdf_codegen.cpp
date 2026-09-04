// Reads a URDF and writes a header that solves its inverse kinematics.
//
// Usage: urdf_codegen <file.urdf> <output.hpp> [--tip L] [--root L]
//                     [--coords xy|xz|yz|xyz] [--name N] [--namespace NS]
//                     [--matrices-only]
//
// This is the whole pipeline as a command. The URDF is read and made exact,
// the chain is posed with the target adjoined to the coefficient field, one
// Grobner basis is computed over Q(p), and the action matrices that basis
// produces are written out as C++. Nothing here knows a pose: the header it
// leaves behind answers every pose the basis describes, and the poses it does
// not describe it refuses.
//
// The expensive step is the Grobner basis over Q(p), and it happens once, here,
// rather than once per request in a control loop. That trade is the reason the
// offline half of this library exists.

#include <array>
#include <cstddef>
#include <cstdio>
#include <fstream>
#include <string>
#include <vector>

#include <urdf/model.h>

#include "varietas/codegen/emit.hpp"
#include "varietas/ik/parametric_ik.hpp"
#include "varietas/urdf/urdf_chain.hpp"

namespace {

struct options {
  std::string urdf;
  std::string output;
  std::string tip;
  std::string root;
  std::string coords;
  std::string name = "urdf_ik";
  std::string name_space = "varietas_generated";
  varietas::codegen::runtime_kind runtime = varietas::codegen::runtime_kind::eigen;
};

bool parse(int argc, char** argv, options& out) {
  std::vector<std::string> positional;
  for (int i = 1; i < argc; ++i) {
    const std::string arg = argv[i];
    const auto value = [&](const char* what) -> std::string {
      if (i + 1 >= argc) {
        std::fprintf(stderr, "%s needs a value\n", what);
        return std::string();
      }
      return argv[++i];
    };
    if (arg == "--tip") {
      out.tip = value("--tip");
    } else if (arg == "--root") {
      out.root = value("--root");
    } else if (arg == "--coords") {
      out.coords = value("--coords");
    } else if (arg == "--name") {
      out.name = value("--name");
    } else if (arg == "--namespace") {
      out.name_space = value("--namespace");
    } else if (arg == "--matrices-only") {
      out.runtime = varietas::codegen::runtime_kind::matrices_only;
    } else if (!arg.empty() && arg[0] == '-') {
      std::fprintf(stderr, "unknown option %s\n", arg.c_str());
      return false;
    } else {
      positional.push_back(arg);
    }
  }
  if (positional.size() != 2) {
    return false;
  }
  out.urdf = positional[0];
  out.output = positional[1];
  return true;
}

// "xy" -> {0, 1}. The order is kept, because it is the order the generated
// header will read its pose argument in.
bool parse_coordinates(const std::string& text, std::vector<std::size_t>& out) {
  for (const char c : text) {
    const std::size_t index = c == 'x' ? 0 : c == 'y' ? 1 : c == 'z' ? 2 : 3;
    if (index == 3) {
      std::fprintf(stderr, "--coords takes only the letters x, y and z\n");
      return false;
    }
    out.push_back(index);
  }
  return !out.empty() && out.size() <= 3;
}

bool write(const std::string& path, const std::string& text) {
  std::ofstream out(path);
  if (!out) {
    std::fprintf(stderr, "cannot write %s\n", path.c_str());
    return false;
  }
  out << text;
  if (!out) {
    std::fprintf(stderr, "write to %s failed\n", path.c_str());
    return false;
  }
  return true;
}

// The solve and the emission, for one (unknowns, parameters) pair. The two
// counts are template parameters because the polynomial ring's variable count
// is, so the dispatch below enumerates the handful of pairs that can occur for
// a tool position.
template <std::size_t N, std::size_t P>
int run(const varietas::chain<varietas::rational>& robot,
        const std::vector<std::size_t>& coordinates, const options& opts) {
  std::array<std::size_t, P> selected{};
  for (std::size_t k = 0; k < P; ++k) {
    selected[k] = coordinates[k];
  }

  const auto result = varietas::ik::parametric_position_ik<N, P>(robot, selected);
  if (!result.ok()) {
    std::fprintf(stderr, "refused: %s\n", varietas::ik::to_string(result.status));
    if (result.status ==
        varietas::ik::parametric_ik_status::dropped_coordinate_is_not_identically_zero) {
      std::fprintf(stderr,
                   "  coordinate %zu is moved by this arm, so it cannot be left out of "
                   "--coords\n",
                   result.offending_coordinate);
    }
    return 1;
  }

  std::printf("branches         %zu\n", result.branches);
  std::printf("quotient dim     %zu\n", result.solution.dimension());

  varietas::codegen::emit_options emit_options;
  emit_options.name = opts.name;
  emit_options.name_space = opts.name_space;
  emit_options.runtime = opts.runtime;
  emit_options.source_note = "Position inverse kinematics of " + robot.name() +
                             ", solved once over Q(pose) by varietas_ik from " + opts.urdf +
                             ". The unknowns are t = tan(q/2) for revolute joints.";

  if (!write(opts.output, varietas::codegen::emit(result.solution, emit_options))) {
    return 1;
  }
  std::printf("wrote            %s\n", opts.output.c_str());
  return 0;
}

}  // namespace

int main(int argc, char** argv) {
  options opts;
  if (!parse(argc, argv, opts)) {
    std::fprintf(stderr,
                 "usage: %s <file.urdf> <output.hpp> [--tip L] [--root L]\n"
                 "          [--coords xy|xz|yz|xyz] [--name N] [--namespace NS]\n"
                 "          [--matrices-only]\n",
                 argv[0]);
    return 2;
  }

  urdf::Model model;
  if (!model.initFile(opts.urdf)) {
    std::fprintf(stderr, "could not parse %s\n", opts.urdf.c_str());
    return 1;
  }

  const std::string tip =
      opts.tip.empty() ? varietas::urdf_import::sole_tip_link(model) : opts.tip;
  const std::string root = opts.root.empty() ? model.getRoot()->name : opts.root;
  if (tip.empty()) {
    std::fprintf(stderr, "%s is branched; name the tip link with --tip\n", opts.urdf.c_str());
    return 1;
  }

  varietas::chain<varietas::rational> robot;
  const auto report = varietas::urdf_import::chain_from_model(model, root, tip, robot);
  if (!report.ok()) {
    std::fprintf(stderr, "refused: %s (%s)\n", varietas::urdf_import::to_string(report.status),
                 report.detail.c_str());
    return 1;
  }

  const std::size_t dof = robot.degrees_of_freedom();
  std::printf("model            %s\n", model.getName().c_str());
  std::printf("chain            %s -> %s\n", root.c_str(), tip.c_str());
  std::printf("degrees of freedom %zu\n", dof);

  // Defaulted from the joint count rather than guessed from the geometry: an
  // arm with two joints can only be asked for a point in a plane, and one with
  // three for a point in space. Anything else has to be said out loud.
  std::vector<std::size_t> coordinates;
  if (opts.coords.empty()) {
    if (dof == 2) {
      coordinates = {0, 1};
    } else if (dof == 3) {
      coordinates = {0, 1, 2};
    } else {
      std::fprintf(stderr,
                   "a %zu-joint arm needs --coords: a tool position constrains at most three "
                   "unknowns, so only one, two or three joints can be solved for this way\n",
                   dof);
      return 1;
    }
  } else if (!parse_coordinates(opts.coords, coordinates)) {
    return 1;
  }

  std::printf("coordinates      %zu\n", coordinates.size());

  // Only P >= N can give a finite solution set, and P is at most three; the
  // bridge says so itself, but there is no point instantiating the rest.
  const std::size_t p = coordinates.size();
  if (dof == 1 && p == 1) return run<1, 1>(robot, coordinates, opts);
  if (dof == 1 && p == 2) return run<1, 2>(robot, coordinates, opts);
  if (dof == 1 && p == 3) return run<1, 3>(robot, coordinates, opts);
  if (dof == 2 && p == 2) return run<2, 2>(robot, coordinates, opts);
  if (dof == 2 && p == 3) return run<2, 3>(robot, coordinates, opts);
  if (dof == 3 && p == 3) return run<3, 3>(robot, coordinates, opts);

  std::fprintf(stderr,
               "refused: %zu joints against %zu coordinates cannot give a finite solution "
               "set; a tool position constrains at most three unknowns\n",
               dof, p);
  return 1;
}
