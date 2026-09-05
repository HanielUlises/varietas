// Solves the inverse kinematics of a URDF chain at one pose, exactly.
//
// Usage: urdf_solve <file.urdf> --xyz X Y Z [--rpy R P Y | --quat X Y Z W]
//                   [--tip L] [--root L] [--position-only] [--denominator D]
//
// This is the other half of what the library can do, and the half the emitter
// cannot reach. urdf_codegen adjoins the target to the coefficient field so
// that one basis answers every pose, which is what makes generated code
// possible and which also caps the problem at two adjoined parameters. Solving
// one pose at a time gives that up and gets orientation in return: the target
// is a constant, the field is Q, and the twelve equations of a full pose are
// no harder for Buchberger than the three of a position.
//
// The pose has to be exact, and that is not a formality. A target assembled by
// rounding sines and cosines into rationals has denominators near 2^52, and
// Buchberger over coefficients of that size does not finish: measured, not
// supposed. So the target is snapped to a nearby exact pose the same way the
// URDF's own geometry is, with the deviation reported rather than hidden.

#include <array>
#include <cmath>
#include <complex>
#include <cstdio>
#include <string>
#include <vector>

#include <urdf/model.h>

#include "varietas/codegen/rational.hpp"
#include "varietas/core/order/grevlex.hpp"
#include "varietas/core/quotient/quotient_basis.hpp"
#include "varietas/core/solve/spectral.hpp"
#include "varietas/kinematics/evaluate.hpp"
#include "varietas/kinematics/rationalize.hpp"
#include "varietas/urdf/rational_approximation.hpp"
#include "varietas/urdf/urdf_chain.hpp"

namespace {

using varietas::grevlex;
using varietas::rational;

struct options {
  std::string urdf;
  std::string tip;
  std::string root;
  double xyz[3] = {0.0, 0.0, 0.0};
  double quat[4] = {0.0, 0.0, 0.0, 1.0};  // x, y, z, w
  bool position_only = false;
  long denominator = 1000000;
};

void quaternion_from_rpy(double roll, double pitch, double yaw, double* out) {
  const double cr = std::cos(roll * 0.5), sr = std::sin(roll * 0.5);
  const double cp = std::cos(pitch * 0.5), sp = std::sin(pitch * 0.5);
  const double cy = std::cos(yaw * 0.5), sy = std::sin(yaw * 0.5);
  out[0] = sr * cp * cy - cr * sp * sy;  // x
  out[1] = cr * sp * cy + sr * cp * sy;  // y
  out[2] = cr * cp * sy - sr * sp * cy;  // z
  out[3] = cr * cp * cy + sr * sp * sy;  // w
}

bool parse(int argc, char** argv, options& out) {
  std::vector<std::string> positional;
  bool have_xyz = false;
  for (int i = 1; i < argc; ++i) {
    const std::string arg = argv[i];
    const auto number = [&](const char* what) -> double {
      if (i + 1 >= argc) {
        std::fprintf(stderr, "%s needs a value\n", what);
        return 0.0;
      }
      return std::atof(argv[++i]);
    };
    if (arg == "--xyz") {
      if (i + 3 >= argc) {
        std::fprintf(stderr, "--xyz needs three values\n");
        return false;
      }
      for (double& v : out.xyz) {
        v = std::atof(argv[++i]);
      }
      have_xyz = true;
    } else if (arg == "--rpy") {
      if (i + 3 >= argc) {
        std::fprintf(stderr, "--rpy needs three values\n");
        return false;
      }
      const double r = std::atof(argv[++i]);
      const double p = std::atof(argv[++i]);
      const double y = std::atof(argv[++i]);
      quaternion_from_rpy(r, p, y, out.quat);
    } else if (arg == "--quat") {
      if (i + 4 >= argc) {
        std::fprintf(stderr, "--quat needs four values (x y z w)\n");
        return false;
      }
      for (double& v : out.quat) {
        v = std::atof(argv[++i]);
      }
    } else if (arg == "--tip") {
      out.tip = argv[++i];
    } else if (arg == "--root") {
      out.root = argv[++i];
    } else if (arg == "--position-only") {
      out.position_only = true;
    } else if (arg == "--denominator") {
      out.denominator = static_cast<long>(number("--denominator"));
    } else if (!arg.empty() && arg[0] == '-') {
      std::fprintf(stderr, "unknown option %s\n", arg.c_str());
      return false;
    } else {
      positional.push_back(arg);
    }
  }
  if (positional.size() != 1 || !have_xyz) {
    return false;
  }
  out.urdf = positional[0];
  return true;
}

// The solve, for a chain of N actuated joints.
template <std::size_t N>
int run(const varietas::chain<rational>& robot,
        const varietas::rigid_transform<rational>& target, const options& opts) {
  using poly = varietas::polynomial<rational, N, grevlex>;

  if (N >= 6 && !opts.position_only) {
    std::fprintf(stderr,
                 "warning: a six-joint full pose has not been observed to finish here; five "
                 "joints takes seconds and six did not complete in eight minutes. Interrupt "
                 "if it runs away.\n");
  }

  const std::vector<poly> generators =
      opts.position_only
          ? varietas::position_ideal_generators<N, grevlex>(robot, target.translation())
          : varietas::pose_ideal_generators<N, grevlex>(robot, target);

  const auto quotient = varietas::standard_monomials(generators);
  if (!quotient.is_zero_dimensional) {
    std::fprintf(stderr,
                 "refused: the ideal is not zero-dimensional, so the arm reaches this target "
                 "along a curve of configurations rather than at finitely many\n");
    return 1;
  }
  if (quotient.dimension() == 0) {
    std::printf("solutions         0 (the ideal is the unit ideal: this pose is unreachable)\n");
    return 0;
  }

  std::printf("quotient dim      %zu\n", quotient.dimension());

  const auto solutions = varietas::solve_zero_dimensional(generators);
  if (!solutions.ok()) {
    std::fprintf(stderr, "refused: %s\n", varietas::to_string(solutions.status));
    return 1;
  }

  const auto real = solutions.real_points();
  std::printf("real solutions    %zu of %zu complex\n", real.size(), solutions.points.size());
  if (real.empty()) {
    return 0;
  }

  // The half-angle variables are what the variety is in; the angles are what a
  // caller drives. Both are printed, and the residual of the original equations
  // is what certifies the floating point step that produced them.
  std::printf("\n%-4s", "");
  for (std::size_t i = 0; i < N; ++i) {
    std::printf("%12s", robot.joints()[i].name.substr(0, 11).c_str());
  }
  std::printf("%14s%14s\n", "residual", "reached err");

  const auto numeric = varietas::chain_cast<double>(robot);
  for (std::size_t k = 0; k < real.size(); ++k) {
    std::array<std::complex<double>, N> point{};
    std::vector<double> angles(N, 0.0);
    for (std::size_t i = 0; i < N; ++i) {
      point[i] = real[k][i];
      angles[i] = varietas::angle_from_variable(real[k][i]);
    }

    // Independent of the algebra: drive the joints and see where the tool went.
    const auto frame = varietas::forward_kinematics(numeric, angles);
    double worst = 0.0;
    for (std::size_t i = 0; i < 3; ++i) {
      const double wanted =
          varietas::coefficient_traits<rational>::to_double(target.translation()[i]);
      worst = std::max(worst, std::fabs(frame.translation()[i] - wanted));
    }

    std::printf("%-4zu", k);
    for (std::size_t i = 0; i < N; ++i) {
      std::printf("%12.6f", angles[i]);
    }
    std::printf("%14.3e%14.3e\n", varietas::residual(generators, point), worst);
  }
  return 0;
}

}  // namespace

int main(int argc, char** argv) {
  options opts;
  if (!parse(argc, argv, opts)) {
    std::fprintf(stderr,
                 "usage: %s <file.urdf> --xyz X Y Z [--rpy R P Y | --quat X Y Z W]\n"
                 "          [--tip L] [--root L] [--position-only] [--denominator D]\n",
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

  varietas::chain<rational> robot;
  const auto report = varietas::urdf_import::chain_from_model(model, root, tip, robot);
  if (!report.ok()) {
    std::fprintf(stderr, "refused: %s (%s)\n", varietas::urdf_import::to_string(report.status),
                 report.detail.c_str());
    return 1;
  }
  robot = robot.fold_fixed_joints();

  // The target, made exact. A pose that is not exactly representable is moved
  // to one that is, and by how far is reported: an answer about a pose a
  // millimetre away is still an honest answer, but only if it says so.
  varietas::vector3<rational> translation;
  double worst_translation = 0.0;
  for (std::size_t i = 0; i < 3; ++i) {
    const auto snap = varietas::urdf_import::snap_scalar(opts.xyz[i], opts.denominator);
    translation[i] = snap.value;
    worst_translation = std::max(worst_translation, snap.deviation);
  }
  const auto rotation =
      varietas::urdf_import::snap_rotation(opts.quat[0], opts.quat[1], opts.quat[2],
                                           opts.quat[3], opts.denominator);
  const varietas::rigid_transform<rational> target(rotation.exact(), translation);

  const std::size_t dof = robot.degrees_of_freedom();
  std::printf("model             %s\n", model.getName().c_str());
  std::printf("chain             %s -> %s\n", root.c_str(), tip.c_str());
  std::printf("degrees of freedom %zu\n", dof);
  std::printf("asking for        %s\n", opts.position_only ? "position only" : "the full pose");
  std::printf("target moved      %.3e m, %.3e rad to become exact\n", worst_translation,
              rotation.deviation_radians);

  switch (dof) {
    case 1: return run<1>(robot, target, opts);
    case 2: return run<2>(robot, target, opts);
    case 3: return run<3>(robot, target, opts);
    case 4: return run<4>(robot, target, opts);
    case 5: return run<5>(robot, target, opts);
    case 6: return run<6>(robot, target, opts);
    default: break;
  }
  std::fprintf(stderr,
               "refused: this command is instantiated for up to six joints, and this chain "
               "has %zu\n",
               dof);
  return 1;
}
