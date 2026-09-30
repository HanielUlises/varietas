// The three-parameter position solve, on several arms, by both routes.
//
// field_cost.cpp established that the anthropomorphic arm over Q(x, y, z) did
// not finish with the subresultant gcd, and gcd_cost.cpp that the remainder
// sequence was why. This measures what happens after the change, across arms
// of different shapes, including arms the decoupling cannot touch: what the
// symbolic solve over Q(x, y, z) costs with the modular gcd, what the
// reconstruction from fixed poses costs, and whether either answer is right.
//
// Right is checked exactly rather than numerically. A configuration is chosen
// whose half-angle tangents are rational, so that its tool position is
// rational too; the parametric action matrices are specialised at that
// position, over Q; and the vector of standard monomials evaluated at the
// configuration has to be a common eigenvector of their transposes, with the
// configuration's own coordinates as eigenvalues. That is the defining
// property of the matrices the header stores, so an error anywhere in the
// basis, the normal forms or the cancellation would show here as an inequality
// between two rationals, with no tolerance to hide in.
//
// Build (from the repository root):
//   g++ -std=c++17 -O2 doc/experiments/spatial_cost.cpp -o spatial_cost \
//     -Ivarietas_core/include -Ivarietas_codegen/include \
//     -Ivarietas_kinematics/include -Ivarietas_ik/include \
//     -I/usr/include/eigen3 -lgmpxx -lgmp
//
//   ./spatial_cost                  # every arm, reconstructed
//   ./spatial_cost <name>           # one of them
//   ./spatial_cost --symbolic ...   # the symbolic solve as well; on the arms
//                                   # with offsets it does not finish

#include <array>
#include <chrono>
#include <cstddef>
#include <cstdio>
#include <string>
#include <vector>

#include "varietas/codegen/modular_gcd.hpp"
#include "varietas/codegen/rational.hpp"
#include "varietas/ik/decoupled_ik.hpp"
#include "varietas/ik/parametric_ik.hpp"
#include "varietas/ik/reconstructed_ik.hpp"
#include "varietas/kinematics/chain.hpp"
#include "varietas/kinematics/rationalize.hpp"

namespace {

using varietas::chain;
using varietas::grevlex;
using varietas::make_rational;
using varietas::rational;
using varietas::revolute_joint;
using varietas::rigid_transform;
using varietas::vector3;

using clock_type = std::chrono::steady_clock;

rigid_transform<rational> at(rational x, rational y, rational z) {
  return rigid_transform<rational>::translation_only(vector3<rational>(x, y, z));
}

rational q(long n, long d = 1) { return make_rational(n, d); }

chain<rational> three_joints(const char* name, const std::array<vector3<rational>, 3>& axes,
                             const std::array<rigid_transform<rational>, 3>& origins,
                             const rigid_transform<rational>& tool) {
  chain<rational> robot(name);
  for (std::size_t i = 0; i < 3; ++i) {
    robot.add_joint(
        revolute_joint<rational>("q" + std::to_string(i + 1), axes[i], origins[i]));
  }
  robot.set_tool(tool);
  return robot;
}

const vector3<rational> ex = vector3<rational>::unit(0);
const vector3<rational> ey = vector3<rational>::unit(1);
const vector3<rational> ez = vector3<rational>::unit(2);

struct arm {
  const char* name;
  const char* description;
  chain<rational> robot;
};

// Five arms, ordered roughly by how far each is from the textbook case. Every
// axis is a rational unit vector, which is what an exact chain needs; the last
// arm's axes are points of the unit sphere with rational coordinates, placed so
// that the three are pairwise skew: no two parallel and no two meeting.
std::vector<arm> arms() {
  const auto I = rigid_transform<rational>::identity();
  std::vector<arm> out;

  out.push_back({"anthropomorphic",
                 "base yaw about z, shoulder and elbow about y, unit links",
                 three_joints("anthropomorphic", {ez, ey, ey}, {I, I, at(q(1), q(0), q(0))},
                              at(q(1), q(0), q(0)))});

  out.push_back({"offset-shoulder",
                 "as above, shoulder displaced (7/25, 0, 8/25) from the base axis",
                 three_joints("offset_shoulder", {ez, ey, ey},
                              {I, at(q(7, 25), q(0), q(8, 25)), at(q(1), q(0), q(0))},
                              at(q(1), q(0), q(0)))});

  out.push_back({"lateral-offset",
                 "elbow displaced 1/4 sideways, out of the plane of the arm",
                 three_joints("lateral_offset", {ez, ey, ey}, {I, I, at(q(1), q(1, 4), q(0))},
                              at(q(1), q(0), q(0)))});

  out.push_back({"off-axis-base",
                 "base joint displaced off its own axis; admits no sweep",
                 three_joints("off_axis_base", {ez, ey, ey},
                              {at(q(1), q(0), q(0)), I, at(q(1), q(0), q(0))},
                              at(q(1), q(0), q(0)))});

  out.push_back({"general",
                 "axes z, (3/5,4/5,0), (2/7,3/7,6/7); offsets in general position",
                 three_joints("general_3r",
                              {ez, vector3<rational>(q(3, 5), q(4, 5), q(0)),
                               vector3<rational>(q(2, 7), q(3, 7), q(6, 7))},
                              {I, at(q(1, 2), q(0), q(1, 3)), at(q(1), q(1, 5), q(0))},
                              at(q(1), q(0), q(1, 2)))});
  return out;
}

// The exact check described at the top.
//
// t is a configuration in half-angle coordinates. The tool position there is
// the rational forward map at t; each parametric matrix entry is a rational
// function of the position and is evaluated there; and w_k = b_k(t) for the
// standard monomials b_k. Then for every unknown i, M_i^T w = t_i w exactly.
bool verify_at(const chain<rational>& robot, const varietas::codegen::parametric_solution<3, 3>& s,
               const std::array<rational, 3>& t, std::string& why) {
  const chain<rational> folded = robot.fold_fixed_joints();
  const auto map = varietas::rational_forward_kinematics<3, grevlex>(folded);
  const rational denominator = map.denominator().evaluate(t);
  if (sgn(denominator) == 0) {
    why = "the configuration is on the pole of the half-angle map";
    return false;
  }
  std::array<rational, 3> pose;
  for (std::size_t k = 0; k < 3; ++k) {
    pose[k] = map.translation(k).evaluate(t) / denominator;
  }

  const std::size_t d = s.dimension();
  std::vector<rational> w(d);
  for (std::size_t k = 0; k < d; ++k) {
    w[k] = s.quotient.monomials[k].evaluate<rational>(t);
  }

  for (std::size_t i = 0; i < 3; ++i) {
    for (std::size_t j = 0; j < d; ++j) {
      rational lhs = 0;
      for (std::size_t k = 0; k < d; ++k) {
        const auto& entry = s.action[i](k, j);
        if (entry.is_zero()) {
          continue;
        }
        if (sgn(entry.denominator().evaluate(pose)) == 0) {
          why = "the chosen pose lies on a pole of the emitted matrices";
          return false;
        }
        lhs += entry.evaluate(pose) * w[k];
      }
      if (lhs != t[i] * w[j]) {
        why = "M^T w != t w for unknown " + std::to_string(i);
        return false;
      }
    }
  }
  return true;
}

// Three configurations, none special, all with rational half-angle tangents.
int configuration_checks(const arm& a, const varietas::codegen::parametric_solution<3, 3>& s) {
  const std::array<std::array<rational, 3>, 3> configurations{{
      {q(1, 2), q(-1, 3), q(2, 5)},
      {q(-3, 7), q(5, 4), q(1, 9)},
      {q(2), q(1, 6), q(-4, 3)},
  }};
  int verified = 0;
  for (const auto& t : configurations) {
    std::string why;
    if (verify_at(a.robot, s, t, why)) {
      ++verified;
    } else {
      std::printf("    check failed: %s\n", why.c_str());
    }
  }
  return verified;
}

void describe(const varietas::codegen::parametric_solution<3, 3>& s) {
  std::size_t entries = 0;
  std::size_t largest_entry = 0;
  for (const auto& m : s.action) {
    for (const auto& e : m.entries) {
      if (e.is_zero()) {
        continue;
      }
      ++entries;
      largest_entry = std::max({largest_entry, e.numerator().size(), e.denominator().size()});
    }
  }
  std::printf("    %zu nonzero matrix entries, the largest %zu terms\n", entries, largest_entry);
}

void run(const arm& a, bool symbolic) {
  std::printf("%-16s %s\n", a.name, a.description);

  // The decoupled solve, where there is one: the only route the arm had
  // before, and the cheapest where it applies.
  const auto decoupled_start = clock_type::now();
  const auto decoupled = varietas::ik::decoupled_position_ik<3>(a.robot);
  const double decoupled_seconds =
      std::chrono::duration<double>(clock_type::now() - decoupled_start).count();
  if (decoupled.ok()) {
    std::printf("    decoupled over Q(r,z): %.3f s, %zu branches\n", decoupled_seconds,
                decoupled.branches);
  } else {
    std::printf("    decoupling refused: %s\n", varietas::ik::to_string(decoupled.status));
  }

  varietas::ik::reconstruction_report report;
  const auto start = clock_type::now();
  const auto reconstructed =
      varietas::ik::reconstructed_position_ik<3, 3>(a.robot, {0, 1, 2}, &report);
  const double seconds = std::chrono::duration<double>(clock_type::now() - start).count();
  if (!reconstructed.ok()) {
    std::printf("    reconstruction refused: %s\n",
                varietas::ik::to_string(reconstructed.status));
  } else {
    const auto& r = report.sampling;
    std::printf("    reconstructed: %.2f s, %zu branches; %zu primes, %zu samples (%zu on the "
                "degree line), degrees %d/%d, largest system %zu unknowns\n",
                seconds, reconstructed.branches, r.primes, r.samples, r.line_samples,
                r.largest_numerator_degree, r.largest_denominator_degree, r.largest_system);
    std::printf("    of which %.2f s sampling and %.2f s in linear systems; exact at %zu rational "
                "poses, and at %d of 3 rational configurations\n",
                r.sampling_seconds, r.solving_seconds, report.exact_checks,
                configuration_checks(a, reconstructed.solution));
    describe(reconstructed.solution);
  }

  if (symbolic) {
    auto& gcd = varietas::modular_gcd_counters();
    gcd = varietas::modular_gcd_statistics{};
    const auto symbolic_start = clock_type::now();
    const auto result = varietas::ik::parametric_position_ik<3, 3>(a.robot, {0, 1, 2});
    const double symbolic_seconds =
        std::chrono::duration<double>(clock_type::now() - symbolic_start).count();
    if (!result.ok()) {
      std::printf("    symbolic refused: %s\n", varietas::ik::to_string(result.status));
    } else {
      std::printf("    symbolic over Q(x,y,z): %.2f s; gcd %zu calls, %.2f s (%.0f%%), largest "
                  "operand %zu terms; exact at %d of 3 rational configurations\n",
                  symbolic_seconds, gcd.calls, gcd.seconds,
                  100.0 * gcd.seconds / symbolic_seconds, gcd.largest_operand,
                  configuration_checks(a, result.solution));
      if (reconstructed.ok()) {
        std::size_t differing = 0;
        for (std::size_t i = 0; i < 3; ++i) {
          for (std::size_t e = 0; e < result.solution.action[i].entries.size(); ++e) {
            differing += result.solution.action[i].entries[e] !=
                         reconstructed.solution.action[i].entries[e];
          }
        }
        std::printf("    symbolic and reconstructed matrices differ in %zu entries\n", differing);
      }
    }
  }
  std::printf("\n");
  std::fflush(stdout);
}

}  // namespace

int main(int argc, char** argv) {
  bool symbolic = false;
  std::string which = "all";
  for (int i = 1; i < argc; ++i) {
    const std::string arg = argv[i];
    if (arg == "--symbolic") {
      symbolic = true;
    } else {
      which = arg;
    }
  }
  for (const auto& a : arms()) {
    if (which == "all" || which == a.name) {
      run(a, symbolic);
    }
  }
  return 0;
}
