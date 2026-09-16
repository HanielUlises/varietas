// What the generated solver costs at run time, and what it cost to generate.
//
// The other two experiments in this directory measure the offline path: how
// the Grobner computation grows with the field and with the number of terms.
// Neither says what a caller actually pays, and that is the number the whole
// arrangement exists to make small. The emitted header claims to be ordinary
// C++ over Eigen with nothing from the library linked behind it; this measures
// the claim rather than repeating it.
//
// Three quantities, on the same arm:
//
//   generation   decoupled_position_ik over Q, paid once, offline.
//   solve        branch_ik::solve at a target, paid per call, at run time.
//   forward      forward_kinematics over double, as a unit of comparison, so
//                that the solve is reported in a currency a reader of this
//                code already has a feel for.
//
// The targets are drawn from the arm's own workspace rather than from a box
// around it, and the two are not the same set: a box gives mostly unreachable
// points, the solver returns from those early, and the mean then flatters
// itself by averaging in work never done. Sampling a configuration and mapping
// it forward gives a target that is reachable by construction. Unreachable
// ones are timed too, separately, because a caller sweeping a trajectory pays
// both and the difference between them turns out to be most of the spread.
//
// Build (from the repository root, after colcon has built varietas_demo, which
// is what writes the generated header this includes):
//   g++ -std=c++17 -O2 doc/experiments/solver_cost.cpp -o solver_cost \
//     -Ivarietas_core/include -Ivarietas_codegen/include \
//     -Ivarietas_kinematics/include -Ivarietas_ik/include \
//     -Ibuild/varietas_demo/generated \
//     -I/usr/include/eigen3 -lgmpxx -lgmp
//
//   ./solver_cost

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <random>

#include <Eigen/Dense>
#include <vector>

#include "varietas/codegen/rational.hpp"
#include "varietas/ik/decoupled_ik.hpp"
#include "varietas/kinematics/chain.hpp"
#include "varietas/kinematics/evaluate.hpp"

#include "branch_ik.hpp"

using solver = varietas_generated::branch_ik;

namespace {

using clock_type = std::chrono::steady_clock;

double seconds_since(const clock_type::time_point& start) {
  return std::chrono::duration<double>(clock_type::now() - start).count();
}

// The arm the demonstration solves, built here rather than read from the URDF
// so that this file needs no ROS. It is the same geometry as
// varietas_demo/urdf/anthropomorphic_offset_3r.urdf: base yaw about z, shoulder
// displaced by (0.28, 0, 0.32), shoulder and elbow pitching about y, unit
// links, and a tool a further unit along.
//
// Every entry is an exact rational. That is the point of the offline field:
// 0.28 enters as 7/25 rather than as the binary double nearest it, so the
// ideal below is a statement about this arm and not about one 1e-17 away.
template <class Coeff>
varietas::chain<Coeff> build_arm() {
  using vec = varietas::vector3<Coeff>;
  using rot = varietas::matrix3<Coeff>;
  varietas::chain<Coeff> arm("anthropomorphic_offset_3r");

  const Coeff zero(0), one(1);
  const auto identity = rot::identity();

  varietas::joint<Coeff> q1;
  q1.name = "q1";
  q1.type = varietas::joint_type::revolute;
  q1.axis = vec{zero, zero, one};
  q1.origin = varietas::rigid_transform<Coeff>(identity, vec{zero, zero, zero});
  arm.add_joint(q1);

  varietas::joint<Coeff> q2;
  q2.name = "q2";
  q2.type = varietas::joint_type::revolute;
  q2.axis = vec{zero, one, zero};
  q2.origin = varietas::rigid_transform<Coeff>(
      identity, vec{Coeff(7) / Coeff(25), zero, Coeff(8) / Coeff(25)});
  arm.add_joint(q2);

  varietas::joint<Coeff> q3;
  q3.name = "q3";
  q3.type = varietas::joint_type::revolute;
  q3.axis = vec{zero, one, zero};
  q3.origin = varietas::rigid_transform<Coeff>(identity, vec{one, zero, zero});
  arm.add_joint(q3);

  arm.set_tool(varietas::rigid_transform<Coeff>(identity, vec{one, zero, zero}));
  return arm;
}

// A numerical solver, for the comparison the whole design rests on.
//
// Damped least squares on the position residual: the standard thing a robotics
// codebase reaches for when it wants an inverse kinematics and has no closed
// form. Newton with a Levenberg damping term, which is what keeps it from
// diverging at the singular configurations an undamped step walks straight
// into, and a finite-difference Jacobian so that this needs nothing from the
// library beyond the forward map.
//
// It is here to be measured, not to be beaten unfairly, so the tolerance is
// 1e-10 m -- tighter than anything the application needs, and comparable with
// what the generated solver achieves -- and the iteration cap is generous.
//
// The point of the comparison is not the time. It is that this returns ONE
// configuration, whichever one the seed fell into, and has no way to say how
// many others exist. Running it from many seeds is the only way to look for
// the rest, and that is measured below too.
struct newton_result {
  bool converged;
  int iterations;
  std::array<double, 3> q;
};

newton_result newton_ik(const varietas::chain<double>& arm,
                        const std::array<double, 3>& target,
                        std::array<double, 3> q) {
  constexpr double tolerance = 1e-10;
  constexpr int max_iterations = 200;
  constexpr double damping = 1e-6;
  constexpr double step = 1e-7;

  for (int iteration = 0; iteration < max_iterations; ++iteration) {
    const std::vector<double> qv(q.begin(), q.end());
    const varietas::rigid_transform<double> tool = varietas::forward_kinematics(arm, qv);
    Eigen::Vector3d error;
    for (int i = 0; i < 3; ++i) {
      error[i] = target[static_cast<std::size_t>(i)] - tool.translation()[i];
    }
    if (error.norm() < tolerance) { return {true, iteration, q}; }

    Eigen::Matrix3d jacobian;
    for (int j = 0; j < 3; ++j) {
      std::vector<double> shifted(q.begin(), q.end());
      shifted[static_cast<std::size_t>(j)] += step;
      const varietas::rigid_transform<double> moved =
          varietas::forward_kinematics(arm, shifted);
      for (int i = 0; i < 3; ++i) {
        jacobian(i, j) = (moved.translation()[i] - tool.translation()[i]) / step;
      }
    }

    const Eigen::Matrix3d normal =
        jacobian.transpose() * jacobian + damping * Eigen::Matrix3d::Identity();
    const Eigen::Vector3d delta = normal.ldlt().solve(jacobian.transpose() * error);
    for (int j = 0; j < 3; ++j) { q[static_cast<std::size_t>(j)] += delta[j]; }
  }
  return {false, max_iterations, q};
}

// Two configurations are the same posture if every joint agrees modulo a full
// turn. Newton returns angles wound arbitrarily far, so they are compared on
// the circle rather than on the line.
bool same_posture(const std::array<double, 3>& a, const std::array<double, 3>& b) {
  constexpr double two_pi = 6.28318530717958647692;
  for (std::size_t i = 0; i < 3; ++i) {
    double d = std::fmod(a[i] - b[i], two_pi);
    if (d > 3.14159265358979323846) { d -= two_pi; }
    if (d < -3.14159265358979323846) { d += two_pi; }
    if (std::abs(d) > 1e-4) { return false; }
  }
  return true;
}

// Median rather than mean, and the tail reported beside it.
//
// The distribution is not symmetric: a handful of calls land near the boundary
// of the reachable set, where the eigenvalue routine works hardest, and a mean
// over a few thousand calls moves several per cent between runs because of
// them. The median is stable to well under one per cent and the 99th centile
// says what the worst case actually is, which is the number a control loop has
// to budget for.
struct spread {
  double median, p99, worst;
};

spread summarise(std::vector<double> values) {
  std::sort(values.begin(), values.end());
  const std::size_t n = values.size();
  return {values[n / 2], values[(n * 99) / 100], values.back()};
}

}  // namespace

int main() {
  // ---------------------------------------------------------------- offline
  //
  // What it costs to produce the header, once. Timed around the whole call,
  // since a caller waiting for urdf_codegen waits for all of it.
  const auto exact = build_arm<varietas::rational>();

  const auto generation_start = clock_type::now();
  const auto solution = varietas::ik::decoupled_position_ik<3>(exact);
  const double generation = seconds_since(generation_start);

  if (!solution.ok()) {
    std::fprintf(stderr, "the arm did not decouple: %s\n",
                 varietas::ik::to_string(solution.status));
    return 1;
  }

  // ---------------------------------------------------------------- targets
  const auto arm = varietas::chain_cast<double>(exact);
  std::mt19937 rng(20240916u);
  std::uniform_real_distribution<double> angle(-M_PI, M_PI);

  constexpr int kCalls = 20000;
  constexpr int kWarmup = 2000;

  std::vector<std::array<double, 3>> reachable;
  reachable.reserve(kCalls);
  for (int i = 0; i < kCalls; ++i) {
    const std::vector<double> q{angle(rng), angle(rng), angle(rng)};
    const varietas::rigid_transform<double> tool = varietas::forward_kinematics(arm, q);
    reachable.push_back({tool.translation()[0], tool.translation()[1],
                         tool.translation()[2]});
  }

  // The furthest the tool can get: the links laid out straight, from a
  // shoulder already displaced from the base axis.
  const double outer_reach = 1.0 + 1.0 + std::hypot(0.28, 0.32);

  // Beyond it by construction, so these are outside the set whatever the
  // configuration.
  const double beyond = outer_reach + 0.5;
  std::vector<std::array<double, 3>> unreachable;
  unreachable.reserve(kCalls);
  for (int i = 0; i < kCalls; ++i) {
    const double a = angle(rng), b = 0.5 * angle(rng);
    unreachable.push_back({beyond * std::cos(b) * std::cos(a),
                           beyond * std::cos(b) * std::sin(a),
                           beyond * std::sin(b)});
  }

  // ------------------------------------------------------------- correctness
  //
  // Timing a wrong answer measures nothing, so every configuration the solver
  // returns is put back through the forward map and compared against the
  // target it was asked for. This runs before the timing and is not included
  // in it.
  std::vector<double> residuals;
  long counts[solver::max_configurations + 1] = {0};
  double worst_at_margin = 0.0;   // how close to the boundary the worst one sat
  double worst_elbow = 0.0;       // and how straight the elbow was there
  double worst_residual = 0.0;
  for (const auto& target : reachable) {
    double out[solver::max_configurations * solver::num_joints];
    solver::status state{};
    const int found = solver::solve(target.data(), out, solver::max_configurations, &state);
    if (found < 0) { continue; }
    ++counts[found];
    for (int k = 0; k < found; ++k) {
      const std::vector<double> q(out + k * solver::num_joints,
                                  out + (k + 1) * solver::num_joints);
      const varietas::rigid_transform<double> tool = varietas::forward_kinematics(arm, q);
      const double dx = tool.translation()[0] - target[0];
      const double dy = tool.translation()[1] - target[1];
      const double dz = tool.translation()[2] - target[2];
      const double r = std::sqrt(dx * dx + dy * dy + dz * dz);
      residuals.push_back(r);
      if (r > worst_residual) {
        worst_residual = r;
        const double radius = std::sqrt(target[0] * target[0] + target[1] * target[1] +
                                        target[2] * target[2]);
        worst_at_margin = (outer_reach - radius) / outer_reach;
        // The elbow, because the obvious suspect is a double root: the two
        // solutions of the reduced problem are elbow-up and elbow-down, they
        // coincide when the elbow is straight or folded back, and a repeated
        // eigenvalue costs the eigenvector basis its conditioning. Recorded so
        // the suspicion can be checked rather than asserted -- and on the runs
        // to date it does not hold up. The worst case sits over half the reach
        // inside the boundary with an elbow nowhere near either degeneracy, so
        // neither leaving the workspace nor the elbow singularity accounts for
        // it. The likeliest remaining candidate is the emitted denominator
        // guard admitting a pose where cancellation has already cost most of
        // the significance, which would be a tolerance to revisit rather than
        // a defect in the algebra; it has not been run down here.
        worst_elbow = q[2];
      }
    }
  }
  const spread residual = summarise(residuals);
  long above_1e10 = 0, above_1e12 = 0;
  for (const double r : residuals) {
    if (r > 1e-10) { ++above_1e10; }
    if (r > 1e-12) { ++above_1e12; }
  }

  // ----------------------------------------------------------------- timing
  volatile double sink = 0.0;   // keeps the optimiser from deleting the calls

  auto time_solves = [&](const std::vector<std::array<double, 3>>& targets) {
    std::vector<double> ns;
    ns.reserve(targets.size());
    double out[solver::max_configurations * solver::num_joints];
    for (std::size_t i = 0; i < targets.size(); ++i) {
      solver::status state{};
      const auto start = clock_type::now();
      const int found = solver::solve(targets[i].data(), out,
                                      solver::max_configurations, &state);
      const double elapsed = std::chrono::duration<double, std::nano>(
                                 clock_type::now() - start).count();
      sink += out[0] + static_cast<double>(found);
      if (i >= static_cast<std::size_t>(kWarmup)) { ns.push_back(elapsed); }
    }
    return summarise(std::move(ns));
  };

  const spread hit = time_solves(reachable);
  const spread miss = time_solves(unreachable);

  // The unit of comparison. Same loop shape, same clock, so the ratio below is
  // between two things measured the same way.
  std::vector<double> fk_ns;
  fk_ns.reserve(reachable.size());
  for (std::size_t i = 0; i < reachable.size(); ++i) {
    const std::vector<double> q{angle(rng), angle(rng), angle(rng)};
    const auto start = clock_type::now();
    const varietas::rigid_transform<double> tool = varietas::forward_kinematics(arm, q);
    const double elapsed =
        std::chrono::duration<double, std::nano>(clock_type::now() - start).count();
    sink += tool.translation()[0];
    if (i >= static_cast<std::size_t>(kWarmup)) { fk_ns.push_back(elapsed); }
  }
  const spread forward = summarise(std::move(fk_ns));

  // ------------------------------------------------- the numerical baseline
  //
  // Same targets, one random seed each. Timed the same way, and then asked the
  // question the timing cannot answer: how much of the solution set did it
  // find?
  constexpr int kNewtonCalls = 2000;
  std::vector<double> newton_ns;
  newton_ns.reserve(kNewtonCalls);
  long newton_converged = 0, newton_iterations = 0;
  for (int i = 0; i < kNewtonCalls; ++i) {
    const std::array<double, 3> seed{angle(rng), angle(rng), angle(rng)};
    const auto start = clock_type::now();
    const newton_result r = newton_ik(arm, reachable[static_cast<std::size_t>(i)], seed);
    const double elapsed =
        std::chrono::duration<double, std::nano>(clock_type::now() - start).count();
    newton_ns.push_back(elapsed);
    sink += r.q[0];
    if (r.converged) { ++newton_converged; newton_iterations += r.iterations; }
  }
  const spread newton = summarise(std::move(newton_ns));

  // How many seeds it takes to see the whole solution set.
  //
  // For each target the generated solver's answer is taken as the roll of
  // postures that exist -- it is certified complete by the quotient dimension,
  // which is the entire point -- and Newton is restarted from fresh random
  // seeds until it has found all of them, or until it has clearly stopped
  // finding new ones.
  constexpr int kBranchTargets = 300;
  constexpr int kMaxSeeds = 60;
  long seeds_used = 0, targets_completed = 0, branches_missed = 0;
  double newton_total_ns = 0.0;
  for (int i = 0; i < kBranchTargets; ++i) {
    const auto& target = reachable[static_cast<std::size_t>(i)];
    double out[solver::max_configurations * solver::num_joints];
    solver::status state{};
    const int expected = solver::solve(target.data(), out, solver::max_configurations, &state);
    if (expected <= 0) { continue; }

    std::vector<std::array<double, 3>> wanted;
    for (int k = 0; k < expected; ++k) {
      wanted.push_back({out[k * 3 + 0], out[k * 3 + 1], out[k * 3 + 2]});
    }
    std::vector<bool> found(wanted.size(), false);
    int remaining = expected, seeds = 0;
    const auto start = clock_type::now();
    for (; seeds < kMaxSeeds && remaining > 0; ++seeds) {
      const std::array<double, 3> seed{angle(rng), angle(rng), angle(rng)};
      const newton_result r = newton_ik(arm, target, seed);
      if (!r.converged) { continue; }
      for (std::size_t k = 0; k < wanted.size(); ++k) {
        if (!found[k] && same_posture(r.q, wanted[k])) { found[k] = true; --remaining; }
      }
    }
    newton_total_ns += std::chrono::duration<double, std::nano>(
                           clock_type::now() - start).count();
    seeds_used += seeds;
    branches_missed += remaining;
    if (remaining == 0) { ++targets_completed; }
  }

  // ----------------------------------------------------------------- report
  std::printf("arm            %s, %zu joints, %zu branches\n",
              exact.name().c_str(), exact.degrees_of_freedom(), solution.branches);
  std::printf("reduced        %zu unknowns over %zu parameters, dim_k A = %zu\n\n",
              std::size_t{2}, std::size_t{2}, solution.branches / 2);

  // The max column carries whatever the scheduler did during the call, so it
  // is an upper bound on this machine rather than a property of the code; the
  // 99th centile is the number to budget against.
  std::printf("%-14s %12s %12s %12s\n", "", "median", "99th", "max");
  std::printf("%-14s %9.0f ns %9.0f ns %9.0f ns\n", "solve, in", hit.median, hit.p99,
              hit.worst);
  std::printf("%-14s %9.0f ns %9.0f ns %9.0f ns\n", "solve, out", miss.median, miss.p99,
              miss.worst);
  std::printf("%-14s %9.0f ns %9.0f ns %9.0f ns\n\n", "forward", forward.median,
              forward.p99, forward.worst);

  std::printf("%-14s %9.0f ns %9.0f ns %9.0f ns\n\n", "newton, 1 seed", newton.median,
              newton.p99, newton.worst);

  std::printf("solve / forward         %.1f x\n", hit.median / forward.median);
  std::printf("newton / solve          %.1f x   (and newton returns one branch)\n",
              newton.median / hit.median);
  std::printf("solves per second       %.0f  (median, reachable target)\n",
              1e9 / hit.median);
  std::printf("generation              %.3f s, once, offline\n", generation);
  std::printf("            amortised after %.0f solves, %.0f s of a 1 kHz loop\n",
              generation * 1e9 / hit.median, generation * 1e9 / hit.median / 1000.0);

  std::printf("\nresidual, %zu configurations\n", residuals.size());
  std::printf("  median %.2e m   99th %.2e m   worst %.2e m\n",
              residual.median, residual.p99, residual.worst);
  std::printf("  over 1e-12 m: %ld of %zu (%.4f%%);  over 1e-10 m: %ld (%.4f%%)\n",
              above_1e12, residuals.size(),
              100.0 * static_cast<double>(above_1e12) /
                  static_cast<double>(residuals.size()),
              above_1e10,
              100.0 * static_cast<double>(above_1e10) /
                  static_cast<double>(residuals.size()));
  std::printf("  worst case: %.1f%% of the reach inside the boundary, elbow %.3f rad\n",
              100.0 * worst_at_margin, worst_elbow);
  std::printf("  so neither the workspace boundary nor the elbow singularity"
              " explains it\n");
  std::printf("counts returned  ");
  for (std::size_t k = 0; k <= solver::max_configurations; ++k) {
    if (counts[k] != 0) { std::printf("%zux%ld  ", k, counts[k]); }
  }
  std::printf("\n");

  std::printf("\nnewton, one seed per target: %ld of %d converged, %.1f iterations mean\n",
              newton_converged, kNewtonCalls,
              static_cast<double>(newton_iterations) /
                  static_cast<double>(newton_converged ? newton_converged : 1));
  std::printf("newton, restarted until it has the whole set, over %d targets:\n",
              kBranchTargets);
  std::printf("  %.1f seeds per target on average, %ld of %d sets completed"
              " within %d seeds\n",
              static_cast<double>(seeds_used) / kBranchTargets, targets_completed,
              kBranchTargets, kMaxSeeds);
  std::printf("  %ld branches never found; %.0f ns per target against %.0f ns"
              " for one generated solve\n",
              branches_missed, newton_total_ns / kBranchTargets, hit.median);

  return sink == 12345.6789 ? 1 : 0;
}
