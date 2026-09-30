// Two generated solvers for one arm: the full solve over Q(x, y, z), and the
// decoupled one over Q(r, z) with the base joint put back by an arctangent.
//
// The arm is the demonstration's, with its shoulder displaced from the base
// axis, which is the arm solver_cost.cpp measures. Until the pose could be
// reconstructed, the decoupled header was the only one it had: the symbolic
// solve over Q(x, y, z) ran for twenty-seven minutes and two gigabytes on it
// without finishing, and reconstruction takes a few seconds. The decoupling
// exists only for arms whose base placement commutes with the base rotation;
// the full solve exists for any arm the counting admits. The question here is
// what a caller gives up by using the full one where both are available: the
// run-time cost per call, the accuracy, and whether the two return the same
// configurations, which they must, both being complete.
//
// Build (from the repository root, after colcon has built varietas_demo, which
// writes the decoupled header, and after writing the full one):
//   ros2 run varietas_urdf urdf_codegen varietas_demo/urdf/anthropomorphic_offset_3r.urdf \
//     full_offset_3r.hpp --name full_offset_3r_ik --reconstruct
//   g++ -std=c++17 -O2 doc/experiments/full_solver_cost.cpp -o full_solver_cost \
//     -I. -Ibuild/varietas_demo/generated -I/usr/include/eigen3
//
//   ./full_solver_cost

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <random>
#include <vector>

#include "branch_ik.hpp"
#include "full_offset_3r.hpp"

using full = varietas_generated::full_offset_3r_ik;
using decoupled = varietas_generated::branch_ik;

namespace {

using clock_type = std::chrono::steady_clock;
using config = std::array<double, 3>;

// The arm, over double: base yaw about z, a shoulder at (0.28, 0, 0.32) and an
// elbow a unit further, both pitching about y, and a tool a unit beyond. It is
// written out rather than taken from the library so that the two headers are
// the only generated code in play.
std::array<double, 3> forward(const config& q) {
  const double reach = 0.28 + std::cos(q[1]) + std::cos(q[1] + q[2]);
  const double height = 0.32 - (std::sin(q[1]) + std::sin(q[1] + q[2]));
  return {std::cos(q[0]) * reach, std::sin(q[0]) * reach, height};
}

double residual(const config& q, const std::array<double, 3>& target) {
  const auto p = forward(q);
  return std::sqrt((p[0] - target[0]) * (p[0] - target[0]) +
                   (p[1] - target[1]) * (p[1] - target[1]) +
                   (p[2] - target[2]) * (p[2] - target[2]));
}

bool same_posture(const config& a, const config& b) {
  constexpr double two_pi = 6.28318530717958647692;
  for (std::size_t i = 0; i < 3; ++i) {
    double d = std::fmod(a[i] - b[i], two_pi);
    if (d > M_PI) d -= two_pi;
    if (d < -M_PI) d += two_pi;
    if (std::abs(d) > 1e-6) return false;
  }
  return true;
}

// The full solver returns t = tan(q/2); the decoupled wrapper returns angles.
int solve_full(const std::array<double, 3>& target, std::vector<config>& out) {
  double raw[full::dimension * full::num_unknowns];
  const int n = full::solve(target.data(), raw, static_cast<int>(full::dimension));
  out.clear();
  for (int k = 0; k < n; ++k) {
    out.push_back({2.0 * std::atan(raw[3 * k]), 2.0 * std::atan(raw[3 * k + 1]),
                   2.0 * std::atan(raw[3 * k + 2])});
  }
  return n;
}

int solve_decoupled(const std::array<double, 3>& target, std::vector<config>& out) {
  double raw[decoupled::max_configurations * decoupled::num_joints];
  const int n = decoupled::solve(target.data(), raw,
                                 static_cast<int>(decoupled::max_configurations));
  out.clear();
  for (int k = 0; k < n; ++k) {
    out.push_back({raw[3 * k], raw[3 * k + 1], raw[3 * k + 2]});
  }
  return n;
}

struct spread {
  double median, p99;
};

spread summarise(std::vector<double> v) {
  std::sort(v.begin(), v.end());
  return {v[v.size() / 2], v[(v.size() * 99) / 100]};
}

template <class Solve>
spread time_calls(const std::vector<std::array<double, 3>>& targets, Solve solve) {
  volatile double sink = 0.0;
  double raw[16];
  std::vector<double> ns;
  ns.reserve(targets.size());
  for (std::size_t i = 0; i < targets.size(); ++i) {
    const auto start = clock_type::now();
    const int n = solve(targets[i].data(), raw);
    ns.push_back(std::chrono::duration<double, std::nano>(clock_type::now() - start).count());
    sink = sink + raw[0] + n;
  }
  ns.erase(ns.begin(), ns.begin() + static_cast<std::ptrdiff_t>(targets.size() / 10));
  return summarise(std::move(ns));
}

}  // namespace

int main() {
  std::mt19937 rng(20260930u);
  std::uniform_real_distribution<double> angle(-M_PI, M_PI);

  constexpr int kTargets = 20000;
  std::vector<std::array<double, 3>> targets;
  targets.reserve(kTargets);
  for (int i = 0; i < kTargets; ++i) {
    targets.push_back(forward({angle(rng), angle(rng), angle(rng)}));
  }

  // Correctness and agreement, before any timing.
  std::vector<double> residual_full, residual_decoupled;
  long agree = 0, disagree = 0;
  std::array<long, 5> counts_full{}, counts_decoupled{};
  std::vector<config> a, b;
  for (const auto& target : targets) {
    const int nf = solve_full(target, a);
    const int nd = solve_decoupled(target, b);
    if (nf >= 0 && nf <= 4) ++counts_full[static_cast<std::size_t>(nf)];
    if (nd >= 0 && nd <= 4) ++counts_decoupled[static_cast<std::size_t>(nd)];
    for (const auto& q : a) residual_full.push_back(residual(q, target));
    for (const auto& q : b) residual_decoupled.push_back(residual(q, target));

    bool same = nf == nd;
    for (std::size_t i = 0; same && i < a.size(); ++i) {
      same = std::any_of(b.begin(), b.end(), [&](const config& q) { return same_posture(a[i], q); });
    }
    same ? ++agree : ++disagree;
  }

  const spread t_full = time_calls(targets, [](const double* p, double* out) {
    return full::solve(p, out, static_cast<int>(full::dimension));
  });
  const spread t_decoupled = time_calls(targets, [](const double* p, double* out) {
    return decoupled::solve(p, out, static_cast<int>(decoupled::max_configurations));
  });

  const spread r_full = summarise(residual_full);
  const spread r_decoupled = summarise(residual_decoupled);
  const auto worst = [](const std::vector<double>& v) { return *std::max_element(v.begin(), v.end()); };

  std::printf("%-22s %12s %12s\n", "", "full Q(x,y,z)", "decoupled");
  std::printf("%-22s %9.0f ns %9.0f ns\n", "solve, median", t_full.median, t_decoupled.median);
  std::printf("%-22s %9.0f ns %9.0f ns\n", "solve, 99th centile", t_full.p99, t_decoupled.p99);
  std::printf("%-22s %12.1e %12.1e\n", "residual, median (m)", r_full.median, r_decoupled.median);
  std::printf("%-22s %12.1e %12.1e\n", "residual, 99th (m)", r_full.p99, r_decoupled.p99);
  std::printf("%-22s %12.1e %12.1e\n", "residual, worst (m)", worst(residual_full),
              worst(residual_decoupled));
  std::printf("%-22s", "configurations returned");
  for (std::size_t k = 0; k < 5; ++k) std::printf("  %zu:%ld/%ld", k, counts_full[k], counts_decoupled[k]);
  std::printf("\nsame configurations on %ld of %d targets (%ld differ)\n", agree, kTargets, disagree);
  return 0;
}
