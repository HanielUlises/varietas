#ifndef VARIETAS_IK_RUNTIME_HPP
#define VARIETAS_IK_RUNTIME_HPP

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <Eigen/LU>

#include "varietas/codegen/parametric_solution.hpp"
#include "varietas/codegen/rational_function.hpp"
#include "varietas/ik/decoupled_ik.hpp"
#include "varietas/ik/reconstructed_ik.hpp"
#include "varietas/ik/spherical_wrist.hpp"
#include "varietas/kinematics/chain.hpp"

// The solvers the emitter writes out, evaluated in process instead.
//
// A generated header is the right thing for a fixed arm: it is compiled once
// and needs nothing at run time. A planner that is handed an arm it has never
// seen, as MoveIt hands a kinematics plugin, cannot compile anything, and for
// that case the same solution is evaluated here, from the same objects the
// emitter consumes, doing the same arithmetic in the same order: the
// separating combination with the same coefficients, the same denominator
// guard, the same Newton steps on the equations as posed, and for six joints
// the same split at the wrist. The tests hold the two to the same answers.
namespace varietas {
namespace ik {

// Joint ranges, and the fitting of an angle into one by whole turns.
//
// A solver returns each angle in (-pi, pi], but a joint's range need not be
// that interval: it can be [0, 2 pi), or wider than a turn, or narrower. An
// angle that lies outside its range may still be reachable a whole turn away,
// and only when no whole turn brings it inside is the configuration out of
// the joint's reach.
struct joint_limits : detail::emitted_limits {
  std::size_t size() const noexcept { return limited.size(); }

  // From the chain's own joints, actuated ones only, in order.
  template <class Coeff>
  static joint_limits of(const chain<Coeff>& robot) {
    joint_limits out;
    static_cast<detail::emitted_limits&>(out) = detail::limits_of(robot);
    return out;
  }

  // Moves each angle of `q` by whole turns into its joint's range, keeping it
  // where it is if it is already inside. False if some joint cannot be fitted,
  // in which case `q` is left partly moved and should be discarded.
  bool fit(double* q) const {
    constexpr double turn = 6.28318530717958647692;
    constexpr double slack = 1e-9;
    for (std::size_t i = 0; i < limited.size(); ++i) {
      if (!limited[i]) {
        continue;
      }
      double v = q[i];
      if (v < lower[i]) {
        v += turn * std::ceil((lower[i] - v - slack) / turn);
      } else if (v > upper[i]) {
        v -= turn * std::ceil((v - upper[i] - slack) / turn);
      }
      if (v < lower[i] - slack || v > upper[i] + slack) {
        return false;
      }
      q[i] = v;
    }
    return true;
  }

  // Fits every configuration and keeps those that fit, in order. Returns how
  // many remain.
  int filter(double* configurations, int count, std::size_t joints) const {
    int kept = 0;
    for (int k = 0; k < count; ++k) {
      double* row = configurations + static_cast<std::size_t>(k) * joints;
      std::vector<double> q(row, row + joints);
      if (!fit(q.data())) {
        continue;
      }
      std::copy(q.begin(), q.end(), configurations + static_cast<std::size_t>(kept) * joints);
      ++kept;
    }
    return kept;
  }
};

namespace detail {

// A polynomial over double, kept as its terms.
template <std::size_t V>
struct compiled_polynomial {
  std::vector<std::pair<std::array<std::uint16_t, V>, double>> terms;

  template <class Order>
  static compiled_polynomial from(const polynomial<rational, V, Order>& p) {
    compiled_polynomial out;
    for (const auto& t : p.terms()) {
      out.terms.push_back({t.mon.exponents(), t.coeff.get_d()});
    }
    return out;
  }

  static double power(double x, unsigned e) {
    double r = 1.0;
    for (unsigned k = 0; k < e; ++k) {
      r *= x;
    }
    return r;
  }

  // The value and the sum of the magnitudes of the terms that made it, which
  // is what a denominator is judged against.
  std::pair<double, double> value_and_scale(const double* x) const {
    double value = 0.0;
    double scale = 0.0;
    for (const auto& [e, c] : terms) {
      double m = c;
      for (std::size_t v = 0; v < V; ++v) {
        m *= power(x[v], e[v]);
      }
      value += m;
      scale += std::abs(m);
    }
    return {value, scale};
  }

  double value(const double* x) const { return value_and_scale(x).first; }

  double derivative(const double* x, std::size_t along) const {
    double value = 0.0;
    for (const auto& [e, c] : terms) {
      if (e[along] == 0) {
        continue;
      }
      double m = c * e[along];
      for (std::size_t v = 0; v < V; ++v) {
        m *= power(x[v], v == along ? e[v] - 1u : e[v]);
      }
      value += m;
    }
    return value;
  }
};

template <std::size_t P>
struct compiled_function {
  compiled_polynomial<P> numerator;
  compiled_polynomial<P> denominator;
  bool zero = true;
  bool constant_denominator = true;

  static compiled_function from(const rational_function<P>& f) {
    compiled_function out;
    out.zero = f.is_zero();
    if (!out.zero) {
      out.numerator = compiled_polynomial<P>::from(f.numerator());
      out.denominator = compiled_polynomial<P>::from(f.denominator());
      out.constant_denominator = f.denominator().degree() == 0;
    }
    return out;
  }

  // False when the denominator vanishes at this pose as far as the arithmetic
  // can tell: when cancellation has eaten all but 1e-12 of the terms that made
  // it, the same guard the generated headers carry.
  bool evaluate(const double* pose, double& out) const {
    if (zero) {
      out = 0.0;
      return true;
    }
    const auto [d, scale] = denominator.value_and_scale(pose);
    if (!constant_denominator && !(std::abs(d) > 1e-12 * std::abs(scale))) {
      return false;
    }
    out = numerator.value(pose) / d;
    return true;
  }
};

}  // namespace detail

// A parametric solution over Q(p), evaluated at a pose. Returns half-angle
// variables, as a generated header's solve() does.
template <std::size_t N, std::size_t P>
class runtime_solver {
 public:
  runtime_solver() = default;

  explicit runtime_solver(const codegen::parametric_solution<N, P>& s)
      : dimension_(s.dimension()), one_index_(s.one_index) {
    for (const auto& m : s.action) {
      std::vector<detail::compiled_function<P>> entries;
      for (const auto& e : m.entries) {
        entries.push_back(detail::compiled_function<P>::from(e));
      }
      action_.push_back(std::move(entries));
    }
    for (const auto& row : s.variable_coordinates) {
      std::vector<detail::compiled_function<P>> entries;
      for (const auto& e : row) {
        entries.push_back(detail::compiled_function<P>::from(e));
      }
      coordinates_.push_back(std::move(entries));
    }
    for (const auto& r : s.residual_numerators) {
      residual_numerators_.push_back(detail::compiled_polynomial<N>::from(r));
    }
    if (!s.residual_numerators.empty()) {
      residual_denominator_ = detail::compiled_polynomial<N>::from(s.residual_denominator);
    }
  }

  std::size_t dimension() const noexcept { return dimension_; }

  // Writes up to `capacity` real solutions, N half-angle values each. Returns
  // how many, or -1 if a denominator vanished at this pose.
  int solve(const double* pose, double* out, int capacity) const {
    using matrix = Eigen::MatrixXd;
    const auto d = static_cast<Eigen::Index>(dimension_);
    matrix separating = matrix::Zero(d, d);
    for (std::size_t i = 0; i < N; ++i) {
      const double c = 1.0 + 0.37 * static_cast<double>(i) + 0.11 * static_cast<double>(i * i);
      for (Eigen::Index col = 0; col < d; ++col) {
        for (Eigen::Index row = 0; row < d; ++row) {
          double v = 0.0;
          if (!action_[i][static_cast<std::size_t>(col * d + row)].evaluate(pose, v)) {
            return -1;
          }
          separating(row, col) += c * v;
        }
      }
    }
    std::vector<std::vector<double>> coordinates(N, std::vector<double>(dimension_, 0.0));
    for (std::size_t i = 0; i < N; ++i) {
      for (std::size_t m = 0; m < dimension_; ++m) {
        if (!coordinates_[i][m].evaluate(pose, coordinates[i][m])) {
          return -1;
        }
      }
    }

    Eigen::EigenSolver<matrix> solver(separating.transpose(), true);
    if (solver.info() != Eigen::Success) {
      return -1;
    }
    const auto vectors = solver.eigenvectors();
    constexpr double tolerance = 1e-8;
    int written = 0;
    for (Eigen::Index k = 0; k < vectors.cols() && written < capacity; ++k) {
      const auto v = vectors.col(k);
      const std::complex<double> scale = v(static_cast<Eigen::Index>(one_index_));
      if (std::abs(scale) <= tolerance) {
        continue;
      }
      std::array<double, N> point{};
      bool real = true;
      for (std::size_t i = 0; i < N; ++i) {
        std::complex<double> value{0.0, 0.0};
        for (std::size_t m = 0; m < dimension_; ++m) {
          value += coordinates[i][m] * v(static_cast<Eigen::Index>(m));
        }
        const std::complex<double> x = value / scale;
        if (std::abs(x.imag()) > tolerance) {
          real = false;
        }
        point[i] = x.real();
      }
      if (!real) {
        continue;
      }
      polish(pose, point.data());
      std::copy(point.begin(), point.end(), out + static_cast<std::size_t>(written) * N);
      ++written;
    }
    return written;
  }

 private:
  void residual(const double* t, const double* pose, double* f, double* jacobian) const {
    const double denominator = residual_denominator_.value(t);
    for (std::size_t k = 0; k < N; ++k) {
      f[k] = residual_numerators_[k].value(t) - denominator * pose[k];
      for (std::size_t i = 0; i < N; ++i) {
        jacobian[k * N + i] = residual_numerators_[k].derivative(t, i) -
                              residual_denominator_.derivative(t, i) * pose[k];
      }
    }
  }

  // Up to two Newton steps, each kept only if it lowers the residual.
  void polish(const double* pose, double* t) const {
    if (residual_numerators_.size() != N) {
      return;
    }
    using square = Eigen::Matrix<double, static_cast<int>(N), static_cast<int>(N), Eigen::RowMajor>;
    using vector = Eigen::Matrix<double, static_cast<int>(N), 1>;
    double f[N];
    double jacobian[N * N];
    residual(t, pose, f, jacobian);
    double size = Eigen::Map<const vector>(f).norm();
    for (int step = 0; step < 2 && size > 0.0; ++step) {
      const vector delta =
          Eigen::Map<const square>(jacobian).partialPivLu().solve(-Eigen::Map<const vector>(f));
      double candidate[N];
      for (std::size_t i = 0; i < N; ++i) {
        candidate[i] = t[i] + delta[static_cast<Eigen::Index>(i)];
      }
      double g[N];
      double next_jacobian[N * N];
      residual(candidate, pose, g, next_jacobian);
      const double next = Eigen::Map<const vector>(g).norm();
      if (!(next < size)) {
        break;
      }
      std::copy(candidate, candidate + N, t);
      std::copy(g, g + N, f);
      std::copy(next_jacobian, next_jacobian + N * N, jacobian);
      size = next;
    }
  }

  std::size_t dimension_ = 0;
  std::size_t one_index_ = 0;
  std::vector<std::vector<detail::compiled_function<P>>> action_;
  std::vector<std::vector<detail::compiled_function<P>>> coordinates_;
  std::vector<detail::compiled_polynomial<N>> residual_numerators_;
  detail::compiled_polynomial<N> residual_denominator_;
};

// Three joints placing a point, by whichever route the arm admits. Returns
// joint angles.
class runtime_arm {
 public:
  runtime_arm() = default;

  static runtime_arm decoupled(const decoupled_solution<3>& s) {
    runtime_arm arm;
    arm.decoupled_ = true;
    arm.reduced_ = runtime_solver<2, 2>(s.reduced);
    arm.frame_ = s.frame;
    arm.max_ = 2 * s.reduced.dimension();
    return arm;
  }

  static runtime_arm full(const codegen::parametric_solution<3, 3>& s) {
    runtime_arm arm;
    arm.full_ = runtime_solver<3, 3>(s);
    arm.max_ = s.dimension();
    return arm;
  }

  std::size_t max_configurations() const noexcept { return max_; }
  bool is_decoupled() const noexcept { return decoupled_; }

  int solve(const double* target, double* out, int capacity) const {
    if (!decoupled_) {
      const int found = full_.solve(target, out, capacity);
      for (int k = 0; k < found; ++k) {
        for (int i = 0; i < 3; ++i) {
          out[3 * k + i] = 2.0 * std::atan(out[3 * k + i]);
        }
      }
      return found;
    }
    // The base angle by an arctangent, and each reduced solution twice: facing
    // the target, and half a turn away reaching backwards.
    constexpr double half_turn = 3.14159265358979323846;
    const double a = target[frame_.radial];
    const double b = target[frame_.swept];
    const double radius = std::hypot(a, b);
    const double heading = std::atan2(b, a);
    const double sign = frame_.reversed ? -1.0 : 1.0;
    const double families[2][2] = {{radius, sign * heading},
                                   {-radius, sign * (heading + half_turn)}};
    int written = 0;
    bool described = false;
    std::vector<double> reduced(reduced_.dimension() * 2);
    for (const auto& family : families) {
      const double pose[2] = {family[0], target[frame_.axis]};
      const int found = reduced_.solve(pose, reduced.data(), static_cast<int>(reduced_.dimension()));
      if (found < 0) {
        continue;
      }
      described = true;
      for (int s = 0; s < found && written < capacity; ++s) {
        double* row = out + static_cast<std::size_t>(written) * 3;
        row[0] = std::remainder(family[1], 2.0 * half_turn);
        row[1] = 2.0 * std::atan(reduced[static_cast<std::size_t>(s) * 2]);
        row[2] = 2.0 * std::atan(reduced[static_cast<std::size_t>(s) * 2 + 1]);
        ++written;
      }
    }
    return described ? written : -1;
  }

 private:
  bool decoupled_ = false;
  runtime_solver<2, 2> reduced_;
  runtime_solver<3, 3> full_;
  sweep_frame frame_;
  std::size_t max_ = 0;
};

// Six joints with a spherical wrist, for a full pose. Returns joint angles.
class runtime_wrist_solver {
 public:
  runtime_wrist_solver() = default;

  runtime_wrist_solver(const wrist_decomposition& w, runtime_arm arm)
      : arm_(std::move(arm)) {
    for (std::size_t i = 0; i < 6; ++i) {
      for (std::size_t r = 0; r < 3; ++r) {
        for (std::size_t c = 0; c < 3; ++c) {
          placement_[i][3 * r + c] = w.placement[i](r, c).get_d();
        }
        axis_[i][r] = w.axis[i][r].get_d();
      }
    }
    for (std::size_t r = 0; r < 3; ++r) {
      for (std::size_t c = 0; c < 3; ++c) {
        tool_rotation_[3 * r + c] = w.tool_rotation(r, c).get_d();
      }
      centre_in_tool_[r] = w.centre_in_tool[r].get_d();
    }
  }

  std::size_t max_configurations() const noexcept { return 2 * arm_.max_configurations(); }
  const runtime_arm& arm() const noexcept { return arm_; }

  // Writes up to `capacity` configurations, six angles each, for the tool at
  // `position` with the row-major rotation `rotation`. Returns how many, or -1
  // if the wrist centre falls where the arm solver does not describe.
  int solve(const double* position, const double* rotation, double* out, int capacity) const {
    double centre[3];
    apply(rotation, centre_in_tool_, centre);
    for (int i = 0; i < 3; ++i) {
      centre[i] += position[i];
    }
    std::vector<double> arm_out(arm_.max_configurations() * 3);
    const int arms = arm_.solve(centre, arm_out.data(), static_cast<int>(arm_.max_configurations()));
    if (arms < 0) {
      return -1;
    }

    double tool_t[9], wanted[9];
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        tool_t[3 * i + j] = tool_rotation_[3 * j + i];
      }
    }
    mul(rotation, tool_t, wanted);
    double b5[3], p[3], o6a6[3];
    apply(placement_[4], axis_[4], b5);
    apply(placement_[5], axis_[5], o6a6);
    apply(placement_[4], o6a6, p);
    const double* a4 = axis_[3];

    int written = 0;
    for (int k = 0; k < arms; ++k) {
      double q[6];
      for (int i = 0; i < 3; ++i) {
        q[i] = arm_out[static_cast<std::size_t>(3 * k + i)];
      }
      double r01[9], r12[9], r23[9], r02[9], r03[9], m[9];
      link(0, q[0], r01);
      link(1, q[1], r12);
      link(2, q[2], r23);
      mul(r01, r12, r02);
      mul(r02, r23, r03);
      mul_transposed(r03, wanted, m);
      double ma6[3], t[3];
      apply(m, axis_[5], ma6);
      apply_transposed(placement_[3], ma6, t);

      const double c12 = dot(a4, b5);
      const double denominator = c12 * c12 - 1.0;
      const double alpha = (c12 * dot(b5, p) - dot(a4, t)) / denominator;
      const double beta = (c12 * dot(a4, t) - dot(b5, p)) / denominator;
      double n[3];
      cross(a4, b5, n);
      double gamma2 = (dot(p, p) - alpha * alpha - beta * beta - 2.0 * alpha * beta * c12) / dot(n, n);
      if (gamma2 < -1e-10 * dot(p, p)) {
        continue;  // this arm configuration leaves a rotation the wrist cannot make
      }
      gamma2 = std::max(gamma2, 0.0);
      const double gamma = std::sqrt(gamma2);
      for (int sign = 0; sign < 2 && written < capacity; ++sign) {
        const double g = sign == 0 ? gamma : -gamma;
        double z[3];
        for (int i = 0; i < 3; ++i) {
          z[i] = alpha * a4[i] + beta * b5[i] + g * n[i];
        }
        q[4] = turn_between(b5, p, z);
        q[3] = turn_between(a4, z, t);
        double r34[9], r45[9], r35[9], r35o[9], rest[9];
        link(3, q[3], r34);
        link(4, q[4], r45);
        mul(r34, r45, r35);
        mul(r35, placement_[5], r35o);
        mul_transposed(r35o, m, rest);
        const double* a6 = axis_[5];
        double e[3] = {1.0, 0.0, 0.0};
        if (std::abs(a6[0]) > 0.6) {
          e[0] = 0.0;
          e[1] = 1.0;
        }
        double x[3], nx[3];
        cross(a6, e, x);
        apply(rest, x, nx);
        q[5] = turn_between(a6, x, nx);
        std::copy(q, q + 6, out + static_cast<std::size_t>(written) * 6);
        ++written;
      }
    }
    return written;
  }

 private:
  static double dot(const double* a, const double* b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }
  static void cross(const double* a, const double* b, double* c) {
    c[0] = a[1] * b[2] - a[2] * b[1];
    c[1] = a[2] * b[0] - a[0] * b[2];
    c[2] = a[0] * b[1] - a[1] * b[0];
  }
  static void mul(const double* a, const double* b, double* c) {
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
        c[3 * i + j] = a[3 * i] * b[j] + a[3 * i + 1] * b[3 + j] + a[3 * i + 2] * b[6 + j];
  }
  static void mul_transposed(const double* a, const double* b, double* c) {
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
        c[3 * i + j] = a[i] * b[j] + a[3 + i] * b[3 + j] + a[6 + i] * b[6 + j];
  }
  static void apply(const double* m, const double* v, double* out) {
    for (int i = 0; i < 3; ++i) out[i] = m[3 * i] * v[0] + m[3 * i + 1] * v[1] + m[3 * i + 2] * v[2];
  }
  static void apply_transposed(const double* m, const double* v, double* out) {
    for (int i = 0; i < 3; ++i) out[i] = m[i] * v[0] + m[3 + i] * v[1] + m[6 + i] * v[2];
  }
  static void rotation(const double* w, double angle, double* r) {
    const double c = std::cos(angle), s = std::sin(angle), v = 1.0 - c;
    r[0] = c + w[0] * w[0] * v;        r[1] = w[0] * w[1] * v - w[2] * s; r[2] = w[0] * w[2] * v + w[1] * s;
    r[3] = w[1] * w[0] * v + w[2] * s; r[4] = c + w[1] * w[1] * v;        r[5] = w[1] * w[2] * v - w[0] * s;
    r[6] = w[2] * w[0] * v - w[1] * s; r[7] = w[2] * w[1] * v + w[0] * s; r[8] = c + w[2] * w[2] * v;
  }
  void link(int joint, double angle, double* out) const {
    double turn[9];
    rotation(axis_[joint], angle, turn);
    mul(placement_[joint], turn, out);
  }
  static double turn_between(const double* w, const double* u, const double* v) {
    double up[3], vp[3], c[3];
    const double wu = dot(w, u), wv = dot(w, v);
    for (int i = 0; i < 3; ++i) {
      up[i] = u[i] - w[i] * wu;
      vp[i] = v[i] - w[i] * wv;
    }
    cross(up, vp, c);
    return std::atan2(dot(w, c), dot(up, vp));
  }

  runtime_arm arm_;
  double placement_[6][9] = {};
  double axis_[6][3] = {};
  double tool_rotation_[9] = {};
  double centre_in_tool_[3] = {};
};

// The whole construction, from an exact chain to a solver for the full pose:
// the wrist decomposition, the arm by the decoupling or else by
// reconstruction. Empty when the chain has no spherical wrist or the arm could
// not be solved; `why` then says which.
inline std::optional<runtime_wrist_solver> build_wrist_solver(const chain<rational>& robot,
                                                              std::string* why = nullptr) {
  const auto wrist = decompose_spherical_wrist(robot);
  if (!wrist.ok()) {
    if (why != nullptr) {
      *why = to_string(wrist.status);
    }
    return std::nullopt;
  }
  const auto decoupled = decoupled_position_ik<3>(wrist.arm);
  if (decoupled.ok()) {
    return runtime_wrist_solver(wrist, runtime_arm::decoupled(decoupled));
  }
  const auto reconstructed = reconstructed_position_ik<3, 3>(wrist.arm, {0, 1, 2});
  if (!reconstructed.ok()) {
    if (why != nullptr) {
      *why = std::string("the arm that places the wrist centre was not solved: ") +
             to_string(reconstructed.status);
    }
    return std::nullopt;
  }
  return runtime_wrist_solver(wrist, runtime_arm::full(reconstructed.solution));
}

}  // namespace ik
}  // namespace varietas

#endif
