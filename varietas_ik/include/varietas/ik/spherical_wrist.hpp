#ifndef VARIETAS_IK_SPHERICAL_WRIST_HPP
#define VARIETAS_IK_SPHERICAL_WRIST_HPP

#include <array>
#include <cstddef>
#include <sstream>
#include <string>
#include <vector>

#include "varietas/codegen/emit.hpp"
#include "varietas/codegen/parametric_solution.hpp"
#include "varietas/codegen/rational.hpp"
#include "varietas/core/config.hpp"
#include "varietas/ik/decoupled_ik.hpp"
#include "varietas/ik/emit_decoupled.hpp"
#include "varietas/kinematics/chain.hpp"
#include "varietas/kinematics/rigid_transform.hpp"

// Six joints and a full pose, by splitting position from orientation.
//
// A tool pose is six parameters, and adjoining six parameters is far past what
// a parametric solve can carry. But most industrial arms are built so that it
// never has to: their last three axes meet in a point, the wrist centre, and
// turning those three joints moves the tool about that point without moving
// the point. Where the wrist centre is depends on the first three joints alone.
// So the target pose fixes the wrist centre, finding the first three joints is
// a three-joint position problem, which is what varietas solves, and the last
// three are whatever rotation remains, which comes in closed form. This is
// Pieper's decomposition.
//
// Nothing about it is approximate. Whether the three wrist axes meet is decided
// over Q, on the exact chain, and so is the wrist centre. The arm is solved by
// the decoupling where it applies and by reconstruction where it does not, and
// either way every configuration of the arm is returned, the count certified by
// the quotient dimension. The remaining rotation is split into three joint
// angles by the subproblems of Paden and Kahan, which return every solution:
// two for a wrist away from its singularity. Each arm configuration therefore
// yields two configurations of the whole manipulator, and a general six-joint
// arm with a spherical wrist has at most eight, which is what this returns.
//
// The wrist axes need not be orthogonal to one another. They need only meet.
namespace varietas {
namespace ik {

enum class wrist_status {
  ok,
  // Not six actuated joints, all revolute.
  wrong_joints,
  // Two of the last three axes are parallel, so they cannot fix a point.
  wrist_axes_parallel,
  // The last three axes do not pass through one point: there is no wrist
  // centre, and the position no longer depends on the first three joints alone.
  wrist_axes_do_not_meet,
};

inline const char* to_string(wrist_status status) {
  switch (status) {
    case wrist_status::ok:
      return "ok";
    case wrist_status::wrong_joints:
      return "the decomposition is for six revolute joints";
    case wrist_status::wrist_axes_parallel:
      return "two of the last three axes are parallel, so the wrist has no centre";
    case wrist_status::wrist_axes_do_not_meet:
      return "the last three axes do not meet in a point, so the arm has no spherical wrist";
  }
  return "unknown";
}

// The split of a six-joint arm at its wrist centre, exactly.
struct wrist_decomposition {
  wrist_status status = wrist_status::ok;

  // The first three joints, with the tool replaced by the wrist centre. This is
  // the position problem, and it goes to the three-joint solvers unchanged.
  chain<rational> arm;

  // The wrist centre, in the base frame with every joint at zero, and in the
  // tool frame, where it is constant because it lies on the last axis.
  vector3<rational> centre;
  vector3<rational> centre_in_tool;

  // What the orientation needs at run time: the rotation of every joint's
  // placement and its axis, in the joint's own frame, and the tool's rotation.
  std::array<matrix3<rational>, 6> placement;
  std::array<vector3<rational>, 6> axis;
  matrix3<rational> tool_rotation = matrix3<rational>::identity();

  std::vector<std::string> joint_names;

  // The six joints' ranges, from the robot description.
  detail::emitted_limits limits;

  bool ok() const noexcept { return status == wrist_status::ok; }
};

namespace detail {

inline rational dot(const vector3<rational>& a, const vector3<rational>& b) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

inline vector3<rational> cross(const vector3<rational>& a, const vector3<rational>& b) {
  return vector3<rational>(a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2],
                           a[0] * b[1] - a[1] * b[0]);
}

// Whether the point lies on the line through `through` along `direction`.
inline bool on_line(const vector3<rational>& point, const vector3<rational>& through,
                    const vector3<rational>& direction) {
  return cross(point - through, direction).is_zero();
}

}  // namespace detail

inline wrist_decomposition decompose_spherical_wrist(const chain<rational>& robot) {
  wrist_decomposition result;
  const chain<rational> folded = robot.fold_fixed_joints();
  if (folded.joints().size() != 6 || folded.degrees_of_freedom() != 6) {
    result.status = wrist_status::wrong_joints;
    return result;
  }
  for (const auto& j : folded.joints()) {
    if (j.type != joint_type::revolute) {
      result.status = wrist_status::wrong_joints;
      return result;
    }
  }

  // Every joint's frame with every joint at zero, exactly. A revolute joint's
  // own rotation fixes its axis, so the axis of joint i is the line through
  // frame i's origin along frame i's rotation of the joint's axis, whatever the
  // joint angles before it, up to the motion of the frame itself.
  std::array<rigid_transform<rational>, 6> frame;
  rigid_transform<rational> running = rigid_transform<rational>::identity();
  for (std::size_t i = 0; i < 6; ++i) {
    const auto& j = folded.joints()[i];
    running = running * j.origin;
    frame[i] = running;
    result.placement[i] = j.origin.rotation();
    result.axis[i] = j.axis;
    result.joint_names.push_back(j.name);
  }
  result.tool_rotation = folded.tool().rotation();
  result.limits = detail::limits_of(folded);

  std::array<vector3<rational>, 3> through;
  std::array<vector3<rational>, 3> along;
  for (std::size_t k = 0; k < 3; ++k) {
    through[k] = frame[3 + k].translation();
    along[k] = frame[3 + k].rotation() * folded.joints()[3 + k].axis;
  }

  // Consecutive wrist axes must not be parallel, and that is a fact about the
  // wrist rather than about a configuration: the angle between one axis and the
  // next is fixed by the placement between them. The fourth and sixth axes are
  // another matter. In the commonest wrist, rolling, pitching and rolling
  // again, they are collinear whenever the pitch is zero, which includes the
  // zero configuration these lines are taken in. That is the wrist's
  // singularity, not a defect, and the lines still meet where they should, so
  // the fourth and sixth are allowed to coincide here.
  const vector3<rational> normal = detail::cross(along[0], along[1]);
  const rational norm = detail::dot(normal, normal);
  if (coefficient_traits<rational>::is_zero(norm) ||
      detail::cross(along[1], along[2]).is_zero()) {
    result.status = wrist_status::wrist_axes_parallel;
    return result;
  }

  // The meeting point of the fourth and fifth axes, exactly: the point of the
  // fourth line at parameter s = ((p5 - p4) x d5) . (d4 x d5) / |d4 x d5|^2.
  const rational s = detail::dot(detail::cross(through[1] - through[0], along[1]), normal) / norm;
  const vector3<rational> centre = through[0] + s * along[0];
  if (!detail::on_line(centre, through[1], along[1]) ||
      !detail::on_line(centre, through[2], along[2])) {
    result.status = wrist_status::wrist_axes_do_not_meet;
    return result;
  }
  result.centre = centre;

  // The arm: the first three joints, and a tool at the wrist centre, which
  // lies on the fourth axis and is therefore fixed in the third joint's frame.
  chain<rational> arm(robot.name() + "_arm");
  for (std::size_t i = 0; i < 3; ++i) {
    arm.add_joint(folded.joints()[i]);
  }
  const auto& third = frame[2];
  arm.set_tool(rigid_transform<rational>::translation_only(
      third.rotation().transpose() * (centre - third.translation())));
  result.arm = arm;

  // The wrist centre in the tool frame, where it is fixed because it lies on
  // the sixth axis, which turns with the tool.
  const rigid_transform<rational> tool = frame[5] * folded.tool();
  result.centre_in_tool = tool.rotation().transpose() * (centre - tool.translation());
  return result;
}

namespace detail {

inline std::string matrix_literal(const matrix3<rational>& m) {
  std::string out = "{";
  for (std::size_t i = 0; i < 3; ++i) {
    for (std::size_t j = 0; j < 3; ++j) {
      out += codegen::detail::emit_rational(m(i, j));
      out += (i == 2 && j == 2) ? "" : ", ";
    }
  }
  return out + "}";
}

inline std::string vector_literal(const vector3<rational>& v) {
  return "{" + codegen::detail::emit_rational(v[0]) + ", " + codegen::detail::emit_rational(v[1]) +
         ", " + codegen::detail::emit_rational(v[2]) + "}";
}

}  // namespace detail

// The six-joint solver, as text to follow the arm's solver in one header.
//
// `arm_name` is the struct that solves the arm and `arm_angles` says whether it
// returns angles, as the decoupled wrapper does, or half-angle variables, as a
// solver emitted directly does.
inline std::string wrist_epilogue(const wrist_decomposition& w, const std::string& arm_name,
                                  bool arm_angles, std::size_t arm_max,
                                  const std::string& name) {
  VARIETAS_ASSERT(w.ok());
  std::ostringstream out;
  out << "// The whole manipulator: the arm above, which places the wrist centre, and\n"
      << "// the wrist, whose three axes meet there and which supplies the rotation.\n"
      << "//\n"
      << "// Pieper's decomposition. The target pose fixes the wrist centre; the arm\n"
      << "// solver returns every configuration of the first three joints that puts it\n"
      << "// there; and for each, the rotation the wrist must still supply is split into\n"
      << "// its three joint angles by the subproblems of Paden and Kahan, which give two\n"
      << "// solutions for a wrist away from its singularity. Returns joint angles.\n"
      << "struct " << name << " {\n"
      << "  static constexpr std::size_t num_joints = 6;\n"
      << "  static constexpr std::size_t max_configurations = " << 2 * arm_max << ";\n"
      << "  using status = " << arm_name << "::status;\n\n"
      << "  // Joint placements' rotations, row-major, and axes in each joint's frame.\n"
      << "  static constexpr double placement[6][9] = {\n";
  for (std::size_t i = 0; i < 6; ++i) {
    out << "    " << detail::matrix_literal(w.placement[i]) << ",  // " << w.joint_names[i] << "\n";
  }
  out << "  };\n  static constexpr double axis[6][3] = {\n";
  for (std::size_t i = 0; i < 6; ++i) {
    out << "    " << detail::vector_literal(w.axis[i]) << ",\n";
  }
  out << "  };\n"
      << "  static constexpr double tool_rotation[9] = " << detail::matrix_literal(w.tool_rotation)
      << ";\n"
      << "  // The wrist centre in the tool frame, where it does not move.\n"
      << "  static constexpr double centre_in_tool[3] = " << detail::vector_literal(w.centre_in_tool)
      << ";\n";
  out << R"CODE(
  static double dot(const double* a, const double* b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }
  static void cross(const double* a, const double* b, double* c) {
    c[0] = a[1] * b[2] - a[2] * b[1];
    c[1] = a[2] * b[0] - a[0] * b[2];
    c[2] = a[0] * b[1] - a[1] * b[0];
  }
  // c = a b, and c = a^T b, for row-major 3x3 matrices.
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
  // The rotation by `angle` about the unit vector `w`, by Rodrigues' formula.
  static void rotation(const double* w, double angle, double* r) {
    const double c = std::cos(angle), s = std::sin(angle), v = 1.0 - c;
    r[0] = c + w[0] * w[0] * v;        r[1] = w[0] * w[1] * v - w[2] * s; r[2] = w[0] * w[2] * v + w[1] * s;
    r[3] = w[1] * w[0] * v + w[2] * s; r[4] = c + w[1] * w[1] * v;        r[5] = w[1] * w[2] * v - w[0] * s;
    r[6] = w[2] * w[0] * v - w[1] * s; r[7] = w[2] * w[1] * v + w[0] * s; r[8] = c + w[2] * w[2] * v;
  }
  // The joint's placement followed by its turn: one link of the chain.
  static void link(int joint, double angle, double* out) {
    double turn[9];
    rotation(axis[joint], angle, turn);
    mul(placement[joint], turn, out);
  }
  // Paden-Kahan subproblem 1: the angle about w that carries u to v.
  static double turn_between(const double* w, const double* u, const double* v) {
    double up[3], vp[3], c[3];
    const double wu = dot(w, u), wv = dot(w, v);
    for (int i = 0; i < 3; ++i) { up[i] = u[i] - w[i] * wu; vp[i] = v[i] - w[i] * wv; }
    cross(up, vp, c);
    return std::atan2(dot(w, c), dot(up, vp));
  }

  // Writes up to `capacity` configurations into `out`, six joint angles each,
  // for the tool at `position` with orientation `rotation`, a row-major 3x3
  // rotation matrix in the base frame. Returns how many were written, or -1 if
  // the wrist centre falls on the locus the arm solver does not describe.
  static int solve(const double* position, const double* rotation_matrix, double* out,
                   int capacity, status* state = nullptr) {
    double centre[3];
    apply(rotation_matrix, centre_in_tool, centre);
    for (int i = 0; i < 3; ++i) centre[i] += position[i];

)CODE";
  out << "    constexpr int arm_max = " << arm_max << ";\n"
      << "    double arm_out[arm_max * 3];\n"
      << "    status arm_state{};\n"
      << "    const int arms = " << arm_name << "::solve(centre, arm_out, arm_max, &arm_state);\n"
      << "    const bool arm_angles = " << (arm_angles ? "true" : "false") << ";\n";
  out << R"CODE(    if (arms < 0) {
      if (state != nullptr) { *state = arm_state; }
      return -1;
    }

    // The rotation the wrist must supply, M = R03^T R Rtool^T, is written
    // O4 R4 O5 R5 O6 R6 with Ok the placements and Rk the turns. The sixth turn
    // fixes the sixth axis, so M a6 = O4 R4 O5 R5 O6 a6, which is subproblem 2
    // for the fourth and fifth angles: Rot(a4, q4) Rot(b5, q5) p = t with
    // b5 = O5 a5, p = O5 O6 a6 and t = O4^T M a6. The sixth follows by
    // subproblem 1.
    double tool_t[9], wanted[9];
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j) tool_t[3 * i + j] = tool_rotation[3 * j + i];
    mul(rotation_matrix, tool_t, wanted);

    double b5[3], p[3], o6a6[3];
    apply(placement[4], axis[4], b5);
    apply(placement[5], axis[5], o6a6);
    apply(placement[4], o6a6, p);
    const double* a4 = axis[3];

    int written = 0;
    for (int k = 0; k < arms; ++k) {
      double q[6];
      for (int i = 0; i < 3; ++i) {
        const double v = arm_out[3 * k + i];
        q[i] = arm_angles ? v : 2.0 * std::atan(v);
      }
      double r01[9], r12[9], r23[9], r02[9], r03[9], m[9];
      link(0, q[0], r01);
      link(1, q[1], r12);
      link(2, q[2], r23);
      mul(r01, r12, r02);
      mul(r02, r23, r03);
      mul_transposed(r03, wanted, m);

      double ma6[3], t[3];
      apply(m, axis[5], ma6);
      apply_transposed(placement[3], ma6, t);

      const double c12 = dot(a4, b5);
      const double denominator = c12 * c12 - 1.0;
      const double alpha = (c12 * dot(b5, p) - dot(a4, t)) / denominator;
      const double beta = (c12 * dot(a4, t) - dot(b5, p)) / denominator;
      double n[3];
      cross(a4, b5, n);
      double gamma2 = (dot(p, p) - alpha * alpha - beta * beta - 2.0 * alpha * beta * c12) / dot(n, n);
      // Negative means this wrist cannot supply the rotation this arm
      // configuration leaves it. A wrist whose consecutive axes are not
      // perpendicular reaches only some orientations, so for such a wrist an
      // arm configuration can put the centre in place and still leave a
      // rotation no setting of the last three joints makes, and that
      // configuration of the arm has no completion. Only a value within
      // rounding of zero is the double root, where the two wrist solutions
      // meet, and is taken as zero.
      if (gamma2 < -1e-10 * dot(p, p)) { continue; }
      if (gamma2 < 0.0) { gamma2 = 0.0; }
      const double gamma = std::sqrt(gamma2);

      for (int sign = 0; sign < 2; ++sign) {
        if (written >= capacity) { break; }
        const double g = sign == 0 ? gamma : -gamma;
        double z[3];
        for (int i = 0; i < 3; ++i) z[i] = alpha * a4[i] + beta * b5[i] + g * n[i];
        q[4] = turn_between(b5, p, z);
        q[3] = turn_between(a4, z, t);

        // What is left for the sixth joint: N = (O4 R4 O5 R5 O6)^T M, a turn
        // about a6, read off any vector perpendicular to a6.
        double r34[9], r45[9], r35[9], r35o[9], rest[9];
        link(3, q[3], r34);
        link(4, q[4], r45);
        mul(r34, r45, r35);
        mul(r35, placement[5], r35o);
        mul_transposed(r35o, m, rest);
        const double* a6 = axis[5];
        double e[3] = {1.0, 0.0, 0.0};
        if (std::abs(a6[0]) > 0.6) { e[0] = 0.0; e[1] = 1.0; }
        double x[3], nx[3];
        cross(a6, e, x);
        apply(rest, x, nx);
        q[5] = turn_between(a6, x, nx);

        double* row = out + static_cast<std::size_t>(written) * num_joints;
        for (int i = 0; i < 6; ++i) row[i] = q[i];
        ++written;
      }
    }
    if (state != nullptr) { *state = arm_state; }
    return written;
  }
)CODE";
  out << detail::limits_members(w.limits, "const double* position, const double* rotation_matrix",
                                "solve(position, rotation_matrix");
  out << "};\n";
  return out.str();
}

// A header for the whole manipulator, its arm solved by the decoupling.
inline std::string emit_spherical_wrist(const wrist_decomposition& w,
                                        const decoupled_solution<3>& arm,
                                        codegen::emit_options options) {
  VARIETAS_ASSERT(w.ok() && arm.ok());
  VARIETAS_ASSERT(options.runtime == codegen::runtime_kind::eigen);
  const std::string name = options.name;
  const std::string arm_name = name + "_arm";
  const std::string reduced_name = name + "_arm_reduced";
  if (options.guard.empty()) {
    options.guard = codegen::detail::default_guard(options);
  }
  options.epilogue = decoupled_epilogue(arm, reduced_name, arm_name) + "\n" +
                     wrist_epilogue(w, arm_name, true, arm.branches, name);
  options.name = reduced_name;
  return codegen::emit(arm.reduced, options);
}

// A header for the whole manipulator, its arm solved in full over Q(x, y, z).
inline std::string emit_spherical_wrist(const wrist_decomposition& w,
                                        const codegen::parametric_solution<3, 3>& arm,
                                        codegen::emit_options options) {
  VARIETAS_ASSERT(w.ok());
  VARIETAS_ASSERT(options.runtime == codegen::runtime_kind::eigen);
  const std::string name = options.name;
  const std::string arm_name = name + "_arm";
  if (options.guard.empty()) {
    options.guard = codegen::detail::default_guard(options);
  }
  options.epilogue = wrist_epilogue(w, arm_name, false, arm.dimension(), name);
  options.name = arm_name;
  return codegen::emit(arm, options);
}

}  // namespace ik
}  // namespace varietas

#endif
