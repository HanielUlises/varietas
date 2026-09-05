#ifndef VARIETAS_IK_DECOUPLED_IK_HPP
#define VARIETAS_IK_DECOUPLED_IK_HPP

#include <array>
#include <cstddef>
#include <string>

#include "varietas/codegen/parametric_solution.hpp"
#include "varietas/codegen/rational.hpp"
#include "varietas/ik/parametric_ik.hpp"
#include "varietas/kinematics/chain.hpp"
#include "varietas/kinematics/rigid_transform.hpp"

// Solving an arm by not adjoining its first joint.
//
// The parametric pipeline costs what it costs because of the number of
// parameters adjoined, not because of the arm: the planar two-link arm solves
// over Q(x, y) in about thirty milliseconds, while the anthropomorphic three-
// link arm over Q(x, y, z) produced nothing in fifteen minutes. Since a
// parametric solve needs one coordinate per joint, three joints means three
// parameters, and three parameters is past what the cancellation in
// rational_function can carry.
//
// A base joint that yaws about a fixed axis does not need to be adjoined. If
// the rest of the arm holds the tool in a plane containing that axis, then
// turning the base sweeps that plane around it, and the tool's position is
// determined by where it sits in the plane, a radius and a height, together
// with the angle the plane has been turned through. The angle is recovered by
// an arctangent, not by an eigenvalue, and what is left is a two-joint problem
// in two parameters: the size the pipeline is comfortable with.
//
// This is the decomposition every closed-form treatment of such an arm begins
// with, and the reason it belongs here rather than in the emitter is that
// whether an arm admits it is a question about the arm. The condition is
// checked rather than assumed, and an arm that does not meet it is refused.
//
// What the caller gets back is a parametric_solution for the reduced problem,
// which emit consumes exactly as it consumes any other, plus the small amount
// of trigonometry needed to put the first joint back.
namespace varietas {
namespace ik {

enum class decoupling_status {
  ok,

  // Fewer than two actuated joints, or a joint count that disagrees with the
  // one the call was instantiated for.
  wrong_degrees_of_freedom,

  // The first joint has to turn, since the whole construction is about the
  // plane it sweeps.
  first_joint_not_revolute,

  // The first axis has to be a coordinate direction, positive or negative. An
  // arbitrary axis can be decoupled too, by rotating the frame first, but that
  // rotation is not in general rational and would leave the exact field.
  first_axis_not_a_coordinate_direction,

  // The first joint's placement has to commute with its own rotation, which for
  // a pure translation means lying along the axis. Anything else moves the axis
  // itself and there is no fixed plane to sweep.
  first_joint_is_displaced,

  // The rest of the arm does not hold the tool in a plane through the axis, so
  // there is nothing to sweep and the problem does not reduce. The status
  // carries the reason from the reduced solve, which is where it was found.
  does_not_reduce,

  // The reduced two-joint problem was itself refused.
  reduced_problem_refused,
};

inline const char* to_string(decoupling_status status) {
  switch (status) {
    case decoupling_status::ok:
      return "ok";
    case decoupling_status::wrong_degrees_of_freedom:
      return "chain does not have the expected number of actuated joints";
    case decoupling_status::first_joint_not_revolute:
      return "the first joint must be revolute";
    case decoupling_status::first_axis_not_a_coordinate_direction:
      return "the first joint's axis must be a coordinate direction";
    case decoupling_status::first_joint_is_displaced:
      return "the first joint's placement must be a translation along its own axis";
    case decoupling_status::does_not_reduce:
      return "the rest of the arm does not hold the tool in a plane through the first axis";
    case decoupling_status::reduced_problem_refused:
      return "the reduced problem was refused";
  }
  return "unknown";
}

// Which coordinate plays which part, once the first axis is known.
//
// Turning about axis k moves the two coordinates a and b and fixes k. Writing
// the tool position as (r cos, r sin) in the (a, b) pair and leaving k alone
// says that the reduced problem lives in the (a, k) plane, a radius measured
// along a and a height measured along k, while b is the coordinate the sweep
// generates and must therefore vanish before the sweep is applied.
struct sweep_frame {
  std::size_t axis = 2;    // k: the direction the first joint turns about
  std::size_t radial = 0;  // a: the radius of the reduced problem
  std::size_t swept = 1;   // b: identically zero on the reduced arm
  bool reversed = false;   // the axis pointed the other way

  // The two coordinates the reduced problem is posed in, radius first.
  std::array<std::size_t, 2> plane() const { return {radial, axis}; }
};

template <std::size_t N>
struct decoupled_solution {
  decoupling_status status = decoupling_status::ok;

  // Set when status is reduced_problem_refused or does_not_reduce.
  parametric_ik_status reduced_status = parametric_ik_status::ok;

  sweep_frame frame;

  // The reduced problem: one joint fewer, two parameters, and emitted exactly
  // as any other parametric solution is.
  codegen::parametric_solution<N - 1, 2> reduced;

  // The name of the joint the sweep put back, for the generated comments.
  std::string first_joint_name;

  // Branches of the whole arm: each solution of the reduced problem occurs at
  // the swept angle and again half a turn away, with the radius negated.
  std::size_t branches = 0;

  bool ok() const noexcept { return status == decoupling_status::ok; }
  explicit operator bool() const noexcept { return ok(); }
};

namespace detail {

// The axis as a signed coordinate direction, if it is one.
inline bool coordinate_direction(const vector3<rational>& axis, std::size_t& index,
                                 bool& reversed) {
  for (std::size_t k = 0; k < 3; ++k) {
    const std::size_t p = (k + 1) % 3;
    const std::size_t q = (k + 2) % 3;
    if (axis[p] != 0 || axis[q] != 0 || axis[k] == 0) {
      continue;
    }
    index = k;
    reversed = axis[k] < 0;
    return true;
  }
  return false;
}

// A translation along the given axis, and no rotation, is the only placement
// that commutes with a rotation about that axis.
inline bool commutes_with_rotation_about(const rigid_transform<rational>& t, std::size_t axis) {
  for (std::size_t i = 0; i < 3; ++i) {
    for (std::size_t j = 0; j < 3; ++j) {
      const rational expected(i == j ? 1 : 0);
      if (t.rotation()(i, j) != expected) {
        return false;
      }
    }
  }
  for (std::size_t i = 0; i < 3; ++i) {
    if (i != axis && t.translation()[i] != 0) {
      return false;
    }
  }
  return true;
}

}  // namespace detail

// The arm with its first joint removed, in that joint's own frame.
//
// The first joint's placement is a translation along its axis, so it commutes
// with the rotation and can be pushed past it into the joint that follows. What
// remains is the subchain as seen from the turning frame, which is the arm the
// reduced problem is about.
inline chain<rational> without_first_joint(const chain<rational>& folded) {
  chain<rational> reduced(folded.name() + "_reduced");
  const auto& joints = folded.joints();
  for (std::size_t i = 1; i < joints.size(); ++i) {
    joint<rational> moved = joints[i];
    if (i == 1) {
      moved.origin = joints[0].origin * moved.origin;
    }
    reduced.add_joint(std::move(moved));
  }
  reduced.set_tool(folded.tool());
  return reduced;
}

// Solve an arm by sweeping its first joint and adjoining only the rest.
template <std::size_t N>
decoupled_solution<N> decoupled_position_ik(const chain<rational>& robot,
                                            buchberger_statistics* statistics = nullptr) {
  static_assert(N >= 2, "there is nothing to decouple in a one-joint arm");

  decoupled_solution<N> result;

  const chain<rational> folded = robot.fold_fixed_joints();
  if (folded.degrees_of_freedom() != N) {
    result.status = decoupling_status::wrong_degrees_of_freedom;
    return result;
  }

  const joint<rational>& first = folded.joints().front();
  if (first.type != joint_type::revolute) {
    result.status = decoupling_status::first_joint_not_revolute;
    return result;
  }

  std::size_t axis = 0;
  bool reversed = false;
  if (!detail::coordinate_direction(first.axis, axis, reversed)) {
    result.status = decoupling_status::first_axis_not_a_coordinate_direction;
    return result;
  }
  if (!detail::commutes_with_rotation_about(first.origin, axis)) {
    result.status = decoupling_status::first_joint_is_displaced;
    return result;
  }

  result.frame.axis = axis;
  result.frame.radial = (axis + 1) % 3;
  result.frame.swept = (axis + 2) % 3;
  result.frame.reversed = reversed;
  result.first_joint_name = first.name;

  const chain<rational> reduced_chain = without_first_joint(folded);

  // The condition the whole construction rests on, checked here rather than
  // inferred from whatever the reduced solve happens to complain about first.
  // With the first joint held still the tool must not move in the swept
  // coordinate at all; if it does, turning the base does not sweep a fixed
  // plane and the position is not determined by a radius and a height.
  //
  // It is a question about the arm, so it is asked of the arm, with no Gröbner
  // basis, just the forward map, and a polynomial that is identically zero or
  // is not.
  {
    const auto map = rational_forward_kinematics<N - 1, grevlex>(reduced_chain);
    if (!map.translation(result.frame.swept).is_zero()) {
      result.status = decoupling_status::does_not_reduce;
      return result;
    }
  }

  const auto reduced =
      parametric_position_ik<N - 1, 2>(reduced_chain, result.frame.plane(), statistics);

  if (!reduced.ok()) {
    result.reduced_status = reduced.status;
    result.status = decoupling_status::reduced_problem_refused;
    return result;
  }

  result.reduced = reduced.solution;
  // Each solution of the reduced problem is reached at the swept angle and
  // again half a turn away, where the radius is measured the other way.
  result.branches = 2 * reduced.branches;
  return result;
}

}  // namespace ik
}  // namespace varietas

#endif
