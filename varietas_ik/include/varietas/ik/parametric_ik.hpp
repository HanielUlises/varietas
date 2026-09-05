#ifndef VARIETAS_IK_PARAMETRIC_IK_HPP
#define VARIETAS_IK_PARAMETRIC_IK_HPP

#include <array>
#include <cstddef>
#include <string>
#include <vector>

#include "varietas/codegen/parametric_solution.hpp"
#include "varietas/codegen/rational.hpp"
#include "varietas/codegen/rational_function.hpp"
#include "varietas/core/config.hpp"
#include "varietas/core/monomial.hpp"
#include "varietas/core/order/grevlex.hpp"
#include "varietas/core/order/order_id.hpp"
#include "varietas/core/quotient/quotient_basis.hpp"
#include "varietas/ik/cast_chain.hpp"
#include "varietas/kinematics/chain.hpp"
#include "varietas/kinematics/rationalize.hpp"

// The join between the robot and the emitter.
//
// varietas_kinematics turns a chain and a *given* pose into an ideal in the
// joint variables; varietas_codegen turns a system solved over Q(p) into a
// header. Neither knows about the other, and the step between them is the one
// this package supplies: pose the inverse kinematics with the target adjoined
// to the coefficient field instead of substituted into it, so that Buchberger
// runs once for the arm rather than once for every point it is asked to reach.
//
// That single change is what makes emission worth doing at all. A basis
// computed at a pose answers that pose and is then thrown away, which is fine
// offline and useless online. A basis computed over Q(p) is an object whose
// entries are functions of the pose, and evaluating those functions is
// arithmetic, which is exactly what emit() writes out.
namespace varietas {
namespace ik {

// Why a chain and a choice of coordinates did not produce a solution.
//
// These are reported rather than asserted because every one of them is a
// property of the robot handed in, not a bug: an arm can genuinely fail to be
// zero-dimensional at a general pose, and a caller can genuinely ask for
// coordinates the arm does not control. Refusing with a reason is the useful
// answer.
enum class parametric_ik_status {
  ok,

  // The chain does not have the N actuated joints the call was instantiated for.
  wrong_degrees_of_freedom,

  // The requested pose coordinates were not P distinct indices below three.
  bad_coordinates,

  // Fewer pose coordinates than unknowns, which cannot give a finite solution
  // set and is caught before anything expensive runs.
  //
  // The residuals are P polynomials in N variables, so by Krull's height
  // theorem every component of their zero set has dimension at least N - P.
  // Saturation only removes components, so it cannot raise the dimension of the
  // ones that survive: the saturated ideal is either the unit ideal, meaning the
  // arm reaches no such point at all, or it is still positive-dimensional. There
  // is nothing to emit either way, and saying so here costs nothing, whereas
  // finding out by running Buchberger over Q(p) on an underdetermined system is
  // the most expensive way to learn it.
  underdetermined,

  // More pose coordinates than unknowns, which cannot have a solution at a
  // general pose and is likewise caught before anything expensive runs.
  //
  // The tool positions an N-joint arm can reach form a variety of dimension at
  // most N, and the parameters are transcendentals, a general point of
  // P-space, not a point of that image. With P > N the image is a proper
  // subvariety of the space the pose ranges over, a general pose does not lie
  // on it, and the ideal is the unit ideal: the arm simply cannot reach the
  // point it is being asked about.
  //
  // Taken with the case above this leaves P = N as the only arrangement that
  // can produce a parametric solver, which is worth knowing before choosing
  // --coords rather than after waiting for a Grobner basis to say so.
  overdetermined,

  // A coordinate that was *not* requested is nevertheless moved by the arm.
  //
  // Leaving a coordinate out of the parameter list drops its equation, and
  // dropping an equation is only legitimate when the equation says nothing,
  // when the tool's coordinate there is identically zero as a polynomial in the
  // joint variables, as the z of a planar arm is. If it is not identically
  // zero, the dropped equation constrained the joints, and the system that
  // remains is a different and larger variety than the one asked about. So the
  // condition is checked rather than assumed.
  dropped_coordinate_is_not_identically_zero,

  // The saturated ideal over Q(p) is not zero-dimensional, so there is no
  // finite solution set at a general pose and no action matrix to emit. An arm
  // with more joints than the coordinates constrain has a positive-dimensional
  // fibre, and that is the usual reason to land here.
  not_zero_dimensional,
};

inline const char* to_string(parametric_ik_status status) {
  switch (status) {
    case parametric_ik_status::ok:
      return "ok";
    case parametric_ik_status::wrong_degrees_of_freedom:
      return "chain does not have the expected number of actuated joints";
    case parametric_ik_status::bad_coordinates:
      return "pose coordinates must be distinct indices below three";
    case parametric_ik_status::underdetermined:
      return "fewer pose coordinates than unknowns, so the solution set cannot be finite";
    case parametric_ik_status::overdetermined:
      return "more pose coordinates than unknowns, so a general pose is out of reach";
    case parametric_ik_status::dropped_coordinate_is_not_identically_zero:
      return "an unrequested coordinate is moved by the arm, so its equation cannot be dropped";
    case parametric_ik_status::not_zero_dimensional:
      return "the ideal over Q(p) is not zero-dimensional at a general pose";
  }
  return "unknown";
}

template <std::size_t N, std::size_t P>
struct parametric_ik_result {
  parametric_ik_status status = parametric_ik_status::ok;

  // Meaningful only when status is ok.
  codegen::parametric_solution<N, P> solution;

  // dim_k A: the number of inverse kinematics branches at a general pose,
  // counted with multiplicity and over the complex numbers. Filled in whenever
  // the quotient could be computed at all, since it is informative even when it
  // is zero or the system turned out not to be zero-dimensional.
  std::size_t branches = 0;

  // For dropped_coordinate_is_not_identically_zero and bad_coordinates, the
  // coordinate at fault.
  std::size_t offending_coordinate = 0;

  bool ok() const noexcept { return status == parametric_ik_status::ok; }
  explicit operator bool() const noexcept { return ok(); }
};

namespace detail {

// The name given to the ring variable a joint contributes. The variable is not
// the joint angle, being tan(q/2) for a revolute joint, and the generated
// header's comments should not pretend otherwise.
template <class Coeff>
inline std::string variable_name(const joint<Coeff>& j) {
  switch (j.type) {
    case joint_type::revolute:
      return "t_" + j.name;
    case joint_type::prismatic:
      return "d_" + j.name;
    case joint_type::fixed:
      break;
  }
  return j.name;
}

inline std::string coordinate_name(std::size_t coordinate) {
  static const char* const names[3] = {"x", "y", "z"};
  return coordinate < 3 ? names[coordinate] : "p" + std::to_string(coordinate);
}

}  // namespace detail

// The inverse kinematics of a tool position, solved once over Q(p).
//
// `coordinates` names which of the three position coordinates are adjoined to
// the field, in the order the generated header will expect them: passing
// {0, 1} to a planar arm gives a solver whose pose argument is (x, y), and the
// z equation is dropped after being checked to be vacuous. Passing all three
// to a spatial arm is the ordinary case.
//
// The result is exactly what emit() consumes, which is the point: from here the
// path to a header is one function call, and nothing between a URDF and that
// header involves solving anything at a particular pose.
template <std::size_t N, std::size_t P>
parametric_ik_result<N, P> parametric_position_ik(
    const chain<rational>& robot, const std::array<std::size_t, P>& coordinates,
    buchberger_statistics* statistics = nullptr) {
  static_assert(P >= 1 && P <= 3, "a tool position has three coordinates");

  static_assert(N >= 1, "a chain with no actuated joint has nothing to solve for");

  using field = rational_function<P>;
  using poly = polynomial<field, N, grevlex>;

  parametric_ik_result<N, P> result;

  // Checked first, and by counting rather than by computing. Only P = N can
  // give a finite, nonempty solution set at a general pose; see the comments on
  // the two statuses for why each of the other cases settles itself.
  if (P < N) {
    result.status = parametric_ik_status::underdetermined;
    return result;
  }
  if (P > N) {
    result.status = parametric_ik_status::overdetermined;
    return result;
  }

  std::array<bool, 3> constrained{};
  constrained.fill(false);
  for (std::size_t k = 0; k < P; ++k) {
    const std::size_t coordinate = coordinates[k];
    if (coordinate >= 3 || constrained[coordinate]) {
      result.status = parametric_ik_status::bad_coordinates;
      result.offending_coordinate = coordinate;
      return result;
    }
    constrained[coordinate] = true;
  }

  const chain<field> parametric = cast_chain<field>(robot).fold_fixed_joints();
  if (parametric.degrees_of_freedom() != N) {
    result.status = parametric_ik_status::wrong_degrees_of_freedom;
    return result;
  }

  const auto map = rational_forward_kinematics<N, grevlex>(parametric);
  const poly denominator = map.denominator();

  for (std::size_t i = 0; i < 3; ++i) {
    if (!constrained[i] && !map.translation(i).is_zero()) {
      result.status = parametric_ik_status::dropped_coordinate_is_not_identically_zero;
      result.offending_coordinate = i;
      return result;
    }
  }

  // The residuals, with the target a parameter of the field rather than a
  // constant of the ring. This is the only line that differs from the ordinary
  // position_residuals, and it is the whole idea.
  std::vector<poly> residuals;
  residuals.reserve(P);
  for (std::size_t k = 0; k < P; ++k) {
    residuals.push_back(map.translation(coordinates[k]) -
                        denominator * poly::constant(field::parameter(k)));
  }

  // Saturating by the half-angle denominators is not optional here for the same
  // reason it is not optional at a fixed pose: clearing them attached the loci
  // t_i = ±i to the variety, and dim_k A counts those too. The count is the
  // certificate that the emitted header has a branch for every configuration
  // and none for anything else, so it has to count configurations only.
  const auto basis = kinematic_ideal_generators<N, grevlex>(parametric, residuals, statistics);
  const auto quotient = standard_monomials(basis);
  result.branches = quotient.dimension();

  if (!quotient.is_zero_dimensional || quotient.dimension() == 0) {
    result.status = parametric_ik_status::not_zero_dimensional;
    return result;
  }

  codegen::parametric_solution<N, P>& solution = result.solution;
  solution.order = order_id::grevlex;
  solution.quotient = quotient;
  solution.one_index = quotient.index_of(monomial<N>::one());

  for (const joint<field>& j : parametric.joints()) {
    if (j.is_actuated()) {
      solution.unknown_names.push_back(detail::variable_name(j));
    }
  }
  for (std::size_t k = 0; k < P; ++k) {
    solution.parameter_names.push_back(detail::coordinate_name(coordinates[k]));
  }

  solution.action.reserve(N);
  for (std::size_t variable = 0; variable < N; ++variable) {
    solution.action.push_back(
        codegen::parametric_action_matrix<N, P, grevlex>(variable, basis, quotient));
  }
  solution.variable_coordinates =
      codegen::parametric_variable_coordinates<N, P, grevlex>(basis, quotient);

  // is_well_formed is what emit() asserts on entry. Checking it here as well
  // means a malformed solution is caught where it was built, with the chain
  // still in hand, rather than inside the emitter.
  VARIETAS_ASSERT(solution.is_well_formed());
  return result;
}

}  // namespace ik
}  // namespace varietas

#endif
