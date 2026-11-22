#ifndef VARIETAS_IK_LIFTED_IDEAL_HPP
#define VARIETAS_IK_LIFTED_IDEAL_HPP

#include <cstddef>
#include <vector>

#include "varietas/codegen/modular_groebner.hpp"
#include "varietas/codegen/rational.hpp"
#include "varietas/core/polynomial.hpp"
#include "varietas/kinematics/chain.hpp"
#include "varietas/kinematics/rationalize.hpp"
#include "varietas/kinematics/rigid_transform.hpp"

// The inverse kinematics ideal at a fixed target over Q, lifted from its
// images modulo primes (varietas/codegen/modular_groebner.hpp) rather than
// computed by Buchberger over Q.
//
// These are the functions of varietas/kinematics/rationalize.hpp of the same
// names without the prefix, with the same residuals, the same saturation and
// the same reduced basis as the answer, and with one difference in what is
// known about it. Buchberger's basis is exact. The lifted one has been checked
// over Q to be a Gröbner basis of an ideal containing the true one, which is
// exact, and to be no larger than it with high probability, which is not; see
// the header that lifts it for what that rests on. What it costs is the
// difference between the two paths: a six-joint arm at a full pose, which
// Buchberger over Q did not finish in eight minutes, takes seconds.
//
// Where the check cannot be made, a saturation that is not zero-dimensional
// such as a redundant arm's, the basis is computed by Buchberger, as the
// statistics report.
namespace varietas {
namespace ik {

template <std::size_t N, class Order>
std::vector<polynomial<rational, N, Order>> lifted_kinematic_ideal_generators(
    const chain<rational>& robot, const std::vector<polynomial<rational, N, Order>>& residuals,
    modular::statistics* statistics = nullptr, const modular::options& options = {}) {
  const auto multipliers = half_angle_multipliers<N, Order>(robot);
  if (multipliers.empty()) {
    return modular::groebner_basis(residuals, statistics, options);
  }
  return modular::saturate_by_product(residuals, multipliers, statistics, options);
}

template <std::size_t N, class Order>
std::vector<polynomial<rational, N, Order>> lifted_pose_ideal_generators(
    const chain<rational>& robot, const rigid_transform<rational>& target,
    modular::statistics* statistics = nullptr, const modular::options& options = {}) {
  return lifted_kinematic_ideal_generators<N, Order>(
      robot, split_pose_residuals<N, Order>(robot, target), statistics, options);
}

template <std::size_t N, class Order>
std::vector<polynomial<rational, N, Order>> lifted_position_ideal_generators(
    const chain<rational>& robot, const vector3<rational>& target,
    modular::statistics* statistics = nullptr, const modular::options& options = {}) {
  const auto map = rational_forward_kinematics<N, Order>(robot);
  return lifted_kinematic_ideal_generators<N, Order>(robot, position_residuals(map, target),
                                                     statistics, options);
}

}  // namespace ik
}  // namespace varietas

#endif
