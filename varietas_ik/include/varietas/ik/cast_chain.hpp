#ifndef VARIETAS_IK_CAST_CHAIN_HPP
#define VARIETAS_IK_CAST_CHAIN_HPP

#include <cstddef>

#include "varietas/kinematics/chain.hpp"
#include "varietas/kinematics/rigid_transform.hpp"

namespace varietas {
namespace ik {

// A chain carried into a larger coefficient field, exactly.
//
// The chain a URDF produces is exact but its coefficients are rationals, and
// the parametric pipeline needs the same chain over Q(p) so that the pose can
// be adjoined to the field rather than to the ring. Nothing about the chain
// changes in the process: the placements are the same numbers, sitting in a
// larger field as constants.
//
// This is deliberately not chain_cast from varietas_kinematics/evaluate.hpp,
// which converts through to_double and back and is documented as a rounding.
// A rounding is the right thing when the destination is the numerical reference
// the algebra is checked against, and the wrong thing here twice over: it would
// throw away the exactness the offline half exists to preserve, and Q(p) has no
// to_double to route through in any case, since a rational function has no
// value until a pose is supplied. What is wanted is the field embedding
// Q -> Q(p), which loses nothing because it is injective.
template <class To, class From>
vector3<To> cast_vector3(const vector3<From>& v) {
  return vector3<To>(To(v[0]), To(v[1]), To(v[2]));
}

template <class To, class From>
matrix3<To> cast_matrix3(const matrix3<From>& m) {
  matrix3<To> result;
  for (std::size_t i = 0; i < 3; ++i) {
    for (std::size_t j = 0; j < 3; ++j) {
      result(i, j) = To(m(i, j));
    }
  }
  return result;
}

template <class To, class From>
rigid_transform<To> cast_rigid_transform(const rigid_transform<From>& t) {
  return rigid_transform<To>(cast_matrix3<To>(t.rotation()),
                             cast_vector3<To>(t.translation()));
}

template <class To, class From>
joint<To> cast_joint(const joint<From>& j) {
  joint<To> result;
  result.name = j.name;
  result.type = j.type;
  result.axis = cast_vector3<To>(j.axis);
  result.origin = cast_rigid_transform<To>(j.origin);
  // Limits are inequalities and are carried as doubles in either field, so they
  // survive the cast untouched — there is nothing to convert.
  result.has_limits = j.has_limits;
  result.lower = j.lower;
  result.upper = j.upper;
  return result;
}

template <class To, class From>
chain<To> cast_chain(const chain<From>& robot) {
  chain<To> result(robot.name());
  for (const joint<From>& j : robot.joints()) {
    result.add_joint(cast_joint<To>(j));
  }
  result.set_tool(cast_rigid_transform<To>(robot.tool()));
  return result;
}

}  // namespace ik
}  // namespace varietas

#endif
