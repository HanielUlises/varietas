#ifndef VARIETAS_KINEMATICS_EVALUATE_HPP
#define VARIETAS_KINEMATICS_EVALUATE_HPP

#include <cmath>
#include <cstddef>
#include <vector>

#include "varietas/core/config.hpp"
#include "varietas/kinematics/chain.hpp"
#include "varietas/kinematics/rigid_transform.hpp"

namespace varietas {

// Numerical forward kinematics: the map evaluated at a configuration, as
// against the polynomial system that represents it.
//
// This is deliberately the floating point side of the library. Evaluating the
// pose at a given angle is transcendental — it needs cos and sin of a real
// number, which no rational field supplies — so nothing here pretends to
// exactness. Its purpose is to be the reference the algebra is checked against:
// a solution recovered from the variety is correct exactly when substituting it
// back through this map reproduces the requested pose, and the rationalisation
// to come has to agree with it at every configuration.

// The coefficient field of a chain, changed by way of the traits' bridge. The
// offline pipeline builds a chain over the exact field and evaluates over
// double; this is where the two meet, and, being a rounding, it is the reason
// evaluation is a check on the algebra and not a substitute for it.
template <class To, class From>
chain<To> chain_cast(const chain<From>& source) {
  using from_traits = coefficient_traits<From>;
  using to_traits = coefficient_traits<To>;

  const auto convert_vector = [](const vector3<From>& v) {
    return vector3<To>(to_traits::from_double(from_traits::to_double(v[0])),
                       to_traits::from_double(from_traits::to_double(v[1])),
                       to_traits::from_double(from_traits::to_double(v[2])));
  };
  const auto convert_transform = [&convert_vector](const rigid_transform<From>& t) {
    matrix3<To> r;
    for (std::size_t i = 0; i < 3; ++i) {
      for (std::size_t j = 0; j < 3; ++j) {
        r(i, j) = to_traits::from_double(from_traits::to_double(t.rotation()(i, j)));
      }
    }
    return rigid_transform<To>(r, convert_vector(t.translation()));
  };

  chain<To> result(source.name());
  for (const joint<From>& j : source.joints()) {
    joint<To> converted;
    converted.name = j.name;
    converted.type = j.type;
    converted.axis = convert_vector(j.axis);
    converted.origin = convert_transform(j.origin);
    converted.has_limits = j.has_limits;
    converted.lower = j.lower;
    converted.upper = j.upper;
    result.add_joint(std::move(converted));
  }
  result.set_tool(convert_transform(source.tool()));
  return result;
}

// The displacement a joint contributes at the given value: a rotation about its
// axis, a translation along it, or nothing at all.
inline rigid_transform<double> joint_displacement(const joint<double>& j, double value) {
  switch (j.type) {
    case joint_type::revolute:
      return rigid_transform<double>::rotation_only(
          rotation_about_axis<double>(j.axis, std::cos(value), std::sin(value)));
    case joint_type::prismatic:
      return rigid_transform<double>::translation_only(value * j.axis);
    case joint_type::fixed:
      break;
  }
  return rigid_transform<double>::identity();
}

// The frame of every link, base first, tip last, so that the whole chain can be
// drawn or compared frame by frame rather than only at the end effector. The
// returned vector has one entry per joint plus a leading base frame and a
// trailing tool frame.
inline std::vector<rigid_transform<double>> link_frames(
    const chain<double>& robot, const std::vector<double>& values) {
  VARIETAS_ASSERT(values.size() == robot.degrees_of_freedom());

  std::vector<rigid_transform<double>> frames;
  frames.reserve(robot.size() + 2);

  rigid_transform<double> running = rigid_transform<double>::identity();
  frames.push_back(running);

  std::size_t variable = 0;
  for (const joint<double>& j : robot.joints()) {
    const double value = j.is_actuated() ? values[variable++] : 0.0;
    running = running * j.origin * joint_displacement(j, value);
    frames.push_back(running);
  }

  frames.push_back(running * robot.tool());
  return frames;
}

inline rigid_transform<double> forward_kinematics(const chain<double>& robot,
                                                 const std::vector<double>& values) {
  return link_frames(robot, values).back();
}

}  // namespace varietas

#endif
