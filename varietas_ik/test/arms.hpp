#ifndef VARIETAS_IK_TEST_ARMS_HPP
#define VARIETAS_IK_TEST_ARMS_HPP

#include "varietas/codegen/rational.hpp"
#include "varietas/kinematics/chain.hpp"
#include "varietas/kinematics/rigid_transform.hpp"

namespace varietas_test {

using varietas::chain;
using varietas::rational;
using varietas::revolute_joint;
using varietas::rigid_transform;
using varietas::vector3;

inline rational unit() { return rational(1); }
inline rational nil() { return rational(0); }

inline rigid_transform<rational> along_x(int length) {
  return rigid_transform<rational>::translation_only(
      vector3<rational>(rational(length), nil(), nil()));
}

// The planar two-link arm of unit links: the running example throughout the
// library, and the arm whose inverse kinematics is known without computing
// anything. Two elbow configurations for a reachable point, so dim_k A = 2 is
// the number the parametric solve has to return.
inline chain<rational> planar_two_link() {
  chain<rational> robot("planar_2r");
  robot.add_joint(revolute_joint<rational>("q1", vector3<rational>::unit(2),
                                           rigid_transform<rational>::identity()));
  robot.add_joint(revolute_joint<rational>("q2", vector3<rational>::unit(2), along_x(1)));
  robot.set_tool(along_x(1));
  return robot;
}

// The same arm with a third coplanar joint. Three joints against two position
// equations leaves a curve of configurations for each reachable point, which is
// what the bridge has to refuse rather than emit.
inline chain<rational> planar_three_link() {
  chain<rational> robot("planar_3r");
  robot.add_joint(revolute_joint<rational>("q1", vector3<rational>::unit(2),
                                           rigid_transform<rational>::identity()));
  robot.add_joint(revolute_joint<rational>("q2", vector3<rational>::unit(2), along_x(1)));
  robot.add_joint(revolute_joint<rational>("q3", vector3<rational>::unit(2), along_x(1)));
  robot.set_tool(along_x(1));
  return robot;
}

// A two-joint arm that genuinely leaves the z = 0 plane: the second axis is x
// and the tool offset is along y, so rotating the second joint swings the tool
// out of the plane. Asking for (x, y) alone would drop an equation that
// constrains the joints.
//
// The offset has to be off the second axis for this to be true. A tool offset
// along x would lie on that axis, rotation about it would move nothing, and the
// arm would stay planar after all — which is why the z = 0 check below is worth
// having rather than obvious.
inline chain<rational> out_of_plane_two_link() {
  chain<rational> robot("out_of_plane_2r");
  robot.add_joint(revolute_joint<rational>("q1", vector3<rational>::unit(2),
                                           rigid_transform<rational>::identity()));
  robot.add_joint(revolute_joint<rational>("q2", vector3<rational>::unit(0), along_x(1)));
  robot.set_tool(rigid_transform<rational>::translation_only(
      vector3<rational>(nil(), unit(), nil())));
  return robot;
}

// Two revolute joints about the same axis through the same point. Only their
// sum moves the tool, so the arm has two joints and one degree of freedom in
// the plane: the tool traces a circle, and every point on it is reached by a
// whole curve of configurations. Two unknowns and two equations, and still
// positive-dimensional — which is the case the dimension count has to catch and
// that no amount of counting generators would.
inline chain<rational> coincident_two_link() {
  chain<rational> robot("coincident_2r");
  robot.add_joint(revolute_joint<rational>("q1", vector3<rational>::unit(2),
                                           rigid_transform<rational>::identity()));
  robot.add_joint(revolute_joint<rational>("q2", vector3<rational>::unit(2),
                                           rigid_transform<rational>::identity()));
  robot.set_tool(along_x(1));
  return robot;
}

// The anthropomorphic arm: a base that yaws about z, then a shoulder and an
// elbow that both pitch about y. The textbook three-joint positioning arm, and
// the one that does not solve over Q(x, y, z) in any reasonable time — its
// whole point here is that it decouples.
inline chain<rational> anthropomorphic_three_link() {
  chain<rational> robot("anthropomorphic_3r");
  robot.add_joint(revolute_joint<rational>("q1", vector3<rational>::unit(2),
                                           rigid_transform<rational>::identity()));
  robot.add_joint(revolute_joint<rational>("q2", vector3<rational>::unit(1),
                                           rigid_transform<rational>::identity()));
  robot.add_joint(revolute_joint<rational>("q3", vector3<rational>::unit(1), along_x(1)));
  robot.set_tool(along_x(1));
  return robot;
}

// The same arm with its base axis pointing the other way.
//
// Turning about -z is turning about z backwards, so the decomposition holds
// exactly as before and the reduced problem is identical. What changes is the
// sign of the base angle recovered from the arctangent — one character in the
// emitted wrapper, and the only thing in the decoupling that nothing else
// exercises.
inline chain<rational> anthropomorphic_reversed_base() {
  chain<rational> robot("anthropomorphic_3r_reversed");
  robot.add_joint(revolute_joint<rational>(
      "q1", vector3<rational>(nil(), nil(), rational(-1)),
      rigid_transform<rational>::identity()));
  robot.add_joint(revolute_joint<rational>("q2", vector3<rational>::unit(1),
                                           rigid_transform<rational>::identity()));
  robot.add_joint(revolute_joint<rational>("q3", vector3<rational>::unit(1), along_x(1)));
  robot.set_tool(along_x(1));
  return robot;
}

// The same arm on a pedestal. The base placement is a translation along the
// axis it turns about, so it commutes with that rotation and the decomposition
// still holds — it only moves the height the reduced problem is posed at.
inline chain<rational> anthropomorphic_on_a_pedestal() {
  chain<rational> robot("anthropomorphic_3r_pedestal");
  robot.add_joint(revolute_joint<rational>(
      "q1", vector3<rational>::unit(2),
      rigid_transform<rational>::translation_only(
          vector3<rational>(nil(), nil(), rational(2)))));
  robot.add_joint(revolute_joint<rational>("q2", vector3<rational>::unit(1),
                                           rigid_transform<rational>::identity()));
  robot.add_joint(revolute_joint<rational>("q3", vector3<rational>::unit(1), along_x(1)));
  robot.set_tool(along_x(1));
  return robot;
}

// The same arm with its base shifted sideways, off its own axis. The placement
// no longer commutes with the rotation, the axis itself is carried around, and
// there is no fixed plane to sweep.
inline chain<rational> anthropomorphic_off_axis() {
  chain<rational> robot("anthropomorphic_3r_offset");
  robot.add_joint(revolute_joint<rational>("q1", vector3<rational>::unit(2), along_x(1)));
  robot.add_joint(revolute_joint<rational>("q2", vector3<rational>::unit(1),
                                           rigid_transform<rational>::identity()));
  robot.add_joint(revolute_joint<rational>("q3", vector3<rational>::unit(1), along_x(1)));
  robot.set_tool(along_x(1));
  return robot;
}

}  // namespace varietas_test

#endif
