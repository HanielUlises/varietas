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
// arm would stay planar after all, which is why the z = 0 check below is worth
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
// positive-dimensional, which is the case the dimension count has to catch and
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
// the one whose solve over Q(x, y, z) first finished once the gcd was
// replaced. Its point here is that it decouples.
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
// sign of the base angle recovered from the arctangent, one character in the
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
// still holds; it only moves the height the reduced problem is posed at.
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

inline rigid_transform<rational> at(rational x, rational y, rational z) {
  return rigid_transform<rational>::translation_only(vector3<rational>(x, y, z));
}

// A six-joint arm of the usual industrial shape: a base turning about z, a
// shoulder displaced 25 mm from the base axis and 400 mm up, an upper arm of
// 455 mm and a forearm of 420 mm with a 35 mm offset, all pitching about y,
// and a wrist rolling about x, pitching about y and rolling about x again,
// its three axes meeting in one point 80 mm short of the flange. The
// dimensions are those of a common six-kilogram arm; nothing depends on them
// beyond their being rationals.
inline chain<rational> industrial_six() {
  const auto I = rigid_transform<rational>::identity();
  chain<rational> robot("industrial_6r");
  robot.add_joint(revolute_joint<rational>("a1", vector3<rational>::unit(2), I));
  robot.add_joint(revolute_joint<rational>("a2", vector3<rational>::unit(1),
                                           at(rational(1, 40), nil(), rational(2, 5))));
  robot.add_joint(revolute_joint<rational>("a3", vector3<rational>::unit(1),
                                           at(nil(), nil(), rational(91, 200))));
  robot.add_joint(revolute_joint<rational>("a4", vector3<rational>::unit(0),
                                           at(rational(21, 50), nil(), rational(7, 200))));
  robot.add_joint(revolute_joint<rational>("a5", vector3<rational>::unit(1), I));
  robot.add_joint(revolute_joint<rational>("a6", vector3<rational>::unit(0), I));
  robot.set_tool(at(rational(2, 25), nil(), nil()));
  return robot;
}

// The same arm with the joint ranges of its description, so that the solvers'
// limit handling has ranges to respect, two of them narrower than a turn on
// one side and one wider than a turn.
inline chain<rational> industrial_six_limited() {
  const chain<rational> free = industrial_six();
  const double lower[6] = {-2.97, -3.32, -2.09, -3.23, -2.09, -6.11};
  const double upper[6] = {2.97, 0.79, 2.72, 3.23, 2.09, 6.11};
  chain<rational> robot("industrial_6r_limited");
  for (std::size_t i = 0; i < free.joints().size(); ++i) {
    auto j = free.joints()[i];
    j.has_limits = true;
    j.lower = lower[i];
    j.upper = upper[i];
    robot.add_joint(j);
  }
  robot.set_tool(free.tool());
  return robot;
}

// A wrist whose axes meet but are not orthogonal, on an arm whose base is
// displaced off its own axis, so that nothing about it is the textbook case:
// the decoupling refuses the arm, and the wrist's rotation has to be split
// along three skew directions through one point.
inline chain<rational> skewed_wrist_six() {
  const auto I = rigid_transform<rational>::identity();
  chain<rational> robot("skewed_wrist_6r");
  robot.add_joint(revolute_joint<rational>("j1", vector3<rational>::unit(2),
                                           at(rational(1, 10), nil(), nil())));
  robot.add_joint(revolute_joint<rational>("j2", vector3<rational>::unit(1),
                                           at(nil(), nil(), rational(1, 2))));
  robot.add_joint(revolute_joint<rational>("j3", vector3<rational>::unit(1),
                                           at(nil(), nil(), rational(3, 5))));
  robot.add_joint(revolute_joint<rational>("j4", vector3<rational>::unit(0),
                                           at(rational(1, 2), nil(), rational(1, 10))));
  robot.add_joint(revolute_joint<rational>(
      "j5", vector3<rational>(nil(), rational(3, 5), rational(4, 5)), I));
  robot.add_joint(revolute_joint<rational>(
      "j6", vector3<rational>(rational(2, 3), rational(1, 3), rational(2, 3)), I));
  robot.set_tool(at(rational(1, 10), rational(1, 20), nil()));
  return robot;
}

// The same arm with its last axis moved a centimetre off the wrist centre, so
// that the three wrist axes no longer meet: the decomposition must refuse it.
inline chain<rational> broken_wrist_six() {
  chain<rational> robot = industrial_six();
  chain<rational> moved("broken_wrist_6r");
  for (std::size_t i = 0; i < robot.joints().size(); ++i) {
    auto j = robot.joints()[i];
    if (i == 5) {
      j.origin = at(nil(), rational(1, 100), nil());
    }
    moved.add_joint(j);
  }
  moved.set_tool(robot.tool());
  return moved;
}

}  // namespace varietas_test

#endif
