// Implicitization of the reachable workspace.
//
// The three cases here are chosen because the construction behaves differently
// in each, and conflating them is how implicitization comes to be oversold. A
// one-joint arm traces a curve and elimination returns it exactly. A planar
// two-joint arm fills its plane, and elimination returns only the plane:
// the annulus it actually reaches is semialgebraic and no ideal describes it.
// A two-joint arm with perpendicular axes traces a surface in three-space, and
// there elimination earns its name.

#include <array>
#include <cstdint>
#include <vector>

#include <gtest/gtest.h>

#include "varietas/codegen/rational.hpp"
#include "varietas/core/order/grevlex.hpp"
#include "varietas/kinematics/chain.hpp"
#include "varietas/kinematics/rationalize.hpp"
#include "varietas/kinematics/workspace.hpp"

namespace {

using varietas::chain;
using varietas::grevlex;
using varietas::rational;
using varietas::revolute_joint;
using varietas::rigid_transform;
using varietas::vector3;

using relation = varietas::polynomial<rational, 3, grevlex>;

rational q(std::int64_t n, std::int64_t d = 1) { return varietas::make_rational(n, d); }

vector3<rational> point(std::int64_t x, std::int64_t y, std::int64_t z) {
  return vector3<rational>(q(x), q(y), q(z));
}

relation x() { return relation::variable(0); }
relation y() { return relation::variable(1); }
relation z() { return relation::variable(2); }
relation constant(std::int64_t n) { return relation::constant(q(n)); }

// A single revolute joint about z with a unit tool offset: the tool traces the
// unit circle in the plane z = 0.
chain<rational> one_link() {
  chain<rational> robot("planar_1r");
  robot.add_joint(revolute_joint<rational>("q1", vector3<rational>::unit(2),
                                           rigid_transform<rational>::identity()));
  robot.set_tool(rigid_transform<rational>::translation_only(point(1, 0, 0)));
  return robot;
}

chain<rational> planar_two_link() {
  chain<rational> robot("planar_2r");
  robot.add_joint(revolute_joint<rational>("q1", vector3<rational>::unit(2),
                                           rigid_transform<rational>::identity()));
  robot.add_joint(revolute_joint<rational>(
      "q2", vector3<rational>::unit(2),
      rigid_transform<rational>::translation_only(point(1, 0, 0))));
  robot.set_tool(rigid_transform<rational>::translation_only(point(1, 0, 0)));
  return robot;
}

// First joint about z, then two units along x, then a joint about y with a unit
// tool offset: the tool traces a torus of major radius two and minor radius one.
chain<rational> torus_arm() {
  chain<rational> robot("torus_2r");
  robot.add_joint(revolute_joint<rational>("q1", vector3<rational>::unit(2),
                                           rigid_transform<rational>::identity()));
  robot.add_joint(revolute_joint<rational>(
      "q2", vector3<rational>::unit(1),
      rigid_transform<rational>::translation_only(point(2, 0, 0))));
  robot.set_tool(rigid_transform<rational>::translation_only(point(1, 0, 0)));
  return robot;
}

// Every relation must vanish at every reachable point, which is the defining
// property and the one worth checking directly rather than only comparing
// bases.
void expect_vanishes_on_workspace(const std::vector<relation>& relations,
                                  const chain<rational>& robot, std::size_t joints) {
  const std::array<std::array<rational, 2>, 5> samples = {{{q(0), q(0)},
                                                           {q(1), q(0)},
                                                           {q(0), q(1)},
                                                           {q(1, 2), q(1, 3)},
                                                           {q(-2), q(3, 5)}}};
  for (const auto& sample : samples) {
    vector3<rational> reached;
    if (joints == 1) {
      const auto map = varietas::rational_forward_kinematics<1, grevlex>(robot);
      reached = map.evaluate({sample[0]}).translation();
    } else {
      const auto map = varietas::rational_forward_kinematics<2, grevlex>(robot);
      reached = map.evaluate({sample[0], sample[1]}).translation();
    }

    const std::array<rational, 3> at = {reached[0], reached[1], reached[2]};
    for (const relation& r : relations) {
      EXPECT_TRUE(varietas::coefficient_traits<rational>::is_zero(r.evaluate(at)))
          << "a reachable point must satisfy every relation";
    }
  }
}

TEST(workspace, a_one_joint_arm_traces_a_circle) {
  const chain<rational> robot = one_link();
  const auto relations = varietas::workspace_relations<2>(robot);

  EXPECT_FALSE(varietas::workspace_is_dense(relations));

  // z = 0 and x^2 + y^2 = 1, which is the circle, exactly.
  const varietas::ideal<rational, 3, grevlex> expected(
      {z(), x() * x() + y() * y() - constant(1)});
  EXPECT_EQ(relations, expected.basis());

  expect_vanishes_on_workspace(relations, robot, 1);
}

// The case where the construction has nothing to say, and says so. The arm
// reaches the annulus 0 <= r <= 2, whose Zariski closure is the whole plane, so
// the only relation that survives is the one that holds identically.
TEST(workspace, a_planar_two_joint_arm_yields_only_its_plane) {
  const chain<rational> robot = planar_two_link();
  const auto relations = varietas::workspace_relations<4>(robot);

  const varietas::ideal<rational, 3, grevlex> expected({z()});
  EXPECT_EQ(relations, expected.basis());

  // The boundary of the reachable set is at r = 2, and nothing in the ideal
  // knows it: the unreachable origin-adjacent points and the far ones alike
  // satisfy z = 0. Reachability is not an algebraic condition.
  const std::array<rational, 3> far_outside = {q(100), q(0), q(0)};
  for (const relation& r : relations) {
    EXPECT_TRUE(
        varietas::coefficient_traits<rational>::is_zero(r.evaluate(far_outside)))
        << "a point far outside the reach still satisfies the closure's equations";
  }

  expect_vanishes_on_workspace(relations, robot, 2);
}

// And the case it exists for: the map is not dominant, the closure is a genuine
// surface, and elimination produces the quartic that cuts it out.
TEST(workspace, perpendicular_axes_trace_a_torus) {
  const chain<rational> robot = torus_arm();
  const auto relations = varietas::workspace_relations<4>(robot);

  EXPECT_FALSE(varietas::workspace_is_dense(relations));
  ASSERT_EQ(relations.size(), 1u) << "a hypersurface has one equation";

  // (x^2 + y^2 + z^2 + R^2 - r^2)^2 = 4 R^2 (x^2 + y^2), with R = 2 and r = 1.
  const relation squared_norm = x() * x() + y() * y() + z() * z();
  const relation torus = (squared_norm + constant(3)) * (squared_norm + constant(3)) -
                         constant(16) * (x() * x() + y() * y());

  EXPECT_EQ(relations.front(), torus.monic());
  EXPECT_EQ(relations.front().degree(), 4);

  expect_vanishes_on_workspace(relations, robot, 2);

  // A point off the torus is not in the closure, which is the statement that
  // the relation is not vacuous.
  const std::array<rational, 3> off = {q(0), q(0), q(0)};
  EXPECT_FALSE(
      varietas::coefficient_traits<rational>::is_zero(relations.front().evaluate(off)));
}

}  // namespace
