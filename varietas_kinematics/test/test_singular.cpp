// The singular locus of an arm, as a variety rather than as a condition number.
//
// The arms here are the smallest ones on which the answer is known in closed
// form, so that every assertion is against a fact of kinematics and not against
// whatever the last run happened to print. The planar two-link arm is singular
// exactly when its elbow is straight or folded, and the two cases are the outer
// boundary of its annulus and the hole in the middle of it. The torus arm is
// the interesting one: it has no real singularity at all, and what the ideal
// says about it is a lesson in the difference between a variety and its real
// points.

#include <array>
#include <cmath>
#include <cstdint>
#include <vector>

#include <gtest/gtest.h>

#include "varietas/codegen/rational.hpp"
#include "varietas/core/ideal/dimension.hpp"
#include "varietas/core/ideal/ideal.hpp"
#include "varietas/core/ideal/saturation.hpp"
#include "varietas/core/order/grevlex.hpp"
#include "varietas/kinematics/chain.hpp"
#include "varietas/kinematics/evaluate.hpp"
#include "varietas/kinematics/singular.hpp"

namespace {

using varietas::chain;
using varietas::grevlex;
using varietas::rational;
using varietas::revolute_joint;
using varietas::rigid_transform;
using varietas::vector3;

using relation = varietas::polynomial<rational, 3, grevlex>;

template <std::size_t V>
using joint_polynomial = varietas::polynomial<rational, V, grevlex>;

rational q(std::int64_t n, std::int64_t d = 1) { return varietas::make_rational(n, d); }

vector3<rational> point(std::int64_t x, std::int64_t y, std::int64_t z) {
  return vector3<rational>(q(x), q(y), q(z));
}

relation x() { return relation::variable(0); }
relation y() { return relation::variable(1); }
relation z() { return relation::variable(2); }
relation constant(std::int64_t n) { return relation::constant(q(n)); }

// Two unit links about parallel z axes: the classical elbow.
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

// Three unit links about parallel z axes, which is the smallest arm that is
// redundant for one task and not for another.
chain<rational> planar_three_link() {
  chain<rational> robot("planar_3r");
  robot.add_joint(revolute_joint<rational>("q1", vector3<rational>::unit(2),
                                           rigid_transform<rational>::identity()));
  robot.add_joint(revolute_joint<rational>(
      "q2", vector3<rational>::unit(2),
      rigid_transform<rational>::translation_only(point(1, 0, 0))));
  robot.add_joint(revolute_joint<rational>(
      "q3", vector3<rational>::unit(2),
      rigid_transform<rational>::translation_only(point(1, 0, 0))));
  robot.set_tool(rigid_transform<rational>::translation_only(point(1, 0, 0)));
  return robot;
}

// A joint about z, an offset of `major` along x, a joint about y, and a unit
// tool: the tool traces a torus of the given major radius and minor radius one.
// At major radius two the torus is an ordinary one; at major radius one it is
// pinched, the inner circle having collapsed to the origin.
chain<rational> torus_arm(std::int64_t major) {
  chain<rational> robot("torus_2r");
  robot.add_joint(revolute_joint<rational>("q1", vector3<rational>::unit(2),
                                           rigid_transform<rational>::identity()));
  robot.add_joint(revolute_joint<rational>(
      "q2", vector3<rational>::unit(1),
      rigid_transform<rational>::translation_only(point(major, 0, 0))));
  robot.set_tool(rigid_transform<rational>::translation_only(point(1, 0, 0)));
  return robot;
}

// c_k and s_k of the k-th revolute joint, as ring variables.
template <std::size_t V>
joint_polynomial<V> cosine(std::size_t k) {
  return joint_polynomial<V>::variable(2 * k);
}

template <std::size_t V>
joint_polynomial<V> sine(std::size_t k) {
  return joint_polynomial<V>::variable(2 * k + 1);
}

template <std::size_t V>
std::vector<joint_polynomial<V>> circles(std::size_t joints) {
  std::vector<joint_polynomial<V>> relations;
  for (std::size_t k = 0; k < joints; ++k) {
    relations.push_back(cosine<V>(k) * cosine<V>(k) + sine<V>(k) * sine<V>(k) -
                        joint_polynomial<V>::constant(q(1)));
  }
  return relations;
}

// --- the Jacobian ----------------------------------------------------------

TEST(singular, the_jacobian_of_a_two_link_arm_is_the_textbook_one) {
  // For a planar two-link arm with unit links the determinant of the position
  // Jacobian is l1 l2 sin q2. The rows out of the plane contribute nothing, so
  // the two other maximal minors vanish identically and are dropped.
  const auto jacobian = varietas::trigonometric_jacobian<4, grevlex>(planar_two_link());
  ASSERT_EQ(jacobian.rows(), 6u);
  ASSERT_EQ(jacobian.cols(), 2u);

  const auto minors = varietas::maximal_minors(
      jacobian.submatrix(varietas::position_rows(), {0, 1}));
  ASSERT_EQ(minors.size(), 1u);

  // In the ambient ring the minor is not sin q2 but (c1² + s1²)·s2. The
  // determinant is a polynomial identity and knows nothing of the circle; the
  // textbook formula is what it becomes on the parameter variety, where
  // c1² + s1² is one. Reducing modulo the circle relations — which are
  // generators of the singular ideal precisely so that this reduction is
  // legitimate — returns the sine of the elbow angle exactly.
  const joint_polynomial<4> circle_factor =
      cosine<4>(0) * cosine<4>(0) + sine<4>(0) * sine<4>(0);
  EXPECT_EQ(minors.front(), circle_factor * sine<4>(1));
  EXPECT_EQ(varietas::normal_form(minors.front(), circles<4>(2)), sine<4>(1))
      << "on the parameter variety the minor is the sine of the elbow angle";
}

TEST(singular, the_polynomial_jacobian_agrees_with_finite_differences) {
  // The construction above is a formula, and a formula is worth checking
  // against the map it claims to differentiate. Central differences of the
  // forward kinematics are an independent statement — nothing they use came
  // from trigonometric.hpp — so agreement to the order of the step is real
  // evidence and not a tautology.
  const chain<double> robot = varietas::chain_cast<double>(planar_three_link());
  const auto jacobian = varietas::trigonometric_jacobian<6, grevlex>(robot);

  const std::array<std::vector<double>, 3> configurations = {
      {{0.3, 0.7, -0.4}, {1.2, -0.2, 2.0}, {-0.9, 1.6, 0.5}}};

  for (const auto& angles : configurations) {
    std::array<double, 6> at{};
    for (std::size_t k = 0; k < 3; ++k) {
      at[2 * k] = std::cos(angles[k]);
      at[2 * k + 1] = std::sin(angles[k]);
    }

    const double step = 1e-6;
    for (std::size_t column = 0; column < 3; ++column) {
      std::vector<double> forward = angles;
      std::vector<double> backward = angles;
      forward[column] += step;
      backward[column] -= step;

      const auto ahead = varietas::forward_kinematics(robot, forward).translation();
      const auto behind = varietas::forward_kinematics(robot, backward).translation();

      for (std::size_t row = 0; row < 3; ++row) {
        const double difference = (ahead[row] - behind[row]) / (2.0 * step);
        EXPECT_NEAR(jacobian(row, column).evaluate(at), difference, 1e-7)
            << "row " << row << ", column " << column;
      }
    }
  }
}

// --- the singular locus as a set of configurations -------------------------

TEST(singular, a_two_link_arm_is_singular_exactly_where_its_elbow_is_straight_or_folded) {
  const auto basis = varietas::singular_ideal<4>(planar_two_link(), varietas::position_rows());

  std::vector<joint_polynomial<4>> expected = circles<4>(2);
  expected.push_back(sine<4>(1));
  EXPECT_EQ(basis, (varietas::ideal<rational, 4, grevlex>(expected).basis()));

  // A circle's worth of configurations: the first joint is free to turn along
  // the singular set, the second is pinned to one of two values.
  const auto dimension = varietas::ideal_dimension(basis);
  EXPECT_FALSE(dimension.is_empty);
  EXPECT_EQ(dimension.dimension, 1u);
  EXPECT_FALSE(dimension.independent[2]);
  EXPECT_FALSE(dimension.independent[3])
      << "the elbow variables are constrained on the singular set";
}

TEST(singular, the_task_decides_what_singular_means) {
  // A planar three-link arm cannot move its tool out of its plane at any
  // configuration, so for the three-dimensional position task its Jacobian has
  // rank at most two everywhere and every configuration is singular. That is
  // not a defect of the construction; it is the correct answer to the question
  // asked, and the reason the rows are the caller's choice.
  const auto everywhere =
      varietas::singular_ideal<6>(planar_three_link(), varietas::position_rows());
  EXPECT_EQ(everywhere,
            (varietas::ideal<rational, 6, grevlex>(circles<6>(3)).basis()))
      << "for a task it cannot span, the arm is singular on its whole parameter space";
  EXPECT_EQ(varietas::ideal_dimension(everywhere).dimension, 3u);

  // Asked for the planar pose it can span — position in the plane and heading
  // within it — the same arm is singular on a surface, cut out by the elbow.
  const auto planar =
      varietas::singular_ideal<6>(planar_three_link(), varietas::planar_pose_rows());
  std::vector<joint_polynomial<6>> expected = circles<6>(3);
  expected.push_back(sine<6>(1));
  EXPECT_EQ(planar, (varietas::ideal<rational, 6, grevlex>(expected).basis()))
      << "the wrist adds no rank, so the singularity is still the elbow alone";
  EXPECT_EQ(varietas::ideal_dimension(planar).dimension, 2u);
}

// --- splitting into branches -----------------------------------------------

TEST(singular, splitting_separates_the_straight_elbow_from_the_folded_one) {
  // The singular set of the planar arm is two disjoint pieces, and no dimension
  // or Gröbner basis distinguishes them: the ideal contains s2, and s2 = 0 is
  // c2 = 1 or c2 = -1. Splitting along c2 - 1 separates them exactly, the two
  // branches being the arm straight and the arm folded back on itself.
  const auto basis = varietas::singular_ideal<6>(planar_three_link(),
                                                 varietas::planar_pose_rows());
  const auto split = varietas::split_along(
      basis, cosine<6>(1) - joint_polynomial<6>::constant(q(1)));

  const varietas::ideal<rational, 6, grevlex> straight(split.on);
  const varietas::ideal<rational, 6, grevlex> folded(split.away);

  EXPECT_TRUE(straight.contains(cosine<6>(1) - joint_polynomial<6>::constant(q(1))));
  EXPECT_TRUE(straight.contains(sine<6>(1)));
  EXPECT_TRUE(folded.contains(cosine<6>(1) + joint_polynomial<6>::constant(q(1))));
  EXPECT_TRUE(folded.contains(sine<6>(1)));

  EXPECT_FALSE(straight.contains(cosine<6>(1) + joint_polynomial<6>::constant(q(1))))
      << "the branches are disjoint, not two descriptions of the same set";

  // Each branch is a torus in the two free joints, so each has the dimension
  // the whole singular set had; what splitting recovers is not size but the
  // fact that there are two pieces.
  EXPECT_EQ(straight.dimension().dimension, 2u);
  EXPECT_EQ(folded.dimension().dimension, 2u);

  // Exhaustive as well as exact: every element of the original ideal is in both
  // branches, so V(straight) and V(folded) are subsets of the singular set, and
  // the splitting identity says their union is all of it.
  for (const auto& g : basis) {
    EXPECT_TRUE(straight.contains(g));
    EXPECT_TRUE(folded.contains(g));
  }
}

// --- the image in the workspace --------------------------------------------

TEST(singular, the_singular_image_of_a_two_link_arm_is_the_boundary_of_its_annulus) {
  const auto relations =
      varietas::singular_workspace_relations<4>(planar_two_link(), varietas::position_rows());

  // The arm is singular with the elbow straight, where the tool is on the
  // circle of radius two, and folded, where it is at the origin. The Zariski
  // closure of the union of the two is cut out by z together with the products
  // x·(x² + y² - 4) and y·(x² + y² - 4).
  const relation squared_radius = x() * x() + y() * y();
  const varietas::ideal<rational, 3, grevlex> expected(
      {z(), x() * (squared_radius - constant(4)), y() * (squared_radius - constant(4))});
  EXPECT_EQ(relations, expected.basis());

  // The defining property, checked directly at the two singular poses.
  const std::array<std::array<rational, 3>, 2> singular_poses = {
      {{q(2), q(0), q(0)}, {q(0), q(0), q(0)}}};
  for (const auto& pose : singular_poses) {
    for (const relation& r : relations) {
      EXPECT_TRUE(varietas::coefficient_traits<rational>::is_zero(r.evaluate(pose)));
    }
  }

  // A reachable but nonsingular pose does not satisfy them, which is what makes
  // the answer a boundary rather than a restatement of the workspace.
  const std::array<rational, 3> interior = {q(1), q(1), q(0)};
  bool some_relation_fails = false;
  for (const relation& r : relations) {
    some_relation_fails =
        some_relation_fails ||
        !varietas::coefficient_traits<rational>::is_zero(r.evaluate(interior));
  }
  EXPECT_TRUE(some_relation_fails);
}

TEST(singular, a_torus_arm_has_singularities_only_at_complex_configurations) {
  // The tool of the torus arm is at distance 2 + cos q2 from the z axis, and
  // the arm loses rank exactly when that distance vanishes. Over the reals it
  // never does, so the arm is nowhere singular — but the ideal is not the unit
  // ideal, because over the closure cos q2 = -2 solves it perfectly well, with
  // a sine of -3 under the square root.
  const auto basis = varietas::singular_ideal<4>(torus_arm(2), varietas::position_rows());
  EXPECT_FALSE(varietas::is_unit_ideal(basis))
      << "the variety is nonempty over the algebraic closure";

  const varietas::ideal<rational, 4, grevlex> singular(basis);
  EXPECT_TRUE(singular.contains(cosine<4>(1) + joint_polynomial<4>::constant(q(2))));
  EXPECT_TRUE(singular.contains(sine<4>(1) * sine<4>(1) + joint_polynomial<4>::constant(q(3))))
      << "s2² = -3, which is where the configuration leaves the reals";

  // The image says the same thing in the workspace, and says it in a form that
  // can be read at a glance: x² + y² = 0 and z² = -3 have no real solution, so
  // the arm meets no singularity anywhere in its motion. That an ideal cannot
  // decide this on its own is the same fact met over unreachable poses —
  // reality is a property of points, not of ideals.
  const auto relations =
      varietas::singular_workspace_relations<4>(torus_arm(2), varietas::position_rows());
  const varietas::ideal<rational, 3, grevlex> image(relations);
  EXPECT_TRUE(image.contains(x() * x() + y() * y()));
  EXPECT_TRUE(image.contains(z() * z() + constant(3)));
}

TEST(singular, a_pinched_torus_arm_is_singular_where_its_inner_circle_collapsed) {
  // With the offset equal to the tool length the inner circle of the torus
  // closes to a point, and the configurations that reach it are singular: the
  // first joint moves the tool nowhere, because the tool is on its axis.
  const auto basis = varietas::singular_ideal<4>(torus_arm(1), varietas::position_rows());

  std::vector<joint_polynomial<4>> expected = circles<4>(2);
  expected.push_back(cosine<4>(1) + joint_polynomial<4>::constant(q(1)));
  EXPECT_EQ(basis, (varietas::ideal<rational, 4, grevlex>(expected).basis()));

  // A circle of singular configurations — the first joint is free — all of them
  // mapping to the single point where the torus is pinched.
  EXPECT_EQ(varietas::ideal_dimension(basis).dimension, 1u);

  // Every singular configuration puts the tool at the origin, and elimination
  // returns x, y and z² rather than x, y and z. The variety is the same single
  // point — that is what the Closure Theorem promises and all it promises —
  // but the ideal is not radical, and the square is not noise: the arm reaches
  // the pinch point tangentially in z, folding back along the axis instead of
  // crossing it, and the elimination ideal records the order of contact that
  // the set alone has forgotten.
  const auto relations =
      varietas::singular_workspace_relations<4>(torus_arm(1), varietas::position_rows());
  const varietas::ideal<rational, 3, grevlex> expected_image({x(), y(), z() * z()});
  EXPECT_EQ(relations, expected_image.basis());

  const varietas::ideal<rational, 3, grevlex> image(relations);
  EXPECT_FALSE(image.contains(z()))
      << "the image ideal is not radical, and saying so is more honest than "
         "quietly taking a radical the library cannot compute";

  const std::array<rational, 3> origin = {q(0), q(0), q(0)};
  for (const relation& r : relations) {
    EXPECT_TRUE(varietas::coefficient_traits<rational>::is_zero(r.evaluate(origin)));
  }
}

}  // namespace
