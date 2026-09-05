// The trigonometric formulation, and its agreement with the half-angle one.
//
// Two rationalisations of the same map must produce the same geometry, and the
// tests here are arranged so that a discrepancy cannot hide behind a tolerance:
// everything runs over the rationals, and the correspondence between the two
// parameterisations is itself rational, since
//
//     cos q = (1 - t^2) / (1 + t^2),   sin q = 2t / (1 + t^2)
//
// carries a rational t to a rational point of the circle. So a configuration
// can be written in both formulations exactly, and the two forward kinematics
// maps must agree there identically rather than nearly.

#include <array>
#include <cstdint>
#include <vector>

#include <gtest/gtest.h>

#include "varietas/codegen/rational.hpp"
#include "varietas/core/order/grevlex.hpp"
#include "varietas/kinematics/chain.hpp"
#include "varietas/kinematics/rationalize.hpp"
#include "varietas/kinematics/trigonometric.hpp"
#include "varietas/kinematics/workspace.hpp"

namespace {

using varietas::chain;
using varietas::grevlex;
using varietas::rational;
using varietas::revolute_joint;
using varietas::rigid_transform;
using varietas::vector3;

using traits = varietas::coefficient_traits<rational>;

rational q(std::int64_t n, std::int64_t d = 1) { return varietas::make_rational(n, d); }

vector3<rational> point(std::int64_t x, std::int64_t y, std::int64_t z) {
  return vector3<rational>(q(x), q(y), q(z));
}

// The point of the circle the half-angle variable t names, exactly.
std::array<rational, 2> circle_point(const rational& t) {
  const rational denominator = traits::one() + t * t;
  const rational inverse = traits::inverse(denominator);
  return {(traits::one() - t * t) * inverse, (t + t) * inverse};
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

TEST(trigonometric, the_ring_has_two_variables_for_each_revolute_joint) {
  const chain<rational> robot = planar_two_link();
  EXPECT_EQ(varietas::trigonometric_variable_count(robot), 4u);

  const auto indices = varietas::cosine_variable_indices(robot);
  ASSERT_EQ(indices.size(), 2u);
  EXPECT_EQ(indices[0], 0u);
  EXPECT_EQ(indices[1], 2u);

  const auto relations = varietas::circle_relations<4, grevlex>(robot);
  EXPECT_EQ(relations.size(), 2u) << "one identity per revolute joint";
  for (const auto& r : relations) {
    EXPECT_EQ(r.degree(), 2);
  }
}

// The agreement that matters: the same configuration, written in both
// parameterisations, must give the same transform exactly.
TEST(trigonometric, agrees_with_the_half_angle_map_at_every_configuration) {
  const chain<rational> robot = torus_arm();

  const auto rational_map = varietas::rational_forward_kinematics<2, grevlex>(robot);
  const auto trig_map = varietas::trigonometric_forward_kinematics<4, grevlex>(robot);

  const std::array<std::array<rational, 2>, 5> samples = {{{q(0), q(0)},
                                                           {q(1), q(0)},
                                                           {q(0), q(1)},
                                                           {q(1, 2), q(1, 3)},
                                                           {q(-2), q(3, 5)}}};
  for (const auto& sample : samples) {
    const auto first = circle_point(sample[0]);
    const auto second = circle_point(sample[1]);

    const auto from_half_angle = rational_map.evaluate({sample[0], sample[1]});
    const auto from_trig =
        trig_map.evaluate({first[0], first[1], second[0], second[1]});

    for (std::size_t i = 0; i < 3; ++i) {
      for (std::size_t j = 0; j < 3; ++j) {
        EXPECT_EQ(from_half_angle.rotation()(i, j), from_trig.rotation()(i, j))
            << "rotation entry (" << i << ", " << j << ")";
      }
      EXPECT_EQ(from_half_angle.translation()[i], from_trig.translation()[i])
          << "translation entry " << i;
    }
  }
}

// On the circle the map lands in SE(3), and over the rationals that is an
// identity rather than a residual: the only admissible orthogonality defect is
// zero. Off the circle it does not, which is what the circle relations are for.
TEST(trigonometric, the_rotation_is_exactly_orthogonal_on_the_circle) {
  const chain<rational> robot = torus_arm();
  const auto map = varietas::trigonometric_forward_kinematics<4, grevlex>(robot);

  const auto first = circle_point(q(1, 2));
  const auto second = circle_point(q(-3));
  const auto on_circle = map.evaluate({first[0], first[1], second[0], second[1]});

  for (std::size_t i = 0; i < 3; ++i) {
    for (std::size_t j = 0; j < 3; ++j) {
      rational entry = traits::zero();
      for (std::size_t k = 0; k < 3; ++k) {
        entry = entry + on_circle.rotation()(k, i) * on_circle.rotation()(k, j);
      }
      const rational expected = (i == j) ? traits::one() : traits::zero();
      EXPECT_EQ(entry, expected) << "R^T R entry (" << i << ", " << j << ")";
    }
  }

  // A point off the circle is not a configuration and the result is not a
  // rotation, which is the statement that the relations carry real content.
  const auto off_circle = map.evaluate({q(2), q(2), first[0], first[1]});
  rational norm = traits::zero();
  for (std::size_t k = 0; k < 3; ++k) {
    norm = norm + off_circle.rotation()(k, 0) * off_circle.rotation()(k, 0);
  }
  EXPECT_NE(norm, traits::one());
}

// The circle relations vanish at a configuration and the pose residuals close
// the graph, so a reachable point together with the parameters that reach it is
// a common zero of every generator, exactly, over the rationals.
TEST(trigonometric, the_graph_generators_vanish_at_a_configuration) {
  const chain<rational> robot = torus_arm();
  const auto generators = varietas::trigonometric_workspace_generators<4>(robot);
  const auto map = varietas::trigonometric_forward_kinematics<4, grevlex>(robot);

  const auto first = circle_point(q(1, 3));
  const auto second = circle_point(q(2));
  const std::array<rational, 4> parameters = {first[0], first[1], second[0], second[1]};
  const auto reached = map.evaluate(parameters).translation();

  const std::array<rational, 7> at = {parameters[0], parameters[1], parameters[2],
                                      parameters[3], reached[0],    reached[1],
                                      reached[2]};
  for (const auto& g : generators) {
    EXPECT_TRUE(traits::is_zero(g.evaluate(at)))
        << "a point of the graph satisfies every generator";
  }
}

// The two formulations eliminate to the same ideal. This is the check that
// licenses moving implicitization onto the trigonometric path: the answer is
// unchanged, only the cost is.
TEST(trigonometric, implicitization_agrees_with_the_half_angle_path) {
  const chain<rational> robot = torus_arm();

  const auto from_trig = varietas::workspace_relations<4>(robot);
  const auto from_half_angle = varietas::workspace_relations_half_angle<2>(robot);

  EXPECT_EQ(from_trig, from_half_angle)
      << "the reduced basis of the elimination ideal does not depend on the "
         "parameterisation it was computed through";
}

}  // namespace
