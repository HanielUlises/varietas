// The tangent half-angle substitution: forward kinematics as a rational map,
// and the inverse kinematics ideal it produces.
//
// The forward direction is checked against the numerical map, which is what the
// substitution has to reproduce. The inverse direction is checked on the planar
// two-link arm, whose solutions are known in closed form, and over the exact
// field: a Gröbner basis is only as meaningful as the field it was computed in,
// and the whole point of the offline pipeline is that this one is exact.

#include <array>
#include <cmath>
#include <cstdint>
#include <random>
#include <vector>

#include <gtest/gtest.h>

#include "varietas/codegen/rational.hpp"
#include "varietas/core/ideal/ideal.hpp"
#include "varietas/core/order/grevlex.hpp"
#include "varietas/core/solve/spectral.hpp"
#include "varietas/kinematics/chain.hpp"
#include "varietas/kinematics/evaluate.hpp"
#include "varietas/kinematics/rationalize.hpp"

namespace {

using varietas::chain;
using varietas::grevlex;
using varietas::rational;
using varietas::rational_transform;
using varietas::revolute_joint;
using varietas::rigid_transform;
using varietas::vector3;

template <class Coeff>
Coeff one() {
  return varietas::coefficient_traits<Coeff>::one();
}

// The planar two-link arm of unit links: the running example, and the case
// whose inverse kinematics is known without computing anything.
template <class Coeff>
chain<Coeff> planar_two_link() {
  const Coeff unit = one<Coeff>();
  const Coeff zero = varietas::coefficient_traits<Coeff>::zero();

  chain<Coeff> robot("planar_2r");
  robot.add_joint(revolute_joint<Coeff>("q1", vector3<Coeff>::unit(2),
                                        rigid_transform<Coeff>::identity()));
  robot.add_joint(revolute_joint<Coeff>(
      "q2", vector3<Coeff>::unit(2),
      rigid_transform<Coeff>::translation_only(vector3<Coeff>(unit, zero, zero))));
  robot.set_tool(
      rigid_transform<Coeff>::translation_only(vector3<Coeff>(unit, zero, zero)));
  return robot;
}

// --- the rational map reproduces the numerical one -------------------------

TEST(rationalize, agrees_with_the_numerical_forward_kinematics) {
  const chain<double> robot = planar_two_link<double>();
  const auto map = varietas::rational_forward_kinematics<2, grevlex>(robot);

  std::mt19937 generator(20260817);
  std::uniform_real_distribution<double> uniform(-3.0, 3.0);

  double worst = 0.0;
  for (int trial = 0; trial < 200; ++trial) {
    // Sample the angles, and substitute t = tan(q / 2): the rational map is
    // parameterised by t, the numerical one by q, and the substitution is the
    // only thing relating them.
    const double q0 = uniform(generator);
    const double q1 = uniform(generator);
    const std::array<double, 2> point = {varietas::variable_from_angle(q0),
                                         varietas::variable_from_angle(q1)};

    const rigid_transform<double> ours = map.evaluate(point);
    const rigid_transform<double> reference =
        varietas::forward_kinematics(robot, {q0, q1});

    for (std::size_t i = 0; i < 3; ++i) {
      worst = std::max(worst,
                       std::fabs(ours.translation()[i] - reference.translation()[i]));
      for (std::size_t j = 0; j < 3; ++j) {
        worst = std::max(worst, std::fabs(ours.rotation()(i, j) -
                                          reference.rotation()(i, j)));
      }
    }
  }
  EXPECT_LT(worst, 1e-12);
}

TEST(rationalize, one_denominator_per_revolute_joint) {
  const chain<rational> robot = planar_two_link<rational>();
  const auto map = varietas::rational_forward_kinematics<2, grevlex>(robot);

  // Each revolute joint contributes exactly one factor 1 + t_i^2, whatever the
  // fixed geometry between them: the closure property the representation rests
  // on would be false otherwise.
  const auto& exponents = map.denominator_exponents();
  EXPECT_EQ(exponents[0], 1);
  EXPECT_EQ(exponents[1], 1);

  // (1 + t0^2)(1 + t1^2), of degree four.
  EXPECT_EQ(map.denominator().degree(), 4);
}

TEST(rationalize, the_numerator_is_a_rotation_once_the_denominator_is_cleared) {
  // R R^T = I holds for the rational map at every configuration, which over the
  // exact field is an exact statement and not a numerical one.
  const chain<rational> robot = planar_two_link<rational>();
  const auto map = varietas::rational_forward_kinematics<2, grevlex>(robot);

  const std::array<rational, 2> point = {varietas::make_rational(1, 2),
                                         varietas::make_rational(1, 3)};
  const rigid_transform<rational> value = map.evaluate(point);

  EXPECT_EQ(orthogonality_defect(value.rotation()), 0.0);
  EXPECT_EQ(value.rotation().determinant(), varietas::make_rational(1));
}

// --- the inverse kinematics ideal ------------------------------------------

// A pose reachable by construction: evaluate the map at rational t, and the
// tool position is rational too, so the whole problem stays in the exact field.
TEST(rationalize, recovers_the_two_elbow_solutions) {
  const chain<rational> robot = planar_two_link<rational>();
  const auto map = varietas::rational_forward_kinematics<2, grevlex>(robot);

  const std::array<rational, 2> known = {varietas::make_rational(1, 2),
                                         varietas::make_rational(1, 3)};
  const vector3<rational> target = map.evaluate(known).translation();

  const auto generators =
      varietas::position_ideal_generators<2, grevlex>(robot, target);
  const varietas::ideal<rational, 2, grevlex> solutions(generators);

  ASSERT_TRUE(solutions.is_zero_dimensional());

  // Two elbow configurations, which is what a planar two-link arm has at a
  // reachable point in the interior of its workspace. The count is read off the
  // leading monomials of the basis, and it is the certificate that no branch
  // was missed, not a claim that two were found.
  EXPECT_EQ(solutions.quotient().monomials.size(), 2u);

  // The configuration the target was built from is a solution, exactly. Over
  // the rationals this is a structural test and not a residual below a
  // threshold: the generator evaluates to the zero of the field or it does not.
  for (const auto& g : generators) {
    EXPECT_TRUE(varietas::coefficient_traits<rational>::is_zero(g.evaluate(known)))
        << "the known configuration must satisfy every generator";
  }

  // And solving recovers both branches, one of them the configuration the
  // target was built from. The other is its elbow reflection, which for a
  // planar two-link arm is the mirror of the wrist in the base-to-tool line, so
  // the two joint angles sum to the same tool direction with q2 reversed.
  const auto solved = varietas::solve_zero_dimensional(solutions.basis());
  ASSERT_TRUE(solved.ok()) << to_string(solved.status);

  const auto real = solved.real_points();
  ASSERT_EQ(real.size(), 2u);

  const double known0 = 0.5;
  const double known1 = 1.0 / 3.0;
  bool found_known = false;
  for (const auto& point : real) {
    if (std::fabs(point[0] - known0) < 1e-9 && std::fabs(point[1] - known1) < 1e-9) {
      found_known = true;
    }
  }
  EXPECT_TRUE(found_known) << "the configuration the pose was built from must be found";

  // The two branches are distinct, and both put the tool at the target.
  EXPECT_GT(std::fabs(real[0][1] - real[1][1]), 1e-6);
  const chain<double> numeric = varietas::chain_cast<double>(robot);
  for (const auto& point : real) {
    const auto pose = varietas::forward_kinematics(
        numeric, {varietas::angle_from_variable(point[0]),
                  varietas::angle_from_variable(point[1])});
    for (std::size_t i = 0; i < 3; ++i) {
      EXPECT_NEAR(pose.translation()[i],
                  varietas::coefficient_traits<rational>::to_double(target[i]), 1e-9);
    }
  }
}

// A point off the plane of the arm is unreachable in the strong sense: the tool
// has z = 0 identically, so one residual reads 0 = 1 and the variety is empty
// over the algebraic closure. That is the case the weak Nullstellensatz detects,
// and the ideal is the whole ring.
TEST(rationalize, a_point_off_the_plane_gives_the_unit_ideal) {
  const chain<rational> robot = planar_two_link<rational>();
  const vector3<rational> target(varietas::make_rational(0), varietas::make_rational(0),
                                 varietas::make_rational(1));

  const auto generators =
      varietas::position_ideal_generators<2, grevlex>(robot, target);

  EXPECT_TRUE(varietas::is_unit_ideal(generators));
}

// A point beyond the reach of the arm but in its plane is a different matter,
// and the distinction is the one the library is built to respect. The arm
// cannot reach five units away with two unit links, but the ideal is not the
// unit ideal: over the algebraic closure the equations are perfectly solvable,
// with cos q2 = 23/2 and an angle that is not real. Unreachability is a
// statement about the real points of the variety, not about the ideal, and so
// it appears where the solver separates the real solutions from the complex
// ones and nowhere earlier.
TEST(rationalize, a_point_beyond_reach_has_solutions_but_no_real_ones) {
  const chain<rational> robot = planar_two_link<rational>();
  const vector3<rational> target(varietas::make_rational(5), varietas::make_rational(0),
                                 varietas::make_rational(0));

  const auto generators =
      varietas::position_ideal_generators<2, grevlex>(robot, target);
  const varietas::ideal<rational, 2, grevlex> solutions(generators);

  EXPECT_FALSE(varietas::is_unit_ideal(generators));
  ASSERT_TRUE(solutions.is_zero_dimensional());
  EXPECT_EQ(solutions.quotient().monomials.size(), 2u) << "still two branches";

  const auto solved = varietas::solve_zero_dimensional(solutions.basis());
  ASSERT_TRUE(solved.ok()) << to_string(solved.status);
  EXPECT_TRUE(solved.real_points().empty())
      << "the arm cannot reach it, and that is a fact about the real points";
}

// Saturation is not a tidying step here. Without it the ideal of a reachable
// point carries the components the substitution invented, and the count that
// certifies completeness is a count of the wrong thing.
TEST(rationalize, saturation_changes_the_answer_on_a_real_arm) {
  const chain<rational> robot = planar_two_link<rational>();
  const auto map = varietas::rational_forward_kinematics<2, grevlex>(robot);

  const std::array<rational, 2> known = {varietas::make_rational(1, 2),
                                         varietas::make_rational(1, 3)};
  const vector3<rational> target = map.evaluate(known).translation();

  const auto residuals = varietas::position_residuals(map, target);
  const varietas::ideal<rational, 2, grevlex> unsaturated(residuals);
  const varietas::ideal<rational, 2, grevlex> saturated(
      varietas::position_ideal_generators<2, grevlex>(robot, target));

  ASSERT_TRUE(saturated.is_zero_dimensional());
  EXPECT_EQ(saturated.quotient().monomials.size(), 2u);

  // The unsaturated ideal is strictly larger as a variety: it still contains
  // the two configurations, and something else besides.
  EXPECT_NE(unsaturated.basis(), saturated.basis());
  if (unsaturated.is_zero_dimensional()) {
    EXPECT_GT(unsaturated.quotient().monomials.size(),
              saturated.quotient().monomials.size());
  }
}

}  // namespace
