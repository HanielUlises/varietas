// The bridge: a chain in, a parametric_solution out.
//
// The property that matters is specialisation. A basis computed over Q(x, y)
// claims to describe the arm at every pose, and the way to test that claim is
// to substitute a pose into it and compare against the basis computed over Q
// with that pose substituted first. Those are two genuinely different
// computations — one Buchberger run over a function field, one over the
// rationals — and they are required to produce the same quotient algebra and
// the same multiplication operators. Everything else here is a statement about
// what the bridge refuses.

#include <array>
#include <cstddef>
#include <vector>

#include <gtest/gtest.h>

#include "varietas/codegen/rational.hpp"
#include "varietas/core/order/grevlex.hpp"
#include "varietas/core/quotient/action_matrix.hpp"
#include "varietas/core/quotient/quotient_basis.hpp"
#include "varietas/ik/parametric_ik.hpp"
#include "varietas/kinematics/rationalize.hpp"

#include "arms.hpp"

namespace {

using varietas::grevlex;
using varietas::rational;
using varietas::ik::parametric_ik_status;
using varietas::ik::parametric_position_ik;

// --- the planar two-link arm solves over Q(x, y) ---------------------------

TEST(parametric_ik, planar_two_link_is_zero_dimensional_with_two_branches) {
  const auto result =
      parametric_position_ik<2, 2>(varietas_test::planar_two_link(), {0, 1});

  ASSERT_TRUE(result.ok()) << varietas::ik::to_string(result.status);
  // Two elbow configurations, and no spurious branch from the half-angle
  // denominators, which is what the saturation is there to guarantee.
  EXPECT_EQ(result.branches, 2u);
  EXPECT_EQ(result.solution.dimension(), 2u);
  EXPECT_TRUE(result.solution.is_well_formed());
}

TEST(parametric_ik, records_the_order_and_the_names) {
  const auto result =
      parametric_position_ik<2, 2>(varietas_test::planar_two_link(), {0, 1});
  ASSERT_TRUE(result.ok());

  EXPECT_EQ(result.solution.order, varietas::order_id::grevlex);
  // The unknowns are the half-angle variables, not the joint angles, and the
  // names say so.
  ASSERT_EQ(result.solution.unknown_names.size(), 2u);
  EXPECT_EQ(result.solution.unknown_names[0], "t_q1");
  EXPECT_EQ(result.solution.unknown_names[1], "t_q2");

  ASSERT_EQ(result.solution.parameter_names.size(), 2u);
  EXPECT_EQ(result.solution.parameter_names[0], "x");
  EXPECT_EQ(result.solution.parameter_names[1], "y");
}

TEST(parametric_ik, one_is_a_standard_monomial_and_one_index_points_at_it) {
  const auto result =
      parametric_position_ik<2, 2>(varietas_test::planar_two_link(), {0, 1});
  ASSERT_TRUE(result.ok());

  ASSERT_LT(result.solution.one_index, result.solution.dimension());
  EXPECT_EQ(result.solution.quotient.monomials[result.solution.one_index].degree(), 0u);
}

TEST(parametric_ik, the_parameter_order_follows_the_coordinates_argument) {
  // Asking for (y, x) has to produce a solver whose first pose component is y.
  const auto result =
      parametric_position_ik<2, 2>(varietas_test::planar_two_link(), {1, 0});
  ASSERT_TRUE(result.ok()) << varietas::ik::to_string(result.status);

  ASSERT_EQ(result.solution.parameter_names.size(), 2u);
  EXPECT_EQ(result.solution.parameter_names[0], "y");
  EXPECT_EQ(result.solution.parameter_names[1], "x");
  EXPECT_EQ(result.branches, 2u);
}

// --- specialisation: the parametric basis agrees with the basis at a pose ---

// The ideal of the same arm at a fixed rational pose, computed over Q.
std::vector<varietas::polynomial<rational, 2, grevlex>> exact_basis_at(const rational& x,
                                                                      const rational& y) {
  using poly = varietas::polynomial<rational, 2, grevlex>;

  const auto robot = varietas_test::planar_two_link();
  const auto map = varietas::rational_forward_kinematics<2, grevlex>(robot);
  const poly denominator = map.denominator();

  std::vector<poly> residuals;
  residuals.push_back(map.translation(0) - denominator * poly::constant(x));
  residuals.push_back(map.translation(1) - denominator * poly::constant(y));
  return varietas::kinematic_ideal_generators<2, grevlex>(robot, residuals);
}

TEST(parametric_ik, specialises_to_the_basis_computed_at_that_pose) {
  // A pose in the interior of the annulus and off every symmetry line, so that
  // nothing about the agreement is an accident of a special position.
  const rational x(3, 4);
  const rational y(5, 7);

  const auto parametric =
      parametric_position_ik<2, 2>(varietas_test::planar_two_link(), {0, 1});
  ASSERT_TRUE(parametric.ok()) << varietas::ik::to_string(parametric.status);

  const auto exact = exact_basis_at(x, y);
  const auto exact_quotient = varietas::standard_monomials(exact);
  ASSERT_TRUE(exact_quotient.is_zero_dimensional);

  // The parametric basis specialises well at this pose exactly when the two
  // quotient algebras have the same standard monomials. A drop in dimension
  // here would mean the pose lies on the locus the one basis fails to describe.
  const auto& quotient = parametric.solution.quotient;
  ASSERT_EQ(quotient.dimension(), exact_quotient.dimension());
  for (const auto& m : quotient.monomials) {
    ASSERT_LT(exact_quotient.index_of(m), exact_quotient.dimension())
        << "a standard monomial over Q(x, y) is not one over Q at this pose";
  }

  const std::array<rational, 2> pose{x, y};
  for (std::size_t variable = 0; variable < 2; ++variable) {
    const auto expected = varietas::variable_action_matrix(variable, exact, exact_quotient);
    const auto& actual = parametric.solution.action[variable];

    for (std::size_t i = 0; i < quotient.dimension(); ++i) {
      for (std::size_t j = 0; j < quotient.dimension(); ++j) {
        // Indexed through the monomials rather than positionally, so that the
        // comparison does not quietly depend on the two enumerations agreeing.
        const std::size_t row = exact_quotient.index_of(quotient.monomials[i]);
        const std::size_t column = exact_quotient.index_of(quotient.monomials[j]);
        const rational entry = actual(i, j).evaluate(pose);
        EXPECT_NEAR(entry.get_d(), expected(row, column), 1e-12)
            << "M_" << variable << "(" << i << ", " << j << ")";
      }
    }
  }
}

TEST(parametric_ik, variable_coordinates_specialise_too) {
  const rational x(3, 4);
  const rational y(5, 7);

  const auto parametric =
      parametric_position_ik<2, 2>(varietas_test::planar_two_link(), {0, 1});
  ASSERT_TRUE(parametric.ok());

  const auto exact = exact_basis_at(x, y);
  const auto exact_quotient = varietas::standard_monomials(exact);
  const auto& quotient = parametric.solution.quotient;

  const std::array<rational, 2> pose{x, y};
  for (std::size_t variable = 0; variable < 2; ++variable) {
    // The normal form of x_i over Q, for comparison with the specialised one.
    const auto reduced = varietas::normal_form(
        varietas::polynomial<rational, 2, grevlex>::variable(variable), exact);

    std::vector<double> expected(exact_quotient.dimension(), 0.0);
    for (const auto& t : reduced.terms()) {
      expected[exact_quotient.index_of(t.mon)] = t.coeff.get_d();
    }

    for (std::size_t k = 0; k < quotient.dimension(); ++k) {
      const std::size_t index = exact_quotient.index_of(quotient.monomials[k]);
      EXPECT_NEAR(parametric.solution.variable_coordinates[variable][k].evaluate(pose).get_d(),
                  expected[index], 1e-12);
    }
  }
}

// --- what the bridge refuses ------------------------------------------------

TEST(parametric_ik, refuses_to_drop_a_coordinate_the_arm_actually_moves) {
  // The second axis is x, so the tool leaves the plane and the z equation is
  // not vacuous. Dropping it would silently answer about a larger variety.
  const auto result =
      parametric_position_ik<2, 2>(varietas_test::out_of_plane_two_link(), {0, 1});

  EXPECT_EQ(result.status, parametric_ik_status::dropped_coordinate_is_not_identically_zero);
  EXPECT_EQ(result.offending_coordinate, 2u);
}

TEST(parametric_ik, refuses_more_unknowns_than_coordinates_without_computing) {
  // Three coplanar joints against two equations: a curve of configurations for
  // each reachable point, and no action matrix to emit. Krull settles this by
  // counting, so it must come back promptly rather than after a Buchberger run
  // over Q(x, y) that has no useful answer to reach.
  const auto result =
      parametric_position_ik<3, 2>(varietas_test::planar_three_link(), {0, 1});

  EXPECT_EQ(result.status, parametric_ik_status::underdetermined);
}

TEST(parametric_ik, refuses_more_coordinates_than_unknowns_without_computing) {
  // A two-joint planar arm asked for a point in space. The tool sweeps a
  // two-dimensional region, a general point of three-space is not in it, and
  // the ideal is the unit ideal — which the solve would discover, expensively,
  // over Q(x, y, z). Counting gets there first.
  const auto result =
      parametric_position_ik<2, 3>(varietas_test::planar_two_link(), {0, 1, 2});

  EXPECT_EQ(result.status, parametric_ik_status::overdetermined);
  EXPECT_EQ(result.branches, 0u);
}

TEST(parametric_ik, refuses_a_positive_dimensional_system_that_counting_cannot_catch) {
  // Two joints about the same axis through the same point, so only their sum
  // moves the tool. The generator count says nothing here — two equations in two
  // unknowns — and only the quotient dimension reveals that the fibre over a
  // reachable point is a curve.
  const auto result =
      parametric_position_ik<2, 2>(varietas_test::coincident_two_link(), {0, 1});

  EXPECT_EQ(result.status, parametric_ik_status::not_zero_dimensional)
      << varietas::ik::to_string(result.status);
}

TEST(parametric_ik, refuses_repeated_or_out_of_range_coordinates) {
  const auto repeated =
      parametric_position_ik<2, 2>(varietas_test::planar_two_link(), {0, 0});
  EXPECT_EQ(repeated.status, parametric_ik_status::bad_coordinates);

  const auto out_of_range =
      parametric_position_ik<2, 2>(varietas_test::planar_two_link(), {0, 3});
  EXPECT_EQ(out_of_range.status, parametric_ik_status::bad_coordinates);
  EXPECT_EQ(out_of_range.offending_coordinate, 3u);
}

TEST(parametric_ik, refuses_a_chain_with_the_wrong_number_of_joints) {
  // Instantiated for three unknowns and three coordinates, handed an arm with
  // two joints. Stated this way the count check above does not fire first, so
  // it is really the degrees of freedom being checked.
  const auto result =
      parametric_position_ik<3, 3>(varietas_test::planar_two_link(), {0, 1, 2});
  EXPECT_EQ(result.status, parametric_ik_status::wrong_degrees_of_freedom);
}

}  // namespace
