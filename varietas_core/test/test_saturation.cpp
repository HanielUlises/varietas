#include <array>
#include <cstdint>
#include <vector>

#include <gtest/gtest.h>

#include "varietas/core/embed.hpp"
#include "varietas/core/ideal/ideal.hpp"
#include "varietas/core/ideal/saturation.hpp"
#include "varietas/core/order/grevlex.hpp"
#include "varietas/core/order/lex.hpp"
#include "varietas/core/polynomial.hpp"

namespace {

using varietas::grevlex;
using varietas::lex;
using varietas::monomial;

template <std::size_t N, class Order = grevlex>
using poly = varietas::polynomial<double, N, Order>;

template <std::size_t N, class Order = grevlex>
poly<N, Order> var(std::size_t i, std::uint8_t power = 1) {
  return poly<N, Order>::variable(i, power);
}

template <std::size_t N, class Order = grevlex>
poly<N, Order> constant(double c) {
  return poly<N, Order>::constant(c);
}

// Membership against a Gröbner basis, which is what a saturation is usually
// asked for: not the shape of the generators but which polynomials they contain.
template <std::size_t N, class Order>
bool contains(const std::vector<poly<N, Order>>& basis, const poly<N, Order>& f) {
  return varietas::is_member(f, basis);
}

// --- embedding -------------------------------------------------------------

TEST(embed, moves_a_polynomial_into_a_larger_ring_and_back) {
  // x^2 y + 3 in two variables, read as a polynomial in three with the
  // variables shifted up by one, then read back.
  const poly<2> f = var<2>(0, 2) * var<2>(1) + constant<2>(3.0);

  const auto lifted = varietas::embed_polynomial<lex, 3>(f, 1);
  EXPECT_EQ(lifted.size(), 2u);
  for (const auto& t : lifted.terms()) {
    EXPECT_EQ(t.mon[0], 0) << "the adjoined variable must not appear";
  }

  const auto returned = varietas::project_polynomial<grevlex, 2>(lifted, 1);
  EXPECT_EQ(returned, f);
}

TEST(embed, re_sorts_terms_under_the_target_order) {
  // Under grevlex the terms of x^2 + y^3 lead with y^3; under lex with x^2.
  // The embedding must sort under the target order rather than copy the list.
  const poly<2> f = var<2>(0, 2) + var<2>(1, 3);
  EXPECT_EQ(f.leading_monomial(), monomial<2>::variable(1, 3));

  const auto lifted = varietas::embed_polynomial<lex, 2>(f, 0);
  EXPECT_EQ(lifted.leading_monomial(), monomial<2>::variable(0, 2));
}

// --- saturation on textbook ideals -----------------------------------------

TEST(saturation, removes_the_component_lying_in_the_hypersurface) {
  // I = (xz, yz) = (z) ∩ (x, y): the plane z = 0 together with the line x = y = 0.
  // Saturating by z discards the plane and leaves the line.
  const poly<3> x = var<3>(0);
  const poly<3> y = var<3>(1);
  const poly<3> z = var<3>(2);

  const auto saturated = varietas::saturate<double, 3, grevlex>({x * z, y * z}, z);

  EXPECT_TRUE(contains(saturated, x));
  EXPECT_TRUE(contains(saturated, y));
  EXPECT_FALSE(contains(saturated, z)) << "the line is not contained in z = 0";

  // (x, y) exactly, not merely something containing it.
  const varietas::ideal<double, 3, grevlex> expected({x, y});
  EXPECT_EQ(saturated, expected.basis());
}

TEST(saturation, a_purely_embedded_component_saturates_to_the_unit_ideal) {
  // V(x^2) is the line x = 0 with multiplicity; removing x = 0 leaves nothing,
  // and the ideal of the empty set is the whole ring.
  const poly<2> x = var<2>(0);

  const auto saturated = varietas::saturate<double, 2, grevlex>({x * x}, x);

  ASSERT_EQ(saturated.size(), 1u);
  EXPECT_TRUE(saturated.front().leading_monomial().is_one());
  EXPECT_TRUE(varietas::is_unit_ideal(saturated));
}

TEST(saturation, leaves_an_ideal_alone_when_the_multiplier_misses_it) {
  // V(x) is not contained in y = 0, so removing that hypersurface and taking
  // the closure gives V(x) back.
  const poly<2> x = var<2>(0);
  const poly<2> y = var<2>(1);

  const auto saturated = varietas::saturate<double, 2, grevlex>({x}, y);

  const varietas::ideal<double, 2, grevlex> expected({x});
  EXPECT_EQ(saturated, expected.basis());
}

TEST(saturation, by_a_nonzero_constant_changes_nothing) {
  const poly<2> x = var<2>(0);
  const poly<2> y = var<2>(1);
  const std::vector<poly<2>> generators = {x * x - y, x * y};

  const auto saturated =
      varietas::saturate<double, 2, grevlex>(generators, constant<2>(3.0));

  const varietas::ideal<double, 2, grevlex> expected(generators);
  EXPECT_EQ(saturated, expected.basis());
}

TEST(saturation, is_idempotent) {
  const poly<3> x = var<3>(0);
  const poly<3> y = var<3>(1);
  const poly<3> z = var<3>(2);

  const auto once = varietas::saturate<double, 3, grevlex>({x * z, y * z}, z);
  const auto twice = varietas::saturate<double, 3, grevlex>(once, z);

  EXPECT_EQ(once, twice);
}

TEST(saturation, of_the_zero_ideal_is_the_zero_ideal) {
  const auto saturated = varietas::saturate<double, 2, grevlex>({}, var<2>(0));
  EXPECT_TRUE(saturated.empty());
}

TEST(saturation, returns_a_basis_under_the_requested_order) {
  // The trailing block of the elimination order is the caller's order, so what
  // comes back is a Gröbner basis under that order and needs no second run.
  // Under lex the basis of (xz, yz) : z^∞ = (x, y) leads with x and y.
  const poly<3, lex> x = var<3, lex>(0);
  const poly<3, lex> y = var<3, lex>(1);
  const poly<3, lex> z = var<3, lex>(2);

  const auto saturated = varietas::saturate<double, 3, lex>({x * z, y * z}, z);

  const varietas::ideal<double, 3, lex> expected({x, y});
  EXPECT_EQ(saturated, expected.basis());
  EXPECT_TRUE(contains(saturated, x * y + x));
}

// --- the case the header exists for ----------------------------------------

// The half-angle substitution clears denominators 1 + t_i^2, and doing so
// attaches to the variety the components on which those factors vanish: the
// points t_i = ±i, images of the configurations q_i = π that the substitution
// sent to infinity. Over the reals they are invisible, which is exactly why
// they have to be removed symbolically: they are counted by dim_k A all the
// same, and the count is the completeness certificate.
// The mild form: the spurious component is a finite set of extra points, and
// the effect is that dim_k A counts solutions the arm does not have.
TEST(saturation, stops_the_quotient_counting_configurations_that_do_not_exist) {
  const poly<2> t0 = var<2>(0);
  const poly<2> t1 = var<2>(1);
  const poly<2> denominator = constant<2>(1.0) + t0 * t0;

  const std::vector<poly<2>> generators = {denominator * (t0 - constant<2>(1.0)),
                                           t1 - constant<2>(1.0)};

  // One real configuration, t0 = 1, and the two points t0 = ±i the substitution
  // invented. All three are counted, and the count is what certifies that no
  // branch of the inverse kinematics was missed.
  const varietas::ideal<double, 2, grevlex> unsaturated(generators);
  ASSERT_TRUE(unsaturated.is_zero_dimensional());
  EXPECT_EQ(unsaturated.quotient().monomials.size(), 3u);

  const auto saturated = varietas::saturate<double, 2, grevlex>(generators, denominator);

  const varietas::ideal<double, 2, grevlex> cleaned(saturated);
  EXPECT_TRUE(cleaned.is_zero_dimensional());
  EXPECT_EQ(cleaned.quotient().monomials.size(), 1u);
  EXPECT_TRUE(contains(saturated, t0 - constant<2>(1.0)));
  EXPECT_TRUE(contains(saturated, t1 - constant<2>(1.0)));
}

// The severe form, and the one that actually arises: the denominator divides
// more than one generator, so the spurious locus is not a few extra points but
// a positive-dimensional component, here the line t0 = ±i with t1 free. The
// Finiteness Theorem then returns no verdict at all, and the solver refuses the
// system rather than the system being solved wrongly. Saturating is what makes
// the ideal zero-dimensional in the first place.
TEST(saturation, restores_the_finiteness_the_substitution_destroyed) {
  const poly<2> t0 = var<2>(0);
  const poly<2> t1 = var<2>(1);
  const poly<2> denominator = constant<2>(1.0) + t0 * t0;

  const std::vector<poly<2>> generators = {denominator * (t0 - constant<2>(1.0)),
                                           denominator * (t1 - constant<2>(1.0))};

  const varietas::ideal<double, 2, grevlex> unsaturated(generators);
  EXPECT_FALSE(unsaturated.is_zero_dimensional())
      << "1 + t0^2 = 0 leaves t1 free, so the variety contains a line";

  const auto saturated = varietas::saturate<double, 2, grevlex>(generators, denominator);

  const varietas::ideal<double, 2, grevlex> expected(
      {t0 - constant<2>(1.0), t1 - constant<2>(1.0)});
  EXPECT_EQ(saturated, expected.basis());

  const varietas::ideal<double, 2, grevlex> cleaned(saturated);
  EXPECT_TRUE(cleaned.is_zero_dimensional());
  EXPECT_EQ(cleaned.quotient().monomials.size(), 1u);
}

TEST(saturation, by_a_product_removes_every_factors_component) {
  // What the kinematic ideal calls for: one denominator per revolute joint,
  // and the product removes the components of all of them at once.
  const poly<2> t0 = var<2>(0);
  const poly<2> t1 = var<2>(1);
  const poly<2> d0 = varietas::half_angle_denominator<double, 2, grevlex>(0);
  const poly<2> d1 = varietas::half_angle_denominator<double, 2, grevlex>(1);

  EXPECT_EQ(d0, constant<2>(1.0) + t0 * t0);
  EXPECT_EQ(d1, constant<2>(1.0) + t1 * t1);

  const std::vector<poly<2>> generators = {d0 * d1 * (t0 - constant<2>(2.0)),
                                           d0 * d1 * (t1 - constant<2>(3.0))};

  const auto saturated =
      varietas::saturate_by_product<double, 2, grevlex>(generators, {d0, d1});

  const varietas::ideal<double, 2, grevlex> expected(
      {t0 - constant<2>(2.0), t1 - constant<2>(3.0)});
  EXPECT_EQ(saturated, expected.basis());
}

}  // namespace
