#include <array>

#include <gtest/gtest.h>

#include "varietas/core/monomial.hpp"
#include "varietas/core/order/grevlex.hpp"
#include "varietas/core/order/grlex.hpp"
#include "varietas/core/order/lex.hpp"
#include "varietas/core/polynomial.hpp"

namespace {

using varietas::grevlex;
using varietas::grlex;
using varietas::lex;

using mon3 = varietas::monomial<3>;

// Written in terms of the monomial's own exponent type rather than a fixed
// width, so that widening it does not have to be followed through the tests.
using exponent = mon3::exponent_type;

mon3 make(exponent a, exponent b, exponent c) {
  return mon3(std::array<exponent, 3>{a, b, c});
}

// Exponents past what a byte holds.
//
// The exponent used to be eight bits, and a product that carried one past 255
// wrapped silently to a different monomial — which is not a loud failure but a
// wrong answer, since every divisibility question asked afterwards is answered
// about the wrong thing. Parameter polynomials over a function field reach this
// range in practice; a degree of 254 in one variable was observed in an
// underdetermined system before this was widened.
TEST(Monomial, ExponentsPastAByteAreExact) {
  const mon3 a = make(200, 0, 0);
  const mon3 b = make(100, 0, 0);
  const mon3 product = a * b;

  EXPECT_EQ(product[0], 300u) << "the exponent wrapped instead of growing";
  EXPECT_EQ(product.degree(), 300u);
  EXPECT_TRUE(mon3::divides(a, product));
  EXPECT_EQ(mon3::divide(product, a), b);
}

TEST(Monomial, DegreeAndProduct) {
  const mon3 m = make(1, 2, 0);
  EXPECT_EQ(m.degree(), 3);
  const mon3 p = m * make(0, 1, 4);
  EXPECT_EQ(p.degree(), 8);
  EXPECT_EQ(p[1], 3);
  EXPECT_EQ(p[2], 4);
}

TEST(Monomial, Divisibility) {
  EXPECT_TRUE(mon3::divides(make(1, 1, 0), make(2, 1, 3)));
  EXPECT_FALSE(mon3::divides(make(1, 2, 0), make(2, 1, 3)));
  EXPECT_EQ(mon3::divide(make(2, 1, 3), make(1, 1, 0)), make(1, 0, 3));
  EXPECT_EQ(mon3::lcm(make(2, 0, 1), make(0, 3, 1)), make(2, 3, 1));
  EXPECT_EQ(mon3::gcd(make(2, 0, 1), make(1, 3, 4)), make(1, 0, 1));
}

TEST(Order, LexRanksFirstVariableHardest) {
  EXPECT_GT(lex::compare(make(1, 0, 0), make(0, 5, 5)), 0);
  EXPECT_EQ(lex::compare(make(1, 2, 3), make(1, 2, 3)), 0);
}

TEST(Order, GrlexUsesDegreeFirst) {
  EXPECT_GT(grlex::compare(make(0, 5, 5), make(1, 0, 0)), 0);
  EXPECT_GT(grlex::compare(make(1, 1, 1), make(1, 0, 2)), 0);
}

// The standard example distinguishing grlex from grevlex in three variables.
TEST(Order, GrevlexDiffersFromGrlex) {
  const mon3 a = make(1, 2, 0);
  const mon3 b = make(0, 0, 3);
  EXPECT_GT(grlex::compare(a, b), 0);
  EXPECT_GT(grevlex::compare(a, b), 0);

  const mon3 c = make(1, 0, 2);
  const mon3 d = make(0, 3, 0);
  EXPECT_GT(grlex::compare(c, d), 0);
  EXPECT_LT(grevlex::compare(c, d), 0);
}

using poly = varietas::polynomial<double, 3, grevlex>;

TEST(Polynomial, TermsAreSortedAndCombined) {
  const poly p({{make(0, 0, 1), 2.0}, {make(2, 0, 0), 3.0}, {make(0, 0, 1), 5.0}});
  ASSERT_EQ(p.size(), 2u);
  EXPECT_EQ(p.leading_monomial(), make(2, 0, 0));
  EXPECT_DOUBLE_EQ(p.leading_coefficient(), 3.0);
  EXPECT_DOUBLE_EQ(p.coefficient_of(make(0, 0, 1)), 7.0);
  EXPECT_DOUBLE_EQ(p.coefficient_of(make(1, 1, 1)), 0.0);
}

TEST(Polynomial, CancellationRemovesTerms) {
  const poly p({{make(1, 0, 0), 1.0}, {make(0, 1, 0), 2.0}});
  const poly q({{make(1, 0, 0), 1.0}, {make(0, 1, 0), -2.0}});
  const poly r = p - q;
  ASSERT_EQ(r.size(), 1u);
  EXPECT_EQ(r.leading_monomial(), make(0, 1, 0));
  EXPECT_DOUBLE_EQ(r.leading_coefficient(), 4.0);
  EXPECT_TRUE((p - p).is_zero());
}

TEST(Polynomial, ProductMatchesPointwiseEvaluation) {
  const poly p({{make(1, 0, 0), 2.0}, {make(0, 1, 0), -1.0}, {make(0, 0, 0), 3.0}});
  const poly q({{make(0, 0, 1), 4.0}, {make(0, 1, 0), 1.0}});
  const poly r = p * q;

  const std::array<double, 3> x{0.5, -1.25, 2.0};
  EXPECT_DOUBLE_EQ(r.evaluate(x), p.evaluate(x) * q.evaluate(x));
  EXPECT_EQ(r.degree(), 2);
}

TEST(Polynomial, MultiplyByMonomialShiftsExponents) {
  const poly p({{make(1, 0, 0), 2.0}, {make(0, 0, 0), 1.0}});
  const poly r = p.multiply_by_monomial(make(0, 2, 0), 3.0);
  EXPECT_EQ(r.leading_monomial(), make(1, 2, 0));
  EXPECT_DOUBLE_EQ(r.leading_coefficient(), 6.0);
  EXPECT_DOUBLE_EQ(r.coefficient_of(make(0, 2, 0)), 3.0);
}

TEST(Polynomial, PruneIsExplicit) {
  poly p({{make(1, 0, 0), 1.0}, {make(0, 1, 0), 1e-18}});
  EXPECT_EQ(p.size(), 2u);
  p.prune(1e-12);
  EXPECT_EQ(p.size(), 1u);
}

}  // namespace
