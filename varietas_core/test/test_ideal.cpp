#include <array>
#include <vector>

#include <gtest/gtest.h>

#include "varietas/core/ideal/buchberger.hpp"
#include "varietas/core/ideal/division.hpp"
#include "varietas/core/ideal/ideal.hpp"
#include "varietas/core/order/block.hpp"
#include "varietas/core/order/grevlex.hpp"
#include "varietas/core/order/grlex.hpp"
#include "varietas/core/order/lex.hpp"
#include "varietas/core/order/weight.hpp"
#include "varietas/core/polynomial.hpp"

namespace {

using varietas::grevlex;
using varietas::grlex;
using varietas::lex;

using mon3 = varietas::monomial<3>;
using poly_lex = varietas::polynomial<double, 3, lex>;
using poly_grlex = varietas::polynomial<double, 3, grlex>;

// Written in terms of the monomial's own exponent type rather than a fixed
// width, so that widening it does not have to be followed through the tests.
using exponent = mon3::exponent_type;

mon3 make(exponent a, exponent b, exponent c) {
  return mon3(std::array<exponent, 3>{a, b, c});
}

// x, y, z as polynomials under lex.
poly_lex x() { return poly_lex::variable(0); }
poly_lex y() { return poly_lex::variable(1); }
poly_lex z() { return poly_lex::variable(2); }
poly_lex constant(double c) { return poly_lex::constant(c); }

TEST(Order, BlockOrderRanksLeadingBlockFirst) {
  using elimination = varietas::block_order<1, lex, grevlex>;

  // Any monomial involving x dominates every monomial free of x.
  EXPECT_GT(elimination::compare(make(1, 0, 0), make(0, 9, 9)), 0);
  // Within the trailing block the tie break decides.
  EXPECT_GT(elimination::compare(make(1, 0, 2), make(1, 0, 1)), 0);
  EXPECT_EQ(elimination::compare(make(2, 1, 1), make(2, 1, 1)), 0);
}

TEST(Order, WeightOrderPrefersHeavyVariables) {
  struct weights {
    static constexpr std::int64_t weight(std::size_t i) noexcept { return i == 0 ? 5 : 1; }
  };
  using heavy_x = varietas::weight_order<weights, grevlex>;

  EXPECT_GT(heavy_x::compare(make(1, 0, 0), make(0, 4, 0)), 0);
  EXPECT_LT(heavy_x::compare(make(1, 0, 0), make(0, 6, 0)), 0);
  // Equal weight falls through to grevlex.
  EXPECT_EQ(heavy_x::compare(make(1, 0, 0), make(0, 5, 0)),
            grevlex::compare(make(1, 0, 0), make(0, 5, 0)));
}

TEST(Division, ReconstructsTheDividend) {
  const poly_lex f = x() * x() * y() + x() * y() * y() + y() * y();
  const std::vector<poly_lex> divisors{x() * y() - constant(1.0), y() * y() - constant(1.0)};

  const auto result = varietas::divide(f, divisors);

  poly_lex reconstructed = result.remainder;
  for (std::size_t i = 0; i < divisors.size(); ++i) {
    reconstructed += result.quotients[i] * divisors[i];
  }
  EXPECT_EQ(reconstructed, f);

  // No monomial of the remainder is divisible by a leading monomial.
  for (const auto& t : result.remainder.terms()) {
    for (const auto& g : divisors) {
      EXPECT_FALSE(mon3::divides(g.leading_monomial(), t.mon));
    }
  }
}

// The classical demonstration that the remainder of division is not an ideal
// membership test unless the divisors form a Gröbner basis.
TEST(Division, RemainderDependsOnDivisorOrderOutsideAGroebnerBasis) {
  const poly_lex f = x() * y() * y() - x();
  const poly_lex g1 = x() * y() + constant(1.0);
  const poly_lex g2 = y() * y() - constant(1.0);

  const auto forward = varietas::divide(f, std::vector<poly_lex>{g1, g2});
  const auto backward = varietas::divide(f, std::vector<poly_lex>{g2, g1});
  EXPECT_NE(forward.remainder, backward.remainder);

  // f does lie in the ideal, and the Gröbner basis says so.
  const varietas::ideal<double, 3, lex> i(std::vector<poly_lex>{g1, g2});
  EXPECT_TRUE(i.contains(f));
}

TEST(Buchberger, SPolynomialCancelsLeadingTerms) {
  const poly_grlex f({{make(2, 1, 0), 1.0}, {make(0, 2, 0), 1.0}});
  const poly_grlex g({{make(1, 2, 0), 3.0}, {make(0, 0, 1), -1.0}});

  const poly_grlex s = varietas::s_polynomial(f, g);
  EXPECT_LT(grlex::compare(s.leading_monomial(), mon3::lcm(f.leading_monomial(),
                                                           g.leading_monomial())),
            0);
}

// Cox, Little and O'Shea, Chapter 2, Section 7, Example 2: the reduced basis of
// this ideal under grlex is known in closed form.
TEST(Buchberger, ReducedBasisOfTheStandardExample) {
  const poly_grlex f({{make(3, 0, 0), 1.0}, {make(0, 2, 0), -2.0}});
  const poly_grlex g({{make(2, 1, 0), 1.0}, {make(0, 0, 0), -2.0}, {make(1, 0, 0), 1.0}});

  const auto basis = varietas::groebner_basis(std::vector<poly_grlex>{f, g});

  // Every element is monic and no leading monomial divides another.
  for (const auto& b : basis) {
    EXPECT_DOUBLE_EQ(b.leading_coefficient(), 1.0);
  }
  for (std::size_t i = 0; i < basis.size(); ++i) {
    for (std::size_t j = 0; j < basis.size(); ++j) {
      if (i != j) {
        EXPECT_FALSE(mon3::divides(basis[j].leading_monomial(), basis[i].leading_monomial()));
      }
    }
  }

  // Every S-polynomial reduces to zero, which is Buchberger's criterion.
  for (std::size_t i = 0; i < basis.size(); ++i) {
    for (std::size_t j = i + 1; j < basis.size(); ++j) {
      EXPECT_TRUE(varietas::normal_form(varietas::s_polynomial(basis[i], basis[j]), basis)
                      .is_zero());
    }
  }

  // The generators lie in the ideal they generate.
  EXPECT_TRUE(varietas::is_member(f, basis));
  EXPECT_TRUE(varietas::is_member(g, basis));
}

TEST(Buchberger, ReducedBasisIsIndependentOfTheGeneratingSet) {
  const poly_grlex f({{make(2, 0, 0), 1.0}, {make(0, 2, 0), 1.0}, {make(0, 0, 0), -1.0}});
  const poly_grlex g({{make(1, 0, 0), 1.0}, {make(0, 1, 0), -1.0}});

  const auto first = varietas::groebner_basis(std::vector<poly_grlex>{f, g});
  // The same ideal, presented by a different and redundant generating set.
  const auto second = varietas::groebner_basis(
      std::vector<poly_grlex>{g, f + g * 3.0, f, f * 2.0 - g});

  ASSERT_EQ(first.size(), second.size());
  for (std::size_t i = 0; i < first.size(); ++i) {
    ASSERT_EQ(first[i].size(), second[i].size());
    for (std::size_t k = 0; k < first[i].size(); ++k) {
      EXPECT_EQ(first[i].terms()[k].mon, second[i].terms()[k].mon);
      EXPECT_NEAR(first[i].terms()[k].coeff, second[i].terms()[k].coeff, 1e-12);
    }
  }
}

TEST(Buchberger, InconsistentSystemYieldsTheUnitIdeal) {
  // x^2 + y^2 = 1 together with x^2 + y^2 = 2 has no solution anywhere.
  const poly_grlex f({{make(2, 0, 0), 1.0}, {make(0, 2, 0), 1.0}, {make(0, 0, 0), -1.0}});
  const poly_grlex g({{make(2, 0, 0), 1.0}, {make(0, 2, 0), 1.0}, {make(0, 0, 0), -2.0}});

  const auto basis = varietas::groebner_basis(std::vector<poly_grlex>{f, g});
  ASSERT_EQ(basis.size(), 1u);
  EXPECT_TRUE(basis.front().leading_monomial().is_one());
  EXPECT_TRUE(varietas::is_unit_ideal(basis));
}

TEST(Buchberger, CriteriaDiscardPairsWithoutChangingTheBasis) {
  // Coprime leading monomials: both S-polynomials are certified in advance.
  const poly_grlex f = poly_grlex::variable(0, 2) - poly_grlex::constant(1.0);
  const poly_grlex g = poly_grlex::variable(1, 3) - poly_grlex::constant(1.0);

  varietas::buchberger_statistics statistics;
  const auto basis = varietas::groebner_basis(std::vector<poly_grlex>{f, g}, &statistics);

  EXPECT_EQ(statistics.pairs_discarded_coprime, 1u);
  EXPECT_EQ(statistics.pairs_reduced, 0u);
  EXPECT_EQ(basis.size(), 2u);
}

// The Elimination Theorem: under lex with x > y > z, the basis elements free of
// x generate the projection of the variety onto the remaining coordinates.
TEST(Elimination, LexBasisContainsTheEliminationIdeal) {
  const varietas::ideal<double, 3, lex> i(std::vector<poly_lex>{
      x() + y() + z() - constant(3.0),
      x() * y() + y() * z() + z() * x() - constant(3.0),
      x() * y() * z() - constant(1.0),
  });

  const auto eliminated = i.eliminate(1);
  ASSERT_FALSE(eliminated.empty());
  for (const auto& g : eliminated) {
    for (const auto& t : g.terms()) {
      EXPECT_EQ(t.mon[0], 0);
    }
    EXPECT_TRUE(i.contains(g));
  }

  // The system is the elementary symmetric system of (1, 1, 1), so the
  // projection is the single point y = z = 1 and every eliminant vanishes there.
  for (const auto& g : eliminated) {
    EXPECT_NEAR(g.evaluate(std::array<double, 3>{0.0, 1.0, 1.0}), 0.0, 1e-9);
  }
}

TEST(Ideal, BasisIsCachedAndInvalidatedByNewGenerators) {
  varietas::ideal<double, 3, grlex> i(std::vector<poly_grlex>{
      poly_grlex::variable(0, 2) - poly_grlex::constant(1.0)});

  EXPECT_EQ(i.basis().size(), 1u);
  EXPECT_FALSE(i.is_unit());

  i.add_generator(poly_grlex::variable(0) - poly_grlex::constant(1.0));
  ASSERT_EQ(i.basis().size(), 1u);
  EXPECT_EQ(i.basis().front().leading_monomial(), make(1, 0, 0));
  EXPECT_FALSE(i.is_unit());
}

}  // namespace
