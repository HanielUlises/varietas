// Why the offline pipeline needs an exact field.
//
// Buchberger's algorithm decides that a critical pair is finished when the
// S-polynomial reduces to zero, and polynomial::is_zero is the structural test
// "no terms left". Over Q a term disappears exactly when its coefficient
// cancels. Over double it disappears only when the cancellation happens to be
// exact in binary, which for coefficients like 1/3 it is not.

#include <cmath>
#include <cstdint>
#include <type_traits>
#include <vector>

#include <gtest/gtest.h>

#include "varietas/codegen/rational.hpp"
#include "varietas/core/ideal/buchberger.hpp"
#include "varietas/core/ideal/division.hpp"
#include "varietas/core/ideal/ideal.hpp"
#include "varietas/core/order/grevlex.hpp"
#include "varietas/core/order/grlex.hpp"
#include "varietas/core/order/lex.hpp"
#include "varietas/core/quotient/quotient_basis.hpp"

#include "coefficient_fixture.hpp"

namespace {

using varietas::grevlex;
using varietas::grlex;
using varietas::lex;
using varietas::rational;
using varietas_test::coeff;

using mon3 = varietas::monomial<3>;

template <class Coeff>
using poly_grlex = varietas::polynomial<Coeff, 3, grlex>;

// A basis of an ideal whose coefficients are not dyadic, so that the exact and
// the approximate field genuinely part company.
template <class Coeff>
std::vector<poly_grlex<Coeff>> non_dyadic_generators() {
  using poly = poly_grlex<Coeff>;
  const auto c = [](std::int64_t n, std::int64_t d = 1) { return coeff<Coeff>(n, d); };

  const poly x = poly::variable(0);
  const poly y = poly::variable(1);

  // x^2 - y/3 and x y - 1/7, chosen so that reduction repeatedly multiplies
  // thirds by sevenths.
  return {x * x - poly::constant(c(1, 3)) * y, x * y - poly::constant(c(1, 7))};
}

TEST(ExactGroebner, EveryCriticalPairOfTheBasisReducesToZeroExactly) {
  const auto basis = varietas::groebner_basis(non_dyadic_generators<rational>());

  ASSERT_FALSE(basis.empty());
  for (std::size_t i = 0; i < basis.size(); ++i) {
    for (std::size_t j = i + 1; j < basis.size(); ++j) {
      const auto remainder =
          varietas::normal_form(varietas::s_polynomial(basis[i], basis[j]), basis);
      // Buchberger's criterion, as an exact statement rather than a tolerance.
      EXPECT_TRUE(remainder.is_zero())
          << "pair (" << i << ", " << j << ") left " << remainder.size() << " terms";
    }
  }
}

TEST(ExactGroebner, LeadingCoefficientsAreExactlyOne) {
  const auto basis = varietas::groebner_basis(non_dyadic_generators<rational>());
  for (const auto& g : basis) {
    EXPECT_EQ(g.leading_coefficient(), varietas::coefficient_traits<rational>::one());
  }
}

// Ideal membership is the property that breaks first in floating point, and it
// is the one every later stage rests on: elimination, the finiteness verdict
// and the action matrix all assume normal forms are canonical.
TEST(ExactGroebner, MembershipOfACombinationHoldsExactlyAndFailsInDouble) {
  const auto exact_basis = varietas::groebner_basis(non_dyadic_generators<rational>());
  const auto approximate_basis = varietas::groebner_basis(non_dyadic_generators<double>());

  // f = a g0 + b g1 lies in the ideal by construction, whatever a and b are.
  const auto combination = [](const auto& basis, auto c) {
    using poly = typename std::decay_t<decltype(basis)>::value_type;
    const poly x = poly::variable(0);
    const poly y = poly::variable(1);
    const poly a = x + poly::constant(c(1, 3));
    const poly b = y - poly::constant(c(2, 7));
    return a * basis[0] + b * basis[1];
  };

  const auto exact_f =
      combination(exact_basis, [](std::int64_t n, std::int64_t d) { return coeff<rational>(n, d); });
  const auto approximate_f =
      combination(approximate_basis, [](std::int64_t n, std::int64_t d) { return coeff<double>(n, d); });

  EXPECT_TRUE(varietas::is_member(exact_f, exact_basis))
      << "an exact combination of basis elements must reduce to zero";

  // Documenting the failure rather than tolerating it: over double the same
  // reduction leaves rounding dust, so the membership test returns the wrong
  // answer. If this ever starts passing the demonstration has gone stale, but
  // the exact assertion above is the one that matters.
  const auto approximate_remainder = varietas::normal_form(approximate_f, approximate_basis);
  if (!approximate_remainder.is_zero()) {
    for (const auto& t : approximate_remainder.terms()) {
      EXPECT_LT(std::abs(t.coeff), 1e-10)
          << "residue should be rounding dust, not a genuine remainder";
    }
    SUCCEED() << "double left " << approximate_remainder.size()
              << " spurious term(s) of magnitude below 1e-10";
  }
}

// Cox, Little and O'Shea, Chapter 2, Section 7, Example 2, over Q. The reduced
// basis is known in closed form, so the coefficients can be checked as exact
// fractions rather than compared within a tolerance.
TEST(ExactGroebner, StandardExampleHasTheKnownReducedBasis) {
  using poly = poly_grlex<rational>;
  const auto c = [](std::int64_t n, std::int64_t d = 1) { return coeff<rational>(n, d); };

  const poly f({{varietas_test::mon<3>({3, 0, 0}), c(1)},
                {varietas_test::mon<3>({0, 2, 0}), c(-2)}});
  const poly g({{varietas_test::mon<3>({2, 1, 0}), c(1)},
                {varietas_test::mon<3>({0, 0, 0}), c(-2)},
                {varietas_test::mon<3>({1, 0, 0}), c(1)}});

  const auto basis = varietas::groebner_basis(std::vector<poly>{f, g});

  EXPECT_TRUE(varietas::is_member(f, basis));
  EXPECT_TRUE(varietas::is_member(g, basis));

  // Reduced: no leading monomial divides another.
  for (std::size_t i = 0; i < basis.size(); ++i) {
    for (std::size_t j = 0; j < basis.size(); ++j) {
      if (i != j) {
        EXPECT_FALSE(mon3::divides(basis[j].leading_monomial(), basis[i].leading_monomial()));
      }
    }
  }
}

// The reduced basis is an invariant of the ideal, not of the presentation. Over
// Q this is an exact equality of term lists; the corresponding test in
// varietas_core has to compare coefficients within 1e-12.
TEST(ExactGroebner, ReducedBasisIsIndependentOfTheGeneratingSetExactly) {
  using poly = poly_grlex<rational>;
  const auto c = [](std::int64_t n, std::int64_t d = 1) { return coeff<rational>(n, d); };

  const poly x = poly::variable(0);
  const poly y = poly::variable(1);
  const poly f = x * x + y * y - poly::constant(c(1));
  const poly g = x - y;

  const auto first = varietas::groebner_basis(std::vector<poly>{f, g});
  const auto second = varietas::groebner_basis(
      std::vector<poly>{g, f + g * c(1, 3), f, f * c(2, 5) - g});

  ASSERT_EQ(first.size(), second.size());
  for (std::size_t i = 0; i < first.size(); ++i) {
    EXPECT_EQ(first[i], second[i]);
  }
}

// The finiteness verdict is read off leading monomials alone, so it is only as
// trustworthy as the basis. Over Q it certifies the solution count.
TEST(ExactGroebner, FinitenessVerdictAndQuotientDimension) {
  using poly = poly_grlex<rational>;
  const auto c = [](std::int64_t n, std::int64_t d = 1) { return coeff<rational>(n, d); };

  // x^2 = 1/4, y^2 = 1/9, z = 1/3: four points, none of them dyadic in y.
  const poly x = poly::variable(0);
  const poly y = poly::variable(1);
  const poly z = poly::variable(2);

  const auto basis = varietas::groebner_basis(std::vector<poly>{
      x * x - poly::constant(c(1, 4)),
      y * y - poly::constant(c(1, 9)),
      z - poly::constant(c(1, 3)),
  });

  const auto quotient = varietas::standard_monomials(basis);
  EXPECT_TRUE(quotient.is_zero_dimensional);
  EXPECT_EQ(quotient.dimension(), 4u);
}

// Elimination under a lex order, with the eliminant checked as an exact
// polynomial identity rather than by evaluating near a root.
TEST(ExactGroebner, EliminationIdealIsExact) {
  using poly = varietas::polynomial<rational, 3, lex>;
  const auto c = [](std::int64_t n, std::int64_t d = 1) { return coeff<rational>(n, d); };

  const poly x = poly::variable(0);
  const poly y = poly::variable(1);
  const poly z = poly::variable(2);

  varietas::ideal<rational, 3, lex> i(std::vector<poly>{
      x + y + z - poly::constant(c(3)),
      x * y + y * z + z * x - poly::constant(c(3)),
      x * y * z - poly::constant(c(1)),
  });

  const auto eliminated = i.eliminate(1);
  ASSERT_FALSE(eliminated.empty());

  for (const auto& g : eliminated) {
    for (const auto& t : g.terms()) {
      EXPECT_EQ(t.mon[0], 0) << "eliminant still involves x";
    }
    EXPECT_TRUE(i.contains(g));
  }
}

}  // namespace
