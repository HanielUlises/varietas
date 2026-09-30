// The modular gcd, against the subresultant one and against itself.
//
// modular_gcd replaces polynomial_gcd on the path where the parametric solve
// spends its time, and it gets there by a chain of steps none of which is
// correct on its own: an image modulo a prime, images at evaluation points,
// interpolation, rational reconstruction. Only the trial division at the end is
// exact. So the tests do not look inside the chain. They check the answer, on
// the cases that are hard for each step: a common factor that exists, one that
// does not, a factor in one variable only, rational coefficients with large
// denominators, and inputs built so that an early prime or an early point is
// unlucky.

#include <array>
#include <cstdint>
#include <random>
#include <vector>

#include <gtest/gtest.h>

#include "varietas/codegen/modular_gcd.hpp"
#include "varietas/codegen/rational.hpp"
#include "varietas/core/gcd.hpp"
#include "varietas/core/ideal/division.hpp"
#include "varietas/core/order/grevlex.hpp"
#include "varietas/core/order/lex.hpp"
#include "varietas/core/polynomial.hpp"

namespace {

using varietas::grevlex;
using varietas::make_rational;
using varietas::rational;

template <std::size_t N>
using poly = varietas::polynomial<rational, N, grevlex>;

template <std::size_t N, class Order = grevlex>
varietas::polynomial<rational, N, Order> random_polynomial(std::mt19937& rng, int degree,
                                                           int terms, int magnitude) {
  using result = varietas::polynomial<rational, N, Order>;
  using mon = varietas::monomial<N>;
  std::uniform_int_distribution<int> exponent(0, degree);
  std::uniform_int_distribution<int> coefficient(-magnitude, magnitude);
  std::vector<typename result::term> out;
  for (int i = 0; i < terms; ++i) {
    std::array<typename mon::exponent_type, N> e{};
    int total = 0;
    for (std::size_t v = 0; v < N; ++v) {
      e[v] = static_cast<typename mon::exponent_type>(exponent(rng));
      total += e[v];
    }
    if (total > degree) {
      continue;
    }
    const int c = coefficient(rng);
    if (c != 0) {
      out.push_back({mon(e), make_rational(c)});
    }
  }
  if (out.empty()) {
    out.push_back({mon::one(), make_rational(1)});
  }
  return result(std::move(out));
}

template <std::size_t N>
bool divides(const poly<N>& d, const poly<N>& p) {
  return varietas::divide(p, std::vector<poly<N>>{d}).remainder.is_zero();
}

TEST(ModularGcd, AgreesWithTheSubresultantGcdOnRandomTrivariateInputs) {
  std::mt19937 rng(17);
  for (int trial = 0; trial < 40; ++trial) {
    const auto g = random_polynomial<3>(rng, 3, 6, 9);
    const auto a = random_polynomial<3>(rng, 3, 6, 9) * g;
    const auto b = random_polynomial<3>(rng, 3, 6, 9) * g;

    const auto expected = varietas::polynomial_gcd(a, b);
    const auto actual = varietas::modular_gcd(a, b);
    EXPECT_EQ(actual, expected) << "trial " << trial;
    EXPECT_TRUE(divides(g.monic(), actual)) << "trial " << trial;
  }
}

TEST(ModularGcd, CoprimeInputsGiveOne) {
  std::mt19937 rng(3);
  for (int trial = 0; trial < 20; ++trial) {
    const auto a = random_polynomial<3>(rng, 4, 8, 20);
    const auto b = random_polynomial<3>(rng, 4, 8, 20);
    const auto expected = varietas::polynomial_gcd(a, b);
    EXPECT_EQ(varietas::modular_gcd(a, b), expected) << "trial " << trial;
  }
}

TEST(ModularGcd, TheResultIsMonicUnderTheOrderOfThePolynomial) {
  const poly<2> x = poly<2>::variable(0);
  const poly<2> y = poly<2>::variable(1);
  const poly<2> three = poly<2>::constant(make_rational(3));
  const auto g = three * x * y + y * y + poly<2>::constant(make_rational(7, 2));
  const auto a = g * (x + y);
  const auto b = g * (x - three);
  const auto result = varietas::modular_gcd(a, b);
  EXPECT_EQ(result, g.monic());
  EXPECT_EQ(result.leading_coefficient(), make_rational(1));
}

// A factor in one variable only lives entirely in the content at one level of
// the recursion, and is invisible at every evaluation of that variable.
TEST(ModularGcd, AFactorInOneVariableIsRecoveredFromTheContent) {
  const poly<3> x = poly<3>::variable(0);
  const poly<3> y = poly<3>::variable(1);
  const poly<3> z = poly<3>::variable(2);
  const poly<3> one = poly<3>::constant(make_rational(1));

  const auto in_z = z * z + z + poly<3>::constant(make_rational(5));
  const auto a = in_z * (x * y + one) * (x + z);
  const auto b = in_z * (x * y - one) * (y + z);
  EXPECT_EQ(varietas::modular_gcd(a, b), in_z.monic());

  const auto in_x = x * x * x - poly<3>::constant(make_rational(2));
  EXPECT_EQ(varietas::modular_gcd(in_x * (y + z), in_x * (y - z)), in_x.monic());
}

// Coefficients whose numerators and denominators are far larger than one prime,
// so that rational reconstruction needs several images before it can succeed.
TEST(ModularGcd, LargeRationalCoefficientsNeedSeveralPrimesAndAreRecoveredExactly) {
  const poly<2> x = poly<2>::variable(0);
  const poly<2> y = poly<2>::variable(1);
  const rational big = varietas::rational_from_string("123456789012345678901/98765432109876543");
  const auto g = x * y + poly<2>::constant(big) * x + poly<2>::constant(make_rational(1, 3));
  const auto a = g * (x + y + poly<2>::constant(make_rational(2)));
  const auto b = g * (x * x - y);
  EXPECT_EQ(varietas::modular_gcd(a, b), g.monic());
}

// An evaluation point where the two images share a spurious factor: at z = 1
// the cofactors x + z and x + 1 coincide, and the image gcd there is too large.
// The largest-leading-monomial rule has to throw that image away.
TEST(ModularGcd, AnUnluckyEvaluationPointIsDiscarded) {
  const poly<3> x = poly<3>::variable(0);
  const poly<3> y = poly<3>::variable(1);
  const poly<3> z = poly<3>::variable(2);
  const poly<3> one = poly<3>::constant(make_rational(1));
  const auto g = x * y + z;
  const auto a = g * (x + z);
  const auto b = g * (x + one);
  EXPECT_EQ(varietas::modular_gcd(a, b), varietas::polynomial_gcd(a, b));
  EXPECT_EQ(varietas::modular_gcd(a, b), g.monic());
}

// A prime that divides the leading coefficient of an operand could hide the
// leading term of the gcd, and has to be skipped. 2^31 - 1 is the first prime
// the algorithm tries.
TEST(ModularGcd, APrimeDividingALeadingCoefficientIsSkipped) {
  const poly<2> x = poly<2>::variable(0);
  const poly<2> y = poly<2>::variable(1);
  const rational p = make_rational(2147483647);
  const auto g = x + y;
  const auto a = g * (poly<2>::constant(p) * x * x + y);
  const auto b = g * (x - y + poly<2>::constant(p));
  EXPECT_EQ(varietas::modular_gcd(a, b), g.monic());
}

// The cofactors are what normalisation keeps, so they are checked as closely as
// the gcd: a = g * (a / g) exactly, in every case including the degenerate ones.
TEST(ModularGcd, CofactorsMultiplyBackToTheOperands) {
  std::mt19937 rng(23);
  for (int trial = 0; trial < 20; ++trial) {
    const auto f = random_polynomial<3>(rng, 3, 5, 9);
    const auto a = random_polynomial<3>(rng, 3, 5, 9) * f;
    const auto b = random_polynomial<3>(rng, 3, 5, 9) *
                   (trial % 2 == 0 ? f : poly<3>::constant(make_rational(1)));
    const auto r = varietas::modular_gcd_with_cofactors(a, b);
    EXPECT_EQ(r.gcd * r.a_over_gcd, a) << "trial " << trial;
    EXPECT_EQ(r.gcd * r.b_over_gcd, b) << "trial " << trial;
    EXPECT_EQ(r.gcd, varietas::polynomial_gcd(a, b)) << "trial " << trial;
  }

  const poly<2> x = poly<2>::variable(0);
  const auto f = poly<2>::constant(make_rational(4)) * x - poly<2>::constant(make_rational(1));
  const auto with_zero = varietas::modular_gcd_with_cofactors(f, poly<2>());
  EXPECT_EQ(with_zero.gcd * with_zero.a_over_gcd, f);
  EXPECT_TRUE(with_zero.b_over_gcd.is_zero());
  const auto zero_first = varietas::modular_gcd_with_cofactors(poly<2>(), f);
  EXPECT_EQ(zero_first.gcd * zero_first.b_over_gcd, f);
  EXPECT_TRUE(zero_first.a_over_gcd.is_zero());
}

TEST(ModularGcd, ZeroAndConstantOperands) {
  const poly<2> x = poly<2>::variable(0);
  const poly<2> y = poly<2>::variable(1);
  const auto f = poly<2>::constant(make_rational(4)) * x * y - y;
  EXPECT_EQ(varietas::modular_gcd(f, poly<2>()), f.monic());
  EXPECT_EQ(varietas::modular_gcd(poly<2>(), f), f.monic());
  EXPECT_TRUE(varietas::modular_gcd(poly<2>(), poly<2>()).is_zero());
  EXPECT_EQ(varietas::modular_gcd(f, poly<2>::constant(make_rational(5))),
            poly<2>::constant(make_rational(1)));
}

TEST(ModularGcd, WorksInOneVariableAndUnderOtherOrders) {
  using lex1 = varietas::polynomial<rational, 1, varietas::lex>;
  const lex1 t = lex1::variable(0);
  const lex1 one = lex1::constant(make_rational(1));
  const auto g = t * t + one;
  EXPECT_EQ(varietas::modular_gcd(g * (t - one), g * (t + one)), g);

  using lex3 = varietas::polynomial<rational, 3, varietas::lex>;
  std::mt19937 rng(11);
  for (int trial = 0; trial < 10; ++trial) {
    const auto gg = random_polynomial<3, varietas::lex>(rng, 2, 5, 7);
    const auto a = random_polynomial<3, varietas::lex>(rng, 2, 5, 7) * gg;
    const auto b = random_polynomial<3, varietas::lex>(rng, 2, 5, 7) * gg;
    EXPECT_EQ(varietas::modular_gcd(a, b), varietas::polynomial_gcd(a, b)) << "trial " << trial;
  }
}

// The size the parametric solve reaches, where the subresultant sequence takes
// seconds. This is here to keep the fast path fast, not to time it: the bound
// is generous and a regression to remainder sequences would blow through it.
TEST(ModularGcd, DenseTrivariateOperandsOfSeveralHundredTermsAreCheap) {
  std::mt19937 rng(7);
  const auto factor = random_polynomial<3>(rng, 6, 200, 50);
  const auto a = random_polynomial<3>(rng, 6, 200, 50) * factor;
  const auto b = random_polynomial<3>(rng, 6, 200, 50) * factor;
  ASSERT_GT(a.size(), 300u);
  const auto result = varietas::modular_gcd(a, b);
  EXPECT_TRUE(divides(result, a));
  EXPECT_TRUE(divides(result, b));
  EXPECT_TRUE(divides(factor.monic(), result));
}

}  // namespace
