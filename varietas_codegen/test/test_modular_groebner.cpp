// Bases over Q lifted from Z/p, against Buchberger over Q.
//
// The reduced basis is unique, so the lifted one has to be Buchberger's,
// coefficient for coefficient, wherever both are computed. Beyond that: a
// saturation is lifted and checked when it is zero-dimensional and computed
// exactly when it is not, the unit ideal is found, and the answer does not
// depend on how many images are computed at once.

#include <array>
#include <cstdint>
#include <random>
#include <vector>

#include <gtest/gtest.h>

#include "varietas/codegen/modular_groebner.hpp"
#include "varietas/codegen/rational.hpp"
#include "varietas/core/ideal/buchberger.hpp"
#include "varietas/core/ideal/saturation.hpp"
#include "varietas/core/order/grevlex.hpp"
#include "varietas/core/order/lex.hpp"
#include "varietas/core/quotient/quotient_basis.hpp"

namespace {

using varietas::grevlex;
using varietas::lex;
using varietas::make_rational;
using varietas::monomial;
using varietas::polynomial;
using varietas::rational;

template <std::size_t N, class Order>
std::vector<polynomial<rational, N, Order>> random_system(std::mt19937_64& rng, int count,
                                                          int degree, int terms) {
  using poly = polynomial<rational, N, Order>;
  std::uniform_int_distribution<long> numerator(-40, 40);
  std::uniform_int_distribution<long> denominator(1, 12);
  std::vector<poly> out;
  for (int k = 0; k < count; ++k) {
    std::vector<typename poly::term> ts;
    for (int t = 0; t < terms; ++t) {
      std::array<std::uint16_t, N> e{};
      int left = std::uniform_int_distribution<int>(0, degree)(rng);
      for (std::size_t v = 0; v < N && left > 0; ++v) {
        const int x = std::uniform_int_distribution<int>(0, left)(rng);
        e[v] = static_cast<std::uint16_t>(x);
        left -= x;
      }
      long n = numerator(rng);
      if (n == 0) {
        n = 1;
      }
      ts.push_back({monomial<N>(e), make_rational(n, denominator(rng))});
    }
    out.emplace_back(std::move(ts));
  }
  return out;
}

TEST(ModularGroebner, LiftsBuchbergersBasisOfZeroDimensionalSystems) {
  std::mt19937_64 rng(11);
  std::size_t lifted = 0;
  for (int k = 0; k < 30; ++k) {
    const auto generators = random_system<3, grevlex>(rng, 3, 3, 4);
    varietas::modular::statistics s;
    const auto basis = varietas::modular::groebner_basis(generators, &s);
    ASSERT_EQ(basis, varietas::groebner_basis(generators)) << "system " << k;
    lifted += s.lifted ? 1 : 0;
  }
  EXPECT_EQ(lifted, 30u);
}

TEST(ModularGroebner, LiftsUnderLex) {
  std::mt19937_64 rng(12);
  for (int k = 0; k < 20; ++k) {
    const auto generators = random_system<3, lex>(rng, 3, 2, 3);
    ASSERT_EQ(varietas::modular::groebner_basis(generators), varietas::groebner_basis(generators))
        << "system " << k;
  }
}

// A plain basis needs nothing zero-dimensional to be checked, so a curve is
// lifted as well.
TEST(ModularGroebner, LiftsPositiveDimensionalSystems) {
  std::mt19937_64 rng(13);
  for (int k = 0; k < 20; ++k) {
    const auto generators = random_system<3, grevlex>(rng, 2, 3, 3);
    varietas::modular::statistics s;
    ASSERT_EQ(varietas::modular::groebner_basis(generators, &s),
              varietas::groebner_basis(generators))
        << "system " << k;
    EXPECT_TRUE(s.lifted);
  }
}

TEST(ModularGroebner, FindsTheUnitIdeal) {
  using poly = polynomial<rational, 2, grevlex>;
  const poly x = poly::variable(0);
  const poly y = poly::variable(1);
  const auto basis = varietas::modular::groebner_basis(
      std::vector<poly>{x * y - poly::constant(make_rational(1, 3)), x});
  EXPECT_TRUE(varietas::is_unit_ideal(basis));
}

// (x^2 + 1)(y - 2) and (x - 3)(y^2 + 1) vanish at (3, 2) and at the four points
// with x = ±i and y = ±i. Saturating by (x^2 + 1)(y^2 + 1) removes the four,
// which is what the half-angle substitution needs done to its own spurious
// points, and leaves the one.
TEST(ModularGroebner, LiftsAZeroDimensionalSaturationAndChecksIt) {
  using poly = polynomial<rational, 2, grevlex>;
  const poly x = poly::variable(0);
  const poly y = poly::variable(1);
  const poly one = poly::constant(rational(1));
  const std::vector<poly> generators{(x * x + one) * (y - poly::constant(rational(2))),
                                     (x - poly::constant(rational(3))) * (y * y + one)};
  const std::vector<poly> multipliers{x * x + one, y * y + one};

  varietas::modular::statistics s;
  const auto lifted = varietas::modular::saturate_by_product(generators, multipliers, &s);
  EXPECT_TRUE(s.lifted);
  EXPECT_EQ(lifted, varietas::saturate_by_product(generators, multipliers));
  EXPECT_EQ(lifted, (std::vector<poly>{x - poly::constant(rational(3)),
                                       y - poly::constant(rational(2))}));
}

// A saturation whose result is a curve cannot be checked by the unit test on a
// finite quotient, so it is computed exactly, and says so.
TEST(ModularGroebner, ComputesAPositiveDimensionalSaturationExactly) {
  using poly = polynomial<rational, 3, grevlex>;
  const poly x = poly::variable(0);
  const poly y = poly::variable(1);
  const poly z = poly::variable(2);
  const poly one = poly::constant(rational(1));
  const std::vector<poly> generators{(x * x + one) * (y - z), (x - poly::constant(rational(3))) * z};
  const std::vector<poly> multipliers{x * x + one};

  varietas::modular::statistics s;
  const auto basis = varietas::modular::saturate_by_product(generators, multipliers, &s);
  EXPECT_FALSE(s.lifted);
  EXPECT_EQ(basis, varietas::saturate_by_product(generators, multipliers));
}

TEST(ModularGroebner, DoesNotDependOnHowManyImagesAreComputedAtOnce) {
  std::mt19937_64 rng(14);
  const auto generators = random_system<3, grevlex>(rng, 3, 3, 4);
  varietas::modular::options one;
  one.threads = 1;
  varietas::modular::options three;
  three.threads = 3;
  varietas::modular::statistics a;
  varietas::modular::statistics b;
  EXPECT_EQ(varietas::modular::groebner_basis(generators, &a, one),
            varietas::modular::groebner_basis(generators, &b, three));
  EXPECT_EQ(a.primes, b.primes);
  EXPECT_EQ(a.reconstructions, b.reconstructions);
}

// Coefficients far wider than one prime: the lifting has to join many of them.
TEST(ModularGroebner, JoinsAsManyPrimesAsTheCoefficientsNeed) {
  using poly = polynomial<rational, 2, grevlex>;
  const poly x = poly::variable(0);
  const poly y = poly::variable(1);
  const rational big = varietas::rational_from_string("1234567890123456789012345678901234567/"
                                            "98765432109876543210987654321");
  const std::vector<poly> generators{x * x - poly::constant(big) * y,
                                     y * y - poly::constant(big) * x - poly::constant(rational(1))};
  varietas::modular::statistics s;
  EXPECT_EQ(varietas::modular::groebner_basis(generators, &s),
            varietas::groebner_basis(generators));
  EXPECT_TRUE(s.lifted);
  EXPECT_GT(s.primes, 8u);
}

}  // namespace
