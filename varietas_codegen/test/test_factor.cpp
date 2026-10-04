// Factorisation over Q, squarefree, univariate and multivariate.
//
// A factorisation has two properties and the tests check both separately. The
// product has to give back the polynomial, which is easy to test and catches
// any loss in the lifting. The factors have to be irreducible, which is not
// easy to test in general, so the inputs are chosen where it is known: the
// cyclotomic polynomials, which split modulo most primes and must be put back
// together; Eisenstein polynomials, irreducible by a criterion that has
// nothing to do with how this code works; and the polynomials of
// Swinnerton-Dyer, which are irreducible over Q and split into linear and
// quadratic factors modulo every prime, so that every subset of the modular
// factors has to be tried and rejected. In several variables the same roles
// are played by x^n - y^m, irreducible exactly when n and m are coprime, and by
// polynomials of degree one in some variable with coprime coefficients.

#include <algorithm>
#include <array>
#include <cstddef>
#include <random>
#include <vector>

#include <gtest/gtest.h>

#include "varietas/codegen/factor.hpp"
#include "varietas/codegen/modular_gcd.hpp"
#include "varietas/codegen/rational.hpp"
#include "varietas/core/order/grevlex.hpp"
#include "varietas/core/order/lex.hpp"
#include "varietas/core/polynomial.hpp"

namespace {

using varietas::grevlex;
using varietas::make_rational;
using varietas::rational;

using upoly = varietas::polynomial<rational, 1, grevlex>;
using tpoly = varietas::polynomial<rational, 3, grevlex>;

upoly x() { return upoly::variable(0); }
upoly c(std::int64_t n, std::int64_t d = 1) { return upoly::constant(make_rational(n, d)); }

// The polynomial with coefficients a_0, a_1, ... in increasing degree.
upoly from_coefficients(const std::vector<std::int64_t>& a) {
  std::vector<upoly::term> terms;
  for (std::size_t k = 0; k < a.size(); ++k) {
    if (a[k] != 0) {
      terms.push_back({varietas::monomial<1>::variable(0, static_cast<std::uint16_t>(k)),
                       make_rational(a[k])});
    }
  }
  return upoly(std::move(terms));
}

upoly power(const upoly& p, unsigned e) {
  upoly r = c(1);
  for (unsigned i = 0; i < e; ++i) {
    r = r * p;
  }
  return r;
}

template <class Poly>
Poly expand(const varietas::factorisation<Poly>& f) {
  Poly r = Poly::constant(f.unit);
  for (const auto& [base, multiplicity] : f.factors) {
    for (unsigned i = 0; i < multiplicity; ++i) {
      r = r * base;
    }
  }
  return r;
}

// Every factor monic, nonconstant, and distinct from every other.
template <class Poly>
void expect_well_formed(const varietas::factorisation<Poly>& f) {
  for (std::size_t i = 0; i < f.factors.size(); ++i) {
    const Poly& g = f.factors[i].base;
    EXPECT_GT(g.degree(), 0u);
    EXPECT_EQ(g.leading_coefficient(), rational(1));
    EXPECT_GE(f.factors[i].multiplicity, 1u);
    for (std::size_t j = i + 1; j < f.factors.size(); ++j) {
      EXPECT_EQ(varietas::modular_gcd(g, f.factors[j].base).degree(), 0u);
    }
  }
}

// The factorisation as a multiset, compared without regard to order.
bool same_factors(const varietas::factorisation<upoly>& actual,
                  std::vector<std::pair<upoly, unsigned>> expected) {
  if (actual.factors.size() != expected.size()) {
    return false;
  }
  for (const auto& [base, multiplicity] : actual.factors) {
    const auto it = std::find_if(expected.begin(), expected.end(), [&](const auto& e) {
      return e.first.monic() == base && e.second == multiplicity;
    });
    if (it == expected.end()) {
      return false;
    }
    expected.erase(it);
  }
  return true;
}

// The n-th cyclotomic polynomial, as x^n - 1 divided by the cyclotomic
// polynomials of the proper divisors of n.
upoly cyclotomic(int n) {
  upoly p = power(x(), static_cast<unsigned>(n)) - c(1);
  for (int d = 1; d < n; ++d) {
    if (n % d == 0) {
      p = varietas::detail::divide_exact(p, cyclotomic(d));
    }
  }
  return p;
}

TEST(FactorUnivariate, ADifferenceOfSquares) {
  const auto f = varietas::factor_univariate(x() * x() - c(1));
  EXPECT_EQ(f.unit, rational(1));
  EXPECT_TRUE(same_factors(f, {{x() - c(1), 1}, {x() + c(1), 1}}));
}

TEST(FactorUnivariate, ConstantsAndLinearPolynomials) {
  const auto constant = varietas::factor_univariate(c(-7, 3));
  EXPECT_EQ(constant.unit, make_rational(-7, 3));
  EXPECT_TRUE(constant.factors.empty());

  const auto linear = varietas::factor_univariate(c(6) * x() + c(4));
  EXPECT_EQ(linear.unit, rational(6));
  EXPECT_TRUE(same_factors(linear, {{x() + c(2, 3), 1}}));
}

// x^4 + 1 is irreducible over Q and reducible modulo every prime, which makes
// it the smallest polynomial on which the recombination has to reject every
// candidate.
TEST(FactorUnivariate, XToTheFourthPlusOneIsIrreducible) {
  const upoly f = power(x(), 4) + c(1);
  const auto factors = varietas::factor_univariate(f);
  EXPECT_TRUE(same_factors(factors, {{f, 1}}));
}

TEST(FactorUnivariate, XToTheNMinusOneSplitsIntoCyclotomicPolynomials) {
  for (int n = 1; n <= 36; ++n) {
    std::vector<std::pair<upoly, unsigned>> expected;
    for (int d = 1; d <= n; ++d) {
      if (n % d == 0) {
        expected.emplace_back(cyclotomic(d), 1);
      }
    }
    const upoly f = power(x(), static_cast<unsigned>(n)) - c(1);
    const auto factors = varietas::factor_univariate(f);
    EXPECT_TRUE(same_factors(factors, expected)) << "n = " << n;
    EXPECT_EQ(expand(factors), f) << "n = " << n;
  }
}

// The product of x + sqrt(2) + sqrt(3) + sqrt(5) over all choices of signs,
// and the same with sqrt(7) adjoined. Each is irreducible over Q, and modulo
// any prime every square root either exists or generates the quadratic
// extension, so the factors there have degree one or two: at least four and
// at least eight of them, every subset of which is a candidate that fails.
TEST(FactorUnivariate, SwinnertonDyerPolynomialsAreIrreducible) {
  const upoly s3 = from_coefficients({576, 0, -960, 0, 352, 0, -40, 0, 1});
  const upoly s4 = from_coefficients({46225, 0, -5596840, 0, 13950764, 0, -7453176, 0, 1513334,
                                      0, -141912, 0, 6476, 0, -136, 0, 1});
  for (const upoly& s : {s3, s4}) {
    const auto before = varietas::factor_counters();
    const auto factors = varietas::factor_univariate(s);
    const auto after = varietas::factor_counters();
    EXPECT_TRUE(same_factors(factors, {{s, 1}}));
    EXPECT_GE(after.modular_factors - before.modular_factors,
              static_cast<std::size_t>(s.degree() / 2));
  }
}

TEST(FactorUnivariate, AProductOfSwinnertonDyerPolynomialsSeparates) {
  const upoly s3 = from_coefficients({576, 0, -960, 0, 352, 0, -40, 0, 1});
  // The same construction for sqrt(2), sqrt(3), sqrt(7).
  const upoly t3 = from_coefficients({400, 0, -1728, 0, 536, 0, -48, 0, 1});
  ASSERT_EQ(varietas::modular_gcd(s3, t3).degree(), 0u);
  const auto factors = varietas::factor_univariate(s3 * t3);
  EXPECT_TRUE(same_factors(factors, {{s3, 1}, {t3, 1}}));
}

// x^n + q(...) with every lower coefficient divisible by q and the constant
// term not by q^2 is irreducible by Eisenstein's criterion, so a product of
// such polynomials has exactly those factors.
TEST(FactorUnivariate, ProductsOfEisensteinPolynomialsSeparateIntoThem) {
  std::mt19937 rng(23);
  const std::array<int, 4> primes{2, 3, 5, 7};
  std::uniform_int_distribution<int> degree(1, 6);
  std::uniform_int_distribution<int> small(-4, 4);
  std::uniform_int_distribution<int> lead(1, 5);
  for (int trial = 0; trial < 30; ++trial) {
    std::vector<std::pair<upoly, unsigned>> expected;
    upoly f = c(1);
    for (int k = 0; k < 3; ++k) {
      const int q = primes[static_cast<std::size_t>(trial + k) % primes.size()];
      const int n = degree(rng);
      std::vector<std::int64_t> a(static_cast<std::size_t>(n) + 1);
      for (int i = 1; i < n; ++i) {
        a[static_cast<std::size_t>(i)] = q * small(rng);
      }
      int constant = 0;
      while (constant % q == 0) {
        constant = small(rng);
      }
      a[0] = q * constant;
      int l = 0;
      while (l % q == 0) {
        l = lead(rng);
      }
      a[static_cast<std::size_t>(n)] = l;
      const upoly e = from_coefficients(a);
      f = f * e;
      auto it = std::find_if(expected.begin(), expected.end(),
                             [&](const auto& x) { return x.first.monic() == e.monic(); });
      if (it == expected.end()) {
        expected.emplace_back(e, 1);
      } else {
        ++it->second;
      }
    }
    const auto factors = varietas::factor_univariate(f);
    EXPECT_TRUE(same_factors(factors, expected)) << "trial " << trial;
    EXPECT_EQ(expand(factors), f) << "trial " << trial;
    expect_well_formed(factors);
  }
}

TEST(FactorUnivariate, MultiplicitiesAndRationalCoefficients) {
  const upoly f = c(3, 2) * power(x() - c(1, 3), 2) * power(x() * x() + c(2), 5) *
                  power(x(), 3) * (c(7) * x() * x() * x() - c(2));
  const auto factors = varietas::factor_univariate(f);
  EXPECT_EQ(factors.unit, make_rational(21, 2));
  EXPECT_TRUE(same_factors(factors, {{x() - c(1, 3), 2},
                                     {x() * x() + c(2), 5},
                                     {x(), 3},
                                     {c(7) * x() * x() * x() - c(2), 1}}));
  EXPECT_EQ(expand(factors), f);
  expect_well_formed(factors);
  for (std::size_t i = 1; i < factors.factors.size(); ++i) {
    EXPECT_LE(factors.factors[i - 1].multiplicity, factors.factors[i].multiplicity);
  }
}

// Roots far apart in size, so that the bound and the precision lifted to are
// large and the symmetric reading of the lifted factors has to be exact.
TEST(FactorUnivariate, LargeCoefficients) {
  const rational big = varietas::rational_from_string("1000000000000000000000000000057");
  const upoly f = (x() - upoly::constant(big)) * (c(7) * x() + c(1)) *
                  (x() * x() - upoly::constant(big * big + 1)) * (x() * x() + c(1, 1000));
  const auto factors = varietas::factor_univariate(f);
  EXPECT_TRUE(same_factors(factors, {{x() - upoly::constant(big), 1},
                                     {x() + c(1, 7), 1},
                                     {x() * x() - upoly::constant(big * big + 1), 1},
                                     {x() * x() + c(1, 1000), 1}}));
  EXPECT_EQ(expand(factors), f);
}

// A polynomial in one of several variables, under an order that is not the
// one the factors are naturally written in.
TEST(FactorUnivariate, OneVariableOfSeveral) {
  using poly = varietas::polynomial<rational, 3, varietas::lex>;
  const poly y = poly::variable(1);
  const poly four = poly::constant(rational(4));
  const auto factors = varietas::factor_univariate(y * y * y * y - four * four);
  ASSERT_EQ(factors.factors.size(), 3u);
  EXPECT_EQ(expand(factors), y * y * y * y - four * four);
  expect_well_formed(factors);
  EXPECT_EQ(factors.factors[0].base.degree(), 1u);
  EXPECT_EQ(factors.factors[1].base.degree(), 1u);
  EXPECT_EQ(factors.factors[2].base, y * y + four);
}

tpoly random_trivariate(std::mt19937& rng, int degree, int terms) {
  using mon = varietas::monomial<3>;
  std::uniform_int_distribution<int> exponent(0, degree);
  std::uniform_int_distribution<int> coefficient(-9, 9);
  std::vector<tpoly::term> out;
  for (int i = 0; i < terms; ++i) {
    std::array<mon::exponent_type, 3> e{};
    for (auto& k : e) {
      k = static_cast<mon::exponent_type>(exponent(rng));
    }
    const int k = coefficient(rng);
    if (k != 0) {
      out.push_back({mon(e), make_rational(k)});
    }
  }
  return tpoly(std::move(out));
}

TEST(SquarefreeDecomposition, SeparatesMultiplicitiesInSeveralVariables) {
  std::mt19937 rng(5);
  for (int trial = 0; trial < 15; ++trial) {
    tpoly a;
    tpoly b;
    tpoly d;
    do {
      a = random_trivariate(rng, 2, 4);
      b = random_trivariate(rng, 2, 4);
      d = random_trivariate(rng, 1, 3);
    } while (a.degree() == 0 || b.degree() == 0 || d.degree() == 0);
    const tpoly f = tpoly::constant(make_rational(-5, 3)) * a * b * b * d * d * d;
    const auto decomposition = varietas::squarefree_decomposition(f);
    EXPECT_EQ(expand(decomposition), f) << "trial " << trial;
    expect_well_formed(decomposition);

    // Each part is squarefree: it shares no factor with any of its partial
    // derivatives that involves the variable differentiated.
    for (const auto& [g, multiplicity] : decomposition.factors) {
      for (std::size_t v = 0; v < 3; ++v) {
        const tpoly dg = varietas::factor_detail::derivative(g, v);
        if (!dg.is_zero()) {
          EXPECT_LE(varietas::detail::degree_in(varietas::modular_gcd(g, dg), v), 0)
              << "trial " << trial;
        }
      }
    }
  }
}

TEST(SquarefreeDecomposition, AFactorFreeOfTheFirstVariableIsFound) {
  const tpoly x0 = tpoly::variable(0);
  const tpoly x1 = tpoly::variable(1);
  const tpoly x2 = tpoly::variable(2);
  const tpoly one = tpoly::constant(rational(1));
  // (x1^2 + 1)^2 is content with respect to x0, x0^3 and (x0 x2 + 1) are not.
  const tpoly f = (x1 * x1 + one) * (x1 * x1 + one) * x0 * x0 * x0 * (x0 * x2 + one) *
                  (x2 - x1) * (x2 - x1);
  const auto decomposition = varietas::squarefree_decomposition(f);
  EXPECT_EQ(expand(decomposition), f);
  ASSERT_EQ(decomposition.factors.size(), 3u);
  EXPECT_EQ(decomposition.factors[0].multiplicity, 1u);
  EXPECT_EQ(decomposition.factors[0].base, (x0 * x2 + one).monic());
  EXPECT_EQ(decomposition.factors[1].multiplicity, 2u);
  EXPECT_EQ(decomposition.factors[1].base, ((x1 * x1 + one) * (x2 - x1)).monic());
  EXPECT_EQ(decomposition.factors[2].multiplicity, 3u);
  EXPECT_EQ(decomposition.factors[2].base, x0);
}

// Multivariate.

template <class Poly>
bool same_factors_multivariate(const varietas::factorisation<Poly>& actual,
                               std::vector<std::pair<Poly, unsigned>> expected) {
  if (actual.factors.size() != expected.size()) {
    return false;
  }
  for (const auto& [base, multiplicity] : actual.factors) {
    const auto it = std::find_if(expected.begin(), expected.end(), [&](const auto& e) {
      return e.first.monic() == base && e.second == multiplicity;
    });
    if (it == expected.end()) {
      return false;
    }
    expected.erase(it);
  }
  return true;
}

using bpoly = varietas::polynomial<rational, 2, grevlex>;

TEST(Factor, XToTheNMinusYToTheM) {
  const bpoly x = bpoly::variable(0);
  const bpoly y = bpoly::variable(1);
  const auto pow = [](const bpoly& p, unsigned e) {
    bpoly r = bpoly::constant(rational(1));
    for (unsigned i = 0; i < e; ++i) {
      r = r * p;
    }
    return r;
  };
  // Coprime exponents: irreducible.
  for (const auto& [n, m] : std::vector<std::pair<unsigned, unsigned>>{{3, 2}, {5, 3}, {7, 4}}) {
    const bpoly f = pow(x, n) - pow(y, m);
    EXPECT_TRUE(same_factors_multivariate(varietas::factor(f), {{f, 1}})) << n << ", " << m;
  }
  const auto fourth = varietas::factor(pow(x, 4) - pow(y, 4));
  EXPECT_TRUE(same_factors_multivariate(fourth, {{x - y, 1}, {x + y, 1}, {x * x + y * y, 1}}));
  const auto sixth = varietas::factor(pow(x, 6) - pow(y, 4));
  EXPECT_TRUE(
      same_factors_multivariate(sixth, {{pow(x, 3) - pow(y, 2), 1}, {pow(x, 3) + pow(y, 2), 1}}));
}

// The boundary of the workspace of a planar arm with links 2 and 1 is the
// pair of circles of radii 3 and 1, and its implicit equation is their
// product.
TEST(Factor, TheWorkspaceBoundaryOfAPlanarArmIsTwoCircles) {
  const bpoly x = bpoly::variable(0);
  const bpoly y = bpoly::variable(1);
  const bpoly r2 = x * x + y * y;
  const bpoly outer = r2 - bpoly::constant(rational(9));
  const bpoly inner = r2 - bpoly::constant(rational(1));
  const bpoly f = bpoly::constant(make_rational(1, 4)) * outer * inner;
  const auto factors = varietas::factor(f);
  EXPECT_EQ(factors.unit, make_rational(1, 4));
  EXPECT_TRUE(same_factors_multivariate(factors, {{outer, 1}, {inner, 1}}));
}

// a x_k + b with a and b coprime and free of x_k is irreducible, whatever a
// and b are. Products of three, each of degree one in a different variable,
// have to come apart into exactly those, with multiplicities and a content.
TEST(Factor, ProductsOfPolynomialsLinearInOneVariable) {
  std::mt19937 rng(41);
  const tpoly one = tpoly::constant(rational(1));
  for (int trial = 0; trial < 20; ++trial) {
    std::vector<std::pair<tpoly, unsigned>> expected;
    tpoly f = tpoly::constant(make_rational(-3, 7));
    for (std::size_t k = 0; k < 3; ++k) {
      tpoly a;
      tpoly b;
      // Coefficients free of x_k, coprime, and b nonzero, or the factor would
      // be x_k times a.
      do {
        a = random_trivariate(rng, 2, 3);
        b = random_trivariate(rng, 2, 4);
        std::vector<tpoly::term> keep_a;
        std::vector<tpoly::term> keep_b;
        for (const auto& t : a.terms()) {
          if (t.mon[k] == 0) {
            keep_a.push_back(t);
          }
        }
        for (const auto& t : b.terms()) {
          if (t.mon[k] == 0) {
            keep_b.push_back(t);
          }
        }
        a = tpoly(std::move(keep_a));
        b = tpoly(std::move(keep_b));
      } while (a.is_zero() || b.is_zero() || varietas::modular_gcd(a, b).degree() != 0);
      const tpoly g = a * tpoly::variable(k) + b;
      const unsigned multiplicity = (trial + static_cast<int>(k)) % 3 == 0 ? 2u : 1u;
      for (unsigned i = 0; i < multiplicity; ++i) {
        f = f * g;
      }
      expected.emplace_back(g, multiplicity);
    }
    f = f * (tpoly::variable(1) * tpoly::variable(1) + one);
    expected.emplace_back(tpoly::variable(1) * tpoly::variable(1) + one, 1u);
    const auto factors = varietas::factor(f);
    EXPECT_EQ(expand(factors), f) << "trial " << trial;
    expect_well_formed(factors);
    // Two of the factors may coincide up to a constant, which merges them.
    std::vector<std::pair<tpoly, unsigned>> merged;
    for (const auto& [g, m] : expected) {
      auto it = std::find_if(merged.begin(), merged.end(),
                             [&](const auto& e) { return e.first.monic() == g.monic(); });
      if (it == merged.end()) {
        merged.emplace_back(g, m);
      } else {
        it->second += m;
      }
    }
    EXPECT_TRUE(same_factors_multivariate(factors, merged)) << "trial " << trial;
  }
}

// Every leading coefficient vanishes at the origin, so the point the lifting
// runs at has to be somewhere else and the factors moved back from it.
TEST(Factor, TheOriginDoesNotQualify) {
  const tpoly x = tpoly::variable(0);
  const tpoly y = tpoly::variable(1);
  const tpoly z = tpoly::variable(2);
  const tpoly one = tpoly::constant(rational(1));
  const tpoly g = x * y + z;
  const tpoly h = x * z * z + y * y + one;
  const tpoly k = y * z * x * x + x + z;
  const auto factors = varietas::factor(g * h * k);
  EXPECT_TRUE(same_factors_multivariate(factors, {{g, 1}, {h, 1}, {k, 1}}));
}

// At y = 0, x^2 - y^3 - 1 is x^2 - 1 and splits, though the polynomial does
// not, so with the point forced to the origin the leading coefficient cannot
// be imposed on a division into halves and the factors have to be lifted as
// power series and recombined. The search for a point would have avoided
// this; the test takes the decision away from it.
TEST(Factor, RecombinationWhenThePointSplitsAnIrreducibleFactor) {
  using varietas::factor_detail::lift_and_recombine;
  const bpoly x = bpoly::variable(0);
  const bpoly y = bpoly::variable(1);
  const bpoly one = bpoly::constant(rational(1));
  const bpoly two = bpoly::constant(rational(2));
  const bpoly g = x * x - y * y * y - one;
  const bpoly h = (y + one) * x - y - two;

  const auto before = varietas::factor_counters();
  const auto alone = lift_and_recombine(g, 0, {x - one, x + one});
  ASSERT_EQ(alone.size(), 1u);
  EXPECT_EQ(alone.front().monic(), g.monic());

  // f(x, 0) = (x - 1)(x + 1)(x - 2): three factors, two true ones.
  const auto both = lift_and_recombine(g * h, 0, {x - one, x + one, x - two});
  ASSERT_EQ(both.size(), 2u);
  const bool order = both[0].monic() == g.monic();
  EXPECT_EQ(both[order ? 0 : 1].monic(), g.monic());
  EXPECT_EQ(both[order ? 1 : 0].monic(), h.monic());
  const auto after = varietas::factor_counters();
  EXPECT_GE(after.series_lifts - before.series_lifts, 2u);
}

// Nothing to lift: a content with respect to the main variable, a power of a
// variable, and a unit, around one irreducible factor.
TEST(Factor, ContentsPowersAndUnits) {
  const tpoly x = tpoly::variable(0);
  const tpoly y = tpoly::variable(1);
  const tpoly z = tpoly::variable(2);
  const tpoly one = tpoly::constant(rational(1));
  const tpoly f = tpoly::constant(make_rational(5, 2)) * z * z * z * (y * y - one) *
                  (y * y - one) * (x * x + y * z + one);
  const auto factors = varietas::factor(f);
  EXPECT_EQ(factors.unit, make_rational(5, 2));
  EXPECT_TRUE(same_factors_multivariate(
      factors, {{z, 3}, {y - one, 2}, {y + one, 2}, {x * x + y * z + one, 1}}));
  EXPECT_EQ(expand(factors), f);
}

}  // namespace
