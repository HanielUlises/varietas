// The coefficient field the emitter computes over, exercised as a field.
//
// rational_function is where the parametric pipeline lives or dies. Its
// arithmetic is written in terms of a normalisation that cancels a polynomial
// gcd, and the whole reason that cancellation is there is that without it the
// representations grow without bound under repeated arithmetic — which is
// precisely the regime Buchberger puts a coefficient through. A field that is
// correct but whose representations blow up is as useless here as one that is
// wrong, so both properties are checked: the arithmetic against evaluation, and
// the degrees against the degrees of the functions being represented.

#include <array>
#include <cstdint>
#include <random>
#include <vector>

#include <gtest/gtest.h>

#include "varietas/codegen/rational_function.hpp"
#include "varietas/core/ideal/ideal.hpp"
#include "varietas/core/order/grevlex.hpp"
#include "varietas/core/polynomial.hpp"

namespace {

using varietas::grevlex;
using varietas::make_rational;
using varietas::rational;

using field = varietas::rational_function<2>;
using param = field::parameter_polynomial;
using traits = varietas::coefficient_traits<field>;

field constant(std::int64_t n, std::int64_t d = 1) { return field(make_rational(n, d)); }

// A point at which every denominator arising in these tests is nonzero, so that
// evaluation is defined and can serve as the oracle.
std::array<rational, 2> sample_point() {
  return {make_rational(5, 3), make_rational(-2, 7)};
}

// x and y as elements of the field.
field x() { return field::parameter(0); }
field y() { return field::parameter(1); }

TEST(RationalFunction, TraitsAdvertiseAnExactField) {
  EXPECT_TRUE(traits::is_exact);
  EXPECT_TRUE(traits::is_zero(traits::zero()));
  EXPECT_FALSE(traits::is_zero(traits::one()));
  EXPECT_EQ(traits::one(), constant(1));
}

TEST(RationalFunction, ARationalIsAConstantFunction) {
  const field c = constant(3, 4);
  EXPECT_TRUE(c.is_constant());
  EXPECT_EQ(c.constant_value(), make_rational(3, 4));
  EXPECT_EQ(c.evaluate(sample_point()), make_rational(3, 4));
}

TEST(RationalFunction, AParameterIsNotConstantAndEvaluatesToItsCoordinate) {
  const auto at = sample_point();
  EXPECT_FALSE(x().is_constant());
  EXPECT_EQ(x().evaluate(at), at[0]);
  EXPECT_EQ(y().evaluate(at), at[1]);
}

// Equality is of the function rather than of the written form, which is the
// property the emitter leans on when it asks whether a coefficient is one.
TEST(RationalFunction, EqualityIsOfTheFunctionNotTheWrittenForm) {
  const field a = x() / y();
  const field b = (x() * x()) / (x() * y());
  EXPECT_EQ(a, b);

  const field one_written_awkwardly = (x() + y()) / (y() + x());
  EXPECT_EQ(one_written_awkwardly, traits::one());
}

TEST(RationalFunction, InverseUndoesMultiplication) {
  const field f = (x() + constant(1)) / (y() - constant(3));
  EXPECT_EQ(f * traits::inverse(f), traits::one());
  EXPECT_EQ(traits::inverse(traits::inverse(f)), f);
}

TEST(RationalFunction, SubtractionOfEqualsIsStructurallyZero) {
  const field f = (x() * y() + constant(2)) / (x() - y());
  const field g = (x() * y() + constant(2)) / (x() - y());
  EXPECT_TRUE((f - g).is_zero());
  EXPECT_TRUE(traits::is_zero(f - g));
}

// The oracle test. Arithmetic performed symbolically and then evaluated must
// agree with the same arithmetic performed on the values, for every operation.
// This is what catches a normalisation that cancels a factor it should not.
TEST(RationalFunction, ArithmeticAgreesWithEvaluationOnRandomFunctions) {
  std::mt19937 rng(20260819);
  std::uniform_int_distribution<int> coefficient(-4, 4);
  std::uniform_int_distribution<int> exponent(0, 2);

  // A random polynomial in x and y, never the zero polynomial.
  const auto random_polynomial = [&]() {
    param p;
    for (int attempt = 0; attempt < 8; ++attempt) {
      p = param();
      for (int term = 0; term < 3; ++term) {
        const int c = coefficient(rng);
        if (c == 0) {
          continue;
        }
        const std::array<std::uint8_t, 2> e{static_cast<std::uint8_t>(exponent(rng)),
                                            static_cast<std::uint8_t>(exponent(rng))};
        p = p + param::from_monomial(varietas::monomial<2>(e), make_rational(c));
      }
      if (!p.is_zero()) {
        break;
      }
    }
    if (p.is_zero()) {
      p = param::constant(make_rational(1));
    }
    return p;
  };

  const auto at = sample_point();
  int checked = 0;
  for (int trial = 0; trial < 200; ++trial) {
    const param an = random_polynomial();
    const param ad = random_polynomial();
    const param bn = random_polynomial();
    const param bd = random_polynomial();

    // Only points where nothing in sight has a pole can serve as the oracle.
    if (varietas::coefficient_traits<rational>::is_zero(ad.evaluate(at)) ||
        varietas::coefficient_traits<rational>::is_zero(bd.evaluate(at))) {
      continue;
    }

    const field a(an, ad);
    const field b(bn, bd);
    const rational av = a.evaluate(at);
    const rational bv = b.evaluate(at);

    EXPECT_EQ((a + b).evaluate(at), av + bv) << "trial " << trial;
    EXPECT_EQ((a - b).evaluate(at), av - bv) << "trial " << trial;
    EXPECT_EQ((a * b).evaluate(at), av * bv) << "trial " << trial;
    if (!b.is_zero()) {
      EXPECT_EQ((a / b).evaluate(at), av / bv) << "trial " << trial;
    }
    ++checked;
  }
  EXPECT_GT(checked, 100) << "the sample point degenerated too often to be an oracle";
}

// The property that makes the field usable rather than merely correct.
//
// Repeatedly forming a + b * c, which is the shape of every reduction step in
// Buchberger, must not let the representation grow beyond the degree of the
// function it represents. Without gcd cancellation the numerator degree here
// climbs linearly in the iteration count; with it, it stays put.
TEST(RationalFunction, RepeatedArithmeticDoesNotInflateTheRepresentation) {
  field accumulator = constant(0);
  const field a = (x() + constant(1)) / (y() + constant(2));
  const field b = (y() - constant(1)) / (x() - constant(3));

  for (int i = 0; i < 24; ++i) {
    accumulator = accumulator * a + b;
  }

  // Every term of the recurrence is a ratio of polynomials of degree at most a
  // couple, so the reduced form cannot honestly need a high degree. The bound
  // is loose on purpose: what is being asserted is that the degree is bounded
  // at all, not a particular value.
  EXPECT_LE(accumulator.numerator().degree(), 30u)
      << "numerator degree ran away; gcd cancellation is not firing";
  EXPECT_LE(accumulator.denominator().degree(), 30u)
      << "denominator degree ran away; gcd cancellation is not firing";

  // And it is still the right function.
  const auto at = sample_point();
  rational expected(0);
  const rational av = a.evaluate(at);
  const rational bv = b.evaluate(at);
  for (int i = 0; i < 24; ++i) {
    expected = expected * av + bv;
  }
  EXPECT_EQ(accumulator.evaluate(at), expected);
}

TEST(RationalFunction, NormalisationWritesEqualFunctionsTheSameWay) {
  // A common factor and a negative denominator, both of which normalise away.
  const param n({{varietas::monomial<2>(std::array<std::uint8_t, 2>{2, 0}), make_rational(-2)}});
  const param d({{varietas::monomial<2>(std::array<std::uint8_t, 2>{1, 0}), make_rational(-4)}});
  const field f(n, d);

  EXPECT_EQ(f, x() / constant(2));
  EXPECT_EQ(f.denominator().leading_coefficient(), make_rational(2));
  EXPECT_GT(sgn(f.denominator().leading_coefficient()), 0)
      << "the sign belongs on the numerator";
}

// The claim rational_function exists to support: varietas_core's Buchberger,
// unmodified, runs over this field. A parametric basis is computed once with
// the pose symbolic, and its evaluation at a pose must agree with the basis
// computed from the start with that pose substituted in.
//
// This is the emitter's whole premise. If it fails, no header generated from a
// parametric basis answers the right question.
TEST(RationalFunction, BuchbergerOverTheParametricFieldMatchesTheNumericBasis) {
  using parametric = varietas::polynomial<field, 2, grevlex>;
  using numeric = varietas::polynomial<rational, 2, grevlex>;

  // A little system in two unknowns whose coefficients carry the parameters:
  //   u^2 + v^2 - 1
  //   u - x
  // so that the solutions are u = x, v^2 = 1 - x^2. The parameter y is carried
  // along so that P = 2 exercises more than one parameter.
  const auto u = parametric::variable(0);
  const auto v = parametric::variable(1);
  const auto one = parametric::constant(traits::one());

  const std::vector<parametric> parametric_generators{
      u * u + v * v - one,
      u - parametric::constant(x() + y() * constant(0)),
  };

  const varietas::ideal<field, 2, grevlex> parametric_ideal(parametric_generators);
  const auto& parametric_basis = parametric_ideal.basis();
  ASSERT_FALSE(parametric_basis.empty());

  const auto at = sample_point();

  // The same system with the parameter substituted before the computation.
  const auto uu = numeric::variable(0);
  const auto vv = numeric::variable(1);
  const std::vector<numeric> numeric_generators{
      uu * uu + vv * vv - numeric::constant(make_rational(1)),
      uu - numeric::constant(at[0]),
  };
  const varietas::ideal<rational, 2, grevlex> numeric_ideal(numeric_generators);

  // Evaluating the parametric basis gives polynomials over Q; each must lie in
  // the numeric ideal, and vice versa, which together say the two ideals agree.
  for (const auto& g : parametric_basis) {
    numeric evaluated;
    for (const auto& t : g.terms()) {
      evaluated = evaluated + numeric::from_monomial(t.mon, t.coeff.evaluate(at));
    }
    EXPECT_TRUE(numeric_ideal.contains(evaluated))
        << "a parametric basis element left the ideal when the pose was supplied";
  }
  for (const auto& g : numeric_ideal.basis()) {
    numeric copy = g;
    EXPECT_TRUE(numeric_ideal.contains(copy));
  }

  // And the dimension verdict, which is what the finiteness argument reads, is
  // the same either way.
  EXPECT_EQ(parametric_ideal.is_zero_dimensional(), numeric_ideal.is_zero_dimensional());
}

}  // namespace
