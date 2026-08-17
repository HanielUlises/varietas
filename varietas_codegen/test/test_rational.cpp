#include <cmath>

#include <gtest/gtest.h>

#include "varietas/codegen/rational.hpp"

namespace {

using varietas::make_rational;
using varietas::rational;
using varietas::rational_from_string;

using traits = varietas::coefficient_traits<rational>;

TEST(Rational, TraitsAdvertiseAnExactField) {
  EXPECT_TRUE(traits::is_exact);
  EXPECT_FALSE(varietas::coefficient_traits<double>::is_exact);
}

TEST(Rational, ValuesAreKeptInLowestTerms) {
  const rational q = make_rational(6, -8);
  EXPECT_EQ(q.get_num(), -3);
  EXPECT_EQ(q.get_den(), 4);

  // Canonical form makes equality structural, which is what the polynomial
  // layer relies on when it decides that a term has cancelled.
  EXPECT_EQ(make_rational(2, 4), make_rational(1, 2));
  EXPECT_EQ(make_rational(-1, 3), make_rational(3, -9));
}

TEST(Rational, ZeroIsDetectedAfterExactCancellation) {
  // Sevenths, not thirds: three thirds happen to round back to exactly 1.0 in
  // binary, which is a good reminder that a field is not exact merely because
  // one example came out right.
  rational sum = traits::zero();
  for (int i = 0; i < 7; ++i) {
    sum = sum + make_rational(1, 7);
  }
  sum = sum - traits::one();

  EXPECT_TRUE(traits::is_zero(sum));
  EXPECT_EQ(sum, 0);

  double approximate = 0.0;
  for (int i = 0; i < 7; ++i) {
    approximate += 1.0 / 7.0;
  }
  approximate -= 1.0;
  EXPECT_FALSE(varietas::coefficient_traits<double>::is_zero(approximate));
  EXPECT_NEAR(approximate, 0.0, 1e-15);
}

TEST(Rational, InverseAndNegateRoundTrip) {
  const rational q = make_rational(-7, 12);

  EXPECT_EQ(q * traits::inverse(q), traits::one());
  EXPECT_EQ(traits::inverse(traits::inverse(q)), q);
  EXPECT_EQ(q + traits::negate(q), traits::zero());
  EXPECT_EQ(traits::negate(traits::negate(q)), q);
}

TEST(Rational, CoefficientsGrowWithoutBound) {
  // A 64 bit numerator is exhausted long before Buchberger is; this is the
  // reason the field is not backed by a fixed width integer.
  rational q = make_rational(3, 7);
  for (int i = 0; i < 12; ++i) {
    q = q * q;
  }

  EXPECT_GT(mpz_sizeinbase(q.get_num().get_mpz_t(), 2), 64u);
  EXPECT_GT(mpz_sizeinbase(q.get_den().get_mpz_t(), 2), 64u);
}

TEST(Rational, ConversionToDoubleIsTheBridgeToEigen) {
  EXPECT_DOUBLE_EQ(traits::to_double(make_rational(1, 4)), 0.25);
  EXPECT_DOUBLE_EQ(traits::to_double(make_rational(-3, 2)), -1.5);
  EXPECT_NEAR(traits::to_double(make_rational(1, 3)), 1.0 / 3.0, 1e-15);
}

TEST(Rational, ConversionFromDoubleIsExactAndThereforeNotWhatYouMeant) {
  // from_double converts the binary value actually stored, so it is exact but
  // it is not 1/10. Anything reading decimal input must go through
  // rational_from_string or make_rational instead.
  const rational converted = traits::from_double(0.1);
  EXPECT_NE(converted, make_rational(1, 10));
  EXPECT_DOUBLE_EQ(traits::to_double(converted), 0.1);

  // Dyadic values do survive the round trip exactly.
  EXPECT_EQ(traits::from_double(0.25), make_rational(1, 4));
}

TEST(Rational, ParsesFractionsAndCanonicalisesThem) {
  EXPECT_EQ(rational_from_string("-7/2"), make_rational(-7, 2));
  EXPECT_EQ(rational_from_string("6/8"), make_rational(3, 4));
  EXPECT_EQ(rational_from_string("5"), make_rational(5));
}

}  // namespace
