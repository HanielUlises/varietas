// The polynomial layer, exercised over both coefficient fields. Every
// assertion here is a statement about the ring k[x, y, z] that holds for any
// field k; where the two instantiations diverge is the subject of
// test_exact_groebner.

#include <array>
#include <vector>

#include <gtest/gtest.h>

#include "varietas/codegen/rational.hpp"
#include "varietas/core/ideal/division.hpp"
#include "varietas/core/order/grevlex.hpp"
#include "varietas/core/order/lex.hpp"
#include "varietas/core/polynomial.hpp"

#include "coefficient_fixture.hpp"

namespace {

using varietas::grevlex;
using varietas::lex;
using varietas_test::coeff;

template <class Coeff>
class PolynomialOverField : public ::testing::Test {
 public:
  using poly = varietas::polynomial<Coeff, 3, lex>;
  using traits = varietas::coefficient_traits<Coeff>;

  static Coeff c(std::int64_t n, std::int64_t d = 1) { return coeff<Coeff>(n, d); }

  static poly x() { return poly::variable(0); }
  static poly y() { return poly::variable(1); }
  static poly z() { return poly::variable(2); }
  static poly constant(std::int64_t n, std::int64_t d = 1) { return poly::constant(c(n, d)); }
};

using field_types = ::testing::Types<double, varietas::rational>;
TYPED_TEST_SUITE(PolynomialOverField, field_types);

TYPED_TEST(PolynomialOverField, TermsAreSortedAndZeroFree) {
  using poly = typename TestFixture::poly;

  const poly f = this->x() * this->y() + this->z() + this->x() * this->x();

  ASSERT_EQ(f.size(), 3u);
  for (std::size_t i = 1; i < f.size(); ++i) {
    EXPECT_GT(lex::compare(f.terms()[i - 1].mon, f.terms()[i].mon), 0);
  }
  for (const auto& t : f.terms()) {
    EXPECT_FALSE(TestFixture::traits::is_zero(t.coeff));
  }
}

TYPED_TEST(PolynomialOverField, AdditiveInverseCancels) {
  const auto f = this->x() * this->y() - this->constant(3) * this->z() + this->constant(1, 2);
  EXPECT_TRUE((f - f).is_zero());
  EXPECT_TRUE((f + (-f)).is_zero());
}

TYPED_TEST(PolynomialOverField, MultiplicationDistributes) {
  const auto f = this->x() + this->constant(1, 3);
  const auto g = this->y() - this->constant(2, 5);
  const auto h = this->z() * this->z() + this->constant(7);

  EXPECT_EQ(f * (g + h), f * g + f * h);
}

TYPED_TEST(PolynomialOverField, MonicRescalingPreservesTheZeroSet) {
  const auto f = this->constant(6) * this->x() * this->x() - this->constant(4) * this->y();
  const auto m = f.monic();

  EXPECT_EQ(m.leading_coefficient(), TestFixture::traits::one());
  EXPECT_EQ(m * f.leading_coefficient(), f);
}

TYPED_TEST(PolynomialOverField, DivisionReconstructsTheDividend) {
  using poly = typename TestFixture::poly;

  const poly f = this->x() * this->x() * this->y() + this->x() * this->y() * this->y() +
                 this->y() * this->y();
  const std::vector<poly> divisors{this->x() * this->y() - this->constant(1),
                                   this->y() * this->y() - this->constant(1)};

  const auto result = varietas::divide(f, divisors);

  poly reconstructed = result.remainder;
  for (std::size_t i = 0; i < divisors.size(); ++i) {
    reconstructed += result.quotients[i] * divisors[i];
  }
  EXPECT_EQ(reconstructed, f);
}

TYPED_TEST(PolynomialOverField, EvaluationAgreesWithSubstitution) {
  using coeff_type = TypeParam;

  const auto f = this->x() * this->y() + this->z() * this->z() - this->constant(1);
  const std::array<coeff_type, 3> point{this->c(2), this->c(3), this->c(4)};

  // 2*3 + 4*4 - 1 = 21
  EXPECT_EQ(f.evaluate(point), this->c(21));
}

}  // namespace
