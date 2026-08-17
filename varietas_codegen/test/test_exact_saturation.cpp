// Saturation over both coefficient fields.
//
// Saturating is a Gröbner basis computation in one more variable than the
// problem has, with the generator 1 - y·h whose whole purpose is to cancel
// against the others. It is therefore exactly the kind of computation that
// depends on reduction to zero being decidable, and the field it runs over is
// not a matter of taste: over the exact field the answer is the ideal, and over
// double it is the ideal up to whatever the rounding left behind.

#include <cstdint>
#include <vector>

#include <gtest/gtest.h>

#include "varietas/codegen/rational.hpp"
#include "varietas/core/ideal/ideal.hpp"
#include "varietas/core/ideal/saturation.hpp"
#include "varietas/core/order/grevlex.hpp"
#include "varietas/core/polynomial.hpp"

#include "coefficient_fixture.hpp"

namespace {

using varietas::grevlex;
using varietas_test::coeff;

template <class Coeff>
class SaturationOverField : public ::testing::Test {
 public:
  using poly = varietas::polynomial<Coeff, 2, grevlex>;
  using traits = varietas::coefficient_traits<Coeff>;

  static Coeff c(std::int64_t n, std::int64_t d = 1) { return coeff<Coeff>(n, d); }

  static poly t0() { return poly::variable(0); }
  static poly t1() { return poly::variable(1); }
  static poly constant(std::int64_t n, std::int64_t d = 1) { return poly::constant(c(n, d)); }

  // 1 + t0^2, the denominator the half-angle substitution clears at the first
  // joint.
  static poly denominator() {
    return varietas::half_angle_denominator<Coeff, 2, grevlex>(0);
  }
};

using field_types = ::testing::Types<double, varietas::rational>;
TYPED_TEST_SUITE(SaturationOverField, field_types);

// The statement that holds over any field: the spurious component introduced by
// clearing the denominator is removed, and the configuration that survives is
// the one the arm has. The coefficients here are non-dyadic on purpose — a
// third and two sevenths are the sort of thing a link length produces once a
// pose has been substituted in.
TYPED_TEST(SaturationOverField, removes_the_denominators_component) {
  using Coeff = TypeParam;
  using fixture = SaturationOverField<Coeff>;

  const auto t0 = fixture::t0();
  const auto t1 = fixture::t1();
  const auto h = fixture::denominator();

  const auto first = h * (fixture::constant(3) * t0 - fixture::constant(1));
  const auto second = h * (fixture::constant(7) * t1 - fixture::constant(2));

  const auto saturated = varietas::saturate<Coeff, 2, grevlex>({first, second}, h);

  const varietas::ideal<Coeff, 2, grevlex> cleaned(saturated);
  EXPECT_TRUE(cleaned.is_zero_dimensional());
  EXPECT_EQ(cleaned.quotient().monomials.size(), 1u)
      << "one configuration, not three";

  // Two generators, each monic and linear in one variable.
  ASSERT_EQ(saturated.size(), 2u);
  for (const auto& g : saturated) {
    EXPECT_EQ(g.degree(), 1);
  }
}

// Where the fields part company. Over the rationals the saturated basis is
// exactly {t0 - 1/3, t1 - 2/7}, and the test can say so with equality; over
// double the same two polynomials carry the nearest representable coefficients
// instead, and only a tolerance can be asserted. The distinction is not
// cosmetic — the exact basis is what a generated header stores, and what its
// order_id and coefficients are later compared against.
TEST(exact_saturation, recovers_the_coefficients_exactly_over_the_rationals) {
  using Coeff = varietas::rational;
  using poly = varietas::polynomial<Coeff, 2, grevlex>;

  const poly t0 = poly::variable(0);
  const poly t1 = poly::variable(1);
  const poly h = varietas::half_angle_denominator<Coeff, 2, grevlex>(0);

  const poly first = h * (poly::constant(varietas::make_rational(3)) * t0 -
                          poly::constant(varietas::make_rational(1)));
  const poly second = h * (poly::constant(varietas::make_rational(7)) * t1 -
                           poly::constant(varietas::make_rational(2)));

  const auto saturated = varietas::saturate<Coeff, 2, grevlex>({first, second}, h);

  const varietas::ideal<Coeff, 2, grevlex> expected(
      {t0 - poly::constant(varietas::make_rational(1, 3)),
       t1 - poly::constant(varietas::make_rational(2, 7))});
  EXPECT_EQ(saturated, expected.basis());
}

TEST(exact_saturation, over_double_the_same_basis_is_only_approximate) {
  using poly = varietas::polynomial<double, 2, grevlex>;

  const poly t0 = poly::variable(0);
  const poly t1 = poly::variable(1);
  const poly h = varietas::half_angle_denominator<double, 2, grevlex>(0);

  const poly first = h * (poly::constant(3.0) * t0 - poly::constant(1.0));
  const poly second = h * (poly::constant(7.0) * t1 - poly::constant(2.0));

  const auto saturated = varietas::saturate<double, 2, grevlex>({first, second}, h);

  const varietas::ideal<double, 2, grevlex> expected(
      {t0 - poly::constant(1.0 / 3.0), t1 - poly::constant(2.0 / 7.0)});

  // Structurally the same basis — the shape of the answer survives rounding
  // here — but the constants are the nearest doubles to a third and two
  // sevenths, and the ideal they generate is not the one asked for.
  ASSERT_EQ(saturated.size(), expected.basis().size());
  for (std::size_t i = 0; i < saturated.size(); ++i) {
    EXPECT_EQ(saturated[i].leading_monomial(), expected.basis()[i].leading_monomial());
  }

  const double recovered = -saturated.back().coefficient_of(varietas::monomial<2>::one());
  EXPECT_NEAR(recovered, 2.0 / 7.0, 1e-15);

  // The exact field is what witnesses the difference. Read the double back as
  // the rational it actually is, and it is not two sevenths: the basis over
  // double generates an ideal whose variety is a point close to the one asked
  // for, and closeness is not what the completeness certificate asserts.
  EXPECT_NE(varietas::rational(recovered), varietas::make_rational(2, 7));
}

}  // namespace
