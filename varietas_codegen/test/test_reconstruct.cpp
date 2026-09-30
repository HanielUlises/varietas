// Rational functions recovered from their values modulo primes.
//
// The black box here is the functions themselves, evaluated at points of Z/p,
// so the answer is known and the test is equality over Q. The functions are
// chosen for the ways reconstruction can go wrong: a numerator and denominator
// of different degrees, a function that is identically zero, a constant, a
// polynomial with no denominator, coefficients with denominators larger than
// one prime, and a denominator whose leading coefficient is not one.

#include <array>
#include <cstddef>
#include <vector>

#include <gtest/gtest.h>

#include "varietas/codegen/prime_field.hpp"
#include "varietas/codegen/rational.hpp"
#include "varietas/codegen/rational_function.hpp"
#include "varietas/codegen/reconstruct.hpp"

namespace {

using varietas::make_rational;
using varietas::rational;
using varietas::rational_from_string;
using varietas::residue;
using field = varietas::rational_function<3>;

field x() { return field::parameter(0); }
field y() { return field::parameter(1); }
field z() { return field::parameter(2); }
field c(long n, long d = 1) { return field(make_rational(n, d)); }

// A rational function's value at a point of Z/p, through its coefficients.
residue value_at(const field& f, const varietas::reconstruction::point<3>& at) {
  const auto evaluate = [&](const field::parameter_polynomial& p) {
    residue sum;
    for (const auto& t : p.terms()) {
      residue term(t.coeff);
      for (std::size_t v = 0; v < 3; ++v) {
        for (unsigned k = 0; k < t.mon[v]; ++k) {
          term = term * at[v];
        }
      }
      sum = sum + term;
    }
    return sum;
  };
  return evaluate(f.numerator()) / evaluate(f.denominator());
}

TEST(Reconstruct, KnownFunctionsAreRecoveredExactly) {
  const rational big = rational_from_string("98765432123456789/1234567");
  const std::vector<field> wanted = {
      (x() * x() + c(3) * y() - z() / c(2)) / (c(1) + x() * y()),
      c(0),
      c(5, 7),
      x() * y() * z() - c(2) * x() + field(big),
      c(1) / (c(3) * x() * x() + z() + c(11, 13)),
      (x() - y()) * (x() + z()) / (c(7) * z() * z() * z() + c(2) * x() * y() - c(1, 5)),
  };

  const varietas::reconstruction::prime_hook prepare = [](varietas::reconstruction::word) {
    return true;
  };
  const varietas::reconstruction::black_box<3> box =
      [&](const varietas::reconstruction::point<3>& at, std::vector<residue>& values) {
        values.clear();
        for (const auto& f : wanted) {
          // A pole of any function makes the point unusable, as a special
          // fibre would in the kinematic black box.
          const auto d = f.denominator();
          residue den;
          for (const auto& t : d.terms()) {
            residue term(t.coeff);
            for (std::size_t v = 0; v < 3; ++v) {
              for (unsigned k = 0; k < t.mon[v]; ++k) {
                term = term * at[v];
              }
            }
            den = den + term;
          }
          if (den.value() == 0) {
            return false;
          }
          values.push_back(value_at(f, at));
        }
        return true;
      };

  std::vector<field> got;
  varietas::reconstruction::statistics stats;
  ASSERT_TRUE(varietas::reconstruction::reconstruct<3>(box, prepare, wanted.size(), got, &stats));
  ASSERT_EQ(got.size(), wanted.size());
  for (std::size_t f = 0; f < wanted.size(); ++f) {
    EXPECT_EQ(got[f], wanted[f]) << "function " << f;
  }
  EXPECT_EQ(stats.largest_numerator_degree, 3);  // x y z
  EXPECT_EQ(stats.largest_denominator_degree, 3);
  // The large constant needs more than one prime to reconstruct.
  EXPECT_GE(stats.primes, 3u);
}

// The same inputs always give the same answer and the same count of samples,
// since a generated header must not depend on when it was generated.
TEST(Reconstruct, TheComputationIsDeterministic) {
  const std::vector<field> wanted = {(x() + c(2) * y()) / (z() + c(3))};
  const varietas::reconstruction::prime_hook prepare = [](varietas::reconstruction::word) {
    return true;
  };
  const varietas::reconstruction::black_box<3> box =
      [&](const varietas::reconstruction::point<3>& at, std::vector<residue>& values) {
        if ((at[2] + residue(3)).value() == 0) {
          return false;
        }
        values = {value_at(wanted[0], at)};
        return true;
      };
  std::vector<field> first;
  std::vector<field> second;
  varietas::reconstruction::statistics a;
  varietas::reconstruction::statistics b;
  ASSERT_TRUE(varietas::reconstruction::reconstruct<3>(box, prepare, 1, first, &a));
  ASSERT_TRUE(varietas::reconstruction::reconstruct<3>(box, prepare, 1, second, &b));
  EXPECT_EQ(first[0], second[0]);
  EXPECT_EQ(a.samples, b.samples);
}

TEST(PrimeField, ArithmeticAndTheImageOfARational) {
  residue::select(1000003);
  const residue a(make_rational(1, 3));
  EXPECT_EQ((a * residue(3)).value(), 1u);
  EXPECT_EQ((residue(2) - residue(5)).value(), 1000003u - 3u);
  EXPECT_EQ((residue(7) / residue(7)).value(), 1u);

  residue::clear_undefined_image();
  const residue undefined(make_rational(1, 1000003));
  (void)undefined;
  EXPECT_TRUE(residue::undefined_image()) << "a denominator the prime divides has no image";
  residue::clear_undefined_image();
  residue::select(2147483647);
}

}  // namespace
