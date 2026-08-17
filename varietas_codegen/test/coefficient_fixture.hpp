#ifndef VARIETAS_CODEGEN_TEST_COEFFICIENT_FIXTURE_HPP
#define VARIETAS_CODEGEN_TEST_COEFFICIENT_FIXTURE_HPP

#include <array>
#include <cstdint>

#include "varietas/codegen/rational.hpp"
#include "varietas/core/monomial.hpp"

namespace varietas_test {

// Builds a coefficient from an integer fraction in whichever field is under
// test, so that one test body can be instantiated over both. Note the
// asymmetry the suites are there to expose: over rational the value is n/d,
// over double it is the nearest representable approximation to it.
template <class Coeff>
struct coefficient_factory;

template <>
struct coefficient_factory<double> {
  static double of(std::int64_t numerator, std::int64_t denominator = 1) {
    return static_cast<double>(numerator) / static_cast<double>(denominator);
  }
};

template <>
struct coefficient_factory<varietas::rational> {
  static varietas::rational of(std::int64_t numerator, std::int64_t denominator = 1) {
    return varietas::make_rational(numerator, denominator);
  }
};

template <class Coeff>
Coeff coeff(std::int64_t numerator, std::int64_t denominator = 1) {
  return coefficient_factory<Coeff>::of(numerator, denominator);
}

template <std::size_t N>
varietas::monomial<N> mon(const std::array<std::uint8_t, N>& exponents) {
  return varietas::monomial<N>(exponents);
}

}  // namespace varietas_test

#endif
