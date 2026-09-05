#ifndef VARIETAS_CORE_MONOMIAL_HPP
#define VARIETAS_CORE_MONOMIAL_HPP

#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>

#include "varietas/core/config.hpp"

namespace varietas {

// Exponent vector of a monomial in N variables. The total degree is cached
// because every graded order needs it on the comparison path.
//
// The exponent is sixteen bits rather than eight, and the reason is not
// generality for its own sake. A Gröbner basis computed over a function field
// carries polynomials in the parameters as its coefficients, and those grow far
// past anything the ring variables reach: an underdetermined system over
// Q(x, y) was observed to produce parameter polynomials of degree 254 in a
// single variable, one step from wrapping an eight-bit exponent around to zero.
// What made that dangerous was not the ceiling but the silence: a wrapped
// exponent is a different monomial, so a polynomial that divided another stops
// dividing it, and the failure surfaces somewhere else entirely as an exact
// division that leaves a remainder.
//
// A wider exponent moves the ceiling out by a factor of two hundred and
// fifty-six; the assertion in operator* below is what makes reaching it loud
// rather than silent, and that is the part that matters. Neither costs anything
// measurable next to the arbitrary-precision coefficients these monomials are
// attached to.
template <std::size_t N>
class monomial {
 public:
  using exponent_type = std::uint16_t;
  using degree_type = std::uint32_t;

  static constexpr std::size_t num_vars = N;

  constexpr monomial() noexcept : exponents_{}, degree_(0) {}

  explicit constexpr monomial(const std::array<exponent_type, N>& exponents) noexcept
      : exponents_(exponents), degree_(0) {
    for (std::size_t i = 0; i < N; ++i) {
      degree_ = static_cast<degree_type>(degree_ + exponents_[i]);
    }
  }

  static constexpr monomial one() noexcept { return monomial(); }

  static constexpr monomial variable(std::size_t i, exponent_type power = 1) noexcept {
    VARIETAS_ASSERT(i < N);
    monomial m;
    m.exponents_[i] = power;
    m.degree_ = power;
    return m;
  }

  constexpr exponent_type operator[](std::size_t i) const noexcept {
    VARIETAS_ASSERT(i < N);
    return exponents_[i];
  }

  constexpr degree_type degree() const noexcept { return degree_; }

  constexpr bool is_one() const noexcept { return degree_ == 0; }

  constexpr const std::array<exponent_type, N>& exponents() const noexcept {
    return exponents_;
  }

  // The only operation that can make an exponent grow, and so the only one that
  // can overflow. lcm and gcd take a maximum and a minimum of exponents already
  // present, and divide subtracts under a precondition that it is exact.
  friend constexpr monomial operator*(const monomial& a, const monomial& b) noexcept {
    monomial r;
    for (std::size_t i = 0; i < N; ++i) {
      const unsigned sum =
          static_cast<unsigned>(a.exponents_[i]) + static_cast<unsigned>(b.exponents_[i]);
      // Wrapping here would not throw anything away loudly: it would produce a
      // perfectly well formed monomial standing for a different one, and every
      // divisibility question asked of it afterwards would be answered wrongly.
      VARIETAS_ASSERT(sum <= static_cast<unsigned>(std::numeric_limits<exponent_type>::max()));
      r.exponents_[i] = static_cast<exponent_type>(sum);
    }
    r.degree_ = static_cast<degree_type>(a.degree_ + b.degree_);
    return r;
  }

  friend constexpr bool operator==(const monomial& a, const monomial& b) noexcept {
    if (a.degree_ != b.degree_) {
      return false;
    }
    for (std::size_t i = 0; i < N; ++i) {
      if (a.exponents_[i] != b.exponents_[i]) {
        return false;
      }
    }
    return true;
  }

  friend constexpr bool operator!=(const monomial& a, const monomial& b) noexcept {
    return !(a == b);
  }

  static constexpr bool divides(const monomial& a, const monomial& b) noexcept {
    if (a.degree_ > b.degree_) {
      return false;
    }
    for (std::size_t i = 0; i < N; ++i) {
      if (a.exponents_[i] > b.exponents_[i]) {
        return false;
      }
    }
    return true;
  }

  // Precondition: divides(a, b).
  static constexpr monomial divide(const monomial& b, const monomial& a) noexcept {
    VARIETAS_ASSERT(divides(a, b));
    monomial r;
    for (std::size_t i = 0; i < N; ++i) {
      r.exponents_[i] = static_cast<exponent_type>(b.exponents_[i] - a.exponents_[i]);
    }
    r.degree_ = static_cast<degree_type>(b.degree_ - a.degree_);
    return r;
  }

  static constexpr monomial lcm(const monomial& a, const monomial& b) noexcept {
    monomial r;
    for (std::size_t i = 0; i < N; ++i) {
      const exponent_type e = a.exponents_[i] > b.exponents_[i] ? a.exponents_[i]
                                                                : b.exponents_[i];
      r.exponents_[i] = e;
      r.degree_ = static_cast<degree_type>(r.degree_ + e);
    }
    return r;
  }

  static constexpr monomial gcd(const monomial& a, const monomial& b) noexcept {
    monomial r;
    for (std::size_t i = 0; i < N; ++i) {
      const exponent_type e = a.exponents_[i] < b.exponents_[i] ? a.exponents_[i]
                                                                : b.exponents_[i];
      r.exponents_[i] = e;
      r.degree_ = static_cast<degree_type>(r.degree_ + e);
    }
    return r;
  }

  template <class Coeff>
  constexpr Coeff evaluate(const std::array<Coeff, N>& point) const {
    Coeff value = coefficient_traits<Coeff>::one();
    for (std::size_t i = 0; i < N; ++i) {
      for (exponent_type k = 0; k < exponents_[i]; ++k) {
        value = value * point[i];
      }
    }
    return value;
  }

 private:
  std::array<exponent_type, N> exponents_;
  degree_type degree_;
};

}  // namespace varietas

#endif
