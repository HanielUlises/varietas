#ifndef VARIETAS_CORE_MONOMIAL_HPP
#define VARIETAS_CORE_MONOMIAL_HPP

#include <array>
#include <cstddef>
#include <cstdint>

#include "varietas/core/config.hpp"

namespace varietas {

// Exponent vector of a monomial in N variables. The total degree is cached
// because every graded order needs it on the comparison path.
template <std::size_t N>
class monomial {
 public:
  using exponent_type = std::uint8_t;
  using degree_type = std::uint16_t;

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

  friend constexpr monomial operator*(const monomial& a, const monomial& b) noexcept {
    monomial r;
    for (std::size_t i = 0; i < N; ++i) {
      r.exponents_[i] = static_cast<exponent_type>(a.exponents_[i] + b.exponents_[i]);
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
