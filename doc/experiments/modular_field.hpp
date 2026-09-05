#ifndef VARIETAS_EXPERIMENTS_MODULAR_FIELD_HPP
#define VARIETAS_EXPERIMENTS_MODULAR_FIELD_HPP

#include <cstddef>
#include <cstdint>
#include <chrono>
#include <cstdio>

#include "varietas/core/config.hpp"
#include "varietas/core/gcd.hpp"
#include "varietas/core/ideal/division.hpp"
#include "varietas/core/order/grevlex.hpp"
#include "varietas/core/polynomial.hpp"

// A prime field and rational functions over it.
//
// This exists to isolate one variable in the cost of a parametric Grobner basis
// computation. The pipeline computes over Q(p), where every coefficient
// operation is arbitrary-precision rational arithmetic and every normalisation
// runs a subresultant gcd over those rationals. Replacing Q by F_p leaves the
// algorithm, the monomial order, the ideal and the number of parameters
// unchanged, and reduces each coefficient operation to one machine-word
// multiplication and one remainder.
//
// The comparison therefore separates two explanations for the observed cost:
// the arithmetic on rational coefficients, and the number and size of the
// polynomials the algorithm produces. The second is unaffected by the change of
// field for a prime that is not unlucky for the input.
namespace varietas {

// F_p for a prime that fits in 31 bits, so that products fit in 63.
class modular {
 public:
  static constexpr long long prime = 2147483647;  // 2^31 - 1

  constexpr modular() noexcept : value_(0) {}
  constexpr modular(long long v) noexcept : value_(((v % prime) + prime) % prime) {}

  constexpr long long value() const noexcept { return value_; }

  friend constexpr modular operator+(modular a, modular b) noexcept {
    return modular::raw((a.value_ + b.value_) % prime);
  }
  friend constexpr modular operator-(modular a, modular b) noexcept {
    return modular::raw((a.value_ - b.value_ + prime) % prime);
  }
  friend constexpr modular operator*(modular a, modular b) noexcept {
    return modular::raw(a.value_ * b.value_ % prime);
  }
  friend modular operator/(modular a, modular b) { return a * inverted(b); }
  friend constexpr modular operator-(modular a) noexcept {
    return modular::raw((prime - a.value_) % prime);
  }

  modular& operator+=(modular o) { return *this = *this + o; }
  modular& operator-=(modular o) { return *this = *this - o; }
  modular& operator*=(modular o) { return *this = *this * o; }
  modular& operator/=(modular o) { return *this = *this / o; }

  friend constexpr bool operator==(modular a, modular b) noexcept {
    return a.value_ == b.value_;
  }
  friend constexpr bool operator!=(modular a, modular b) noexcept {
    return a.value_ != b.value_;
  }

  // Fermat's little theorem. The field is small enough that a table would not
  // help and the exponentiation is not on any hot path that matters here.
  static modular inverted(modular a) {
    VARIETAS_ASSERT(a.value_ != 0);
    long long result = 1;
    long long base = a.value_;
    long long exponent = prime - 2;
    while (exponent > 0) {
      if ((exponent & 1) != 0) {
        result = result * base % prime;
      }
      base = base * base % prime;
      exponent >>= 1;
    }
    return raw(result);
  }

 private:
  static constexpr modular raw(long long v) noexcept {
    modular m;
    m.value_ = v;
    return m;
  }
  long long value_;
};

template <>
struct coefficient_traits<modular> {
  static constexpr bool is_exact = true;

  static modular zero() { return modular(0); }
  static modular one() { return modular(1); }
  static bool is_zero(modular c) { return c == modular(0); }
  static modular inverse(modular c) { return modular::inverted(c); }
  static modular negate(modular c) { return -c; }
  static double to_double(modular c) { return static_cast<double>(c.value()); }
  static modular from_double(double c) { return modular(static_cast<long long>(c)); }
};

// F_p(p_1, ..., p_P), written the same way rational_function writes Q(p).
//
// The normalisation is the same in structure: divide out the monomial content,
// divide out the polynomial gcd, and fix a representative by making the
// denominator monic. Only the coefficient arithmetic underneath differs.
template <std::size_t P>
class modular_function {
 public:
  using parameter_order = grevlex;
  using parameter_polynomial = polynomial<modular, P, parameter_order>;

  static constexpr std::size_t num_parameters = P;

  modular_function() : numerator_(), denominator_(parameter_polynomial::constant(modular(1))) {}

  // NOLINTNEXTLINE(google-explicit-constructor) — a scalar is a constant function
  modular_function(modular c)
      : numerator_(parameter_polynomial::constant(c)),
        denominator_(parameter_polynomial::constant(modular(1))) {
    normalize();
  }

  modular_function(const parameter_polynomial& n, const parameter_polynomial& d)
      : numerator_(n), denominator_(d) {
    VARIETAS_ASSERT(!d.is_zero());
    normalize();
  }

  static modular_function parameter(std::size_t i) {
    return modular_function(parameter_polynomial::variable(i),
                            parameter_polynomial::constant(modular(1)));
  }

  const parameter_polynomial& numerator() const noexcept { return numerator_; }
  const parameter_polynomial& denominator() const noexcept { return denominator_; }
  bool is_zero() const noexcept { return numerator_.is_zero(); }

  modular_function& operator+=(const modular_function& o) {
    numerator_ = numerator_ * o.denominator_ + o.numerator_ * denominator_;
    denominator_ = denominator_ * o.denominator_;
    normalize();
    return *this;
  }
  modular_function& operator-=(const modular_function& o) {
    numerator_ = numerator_ * o.denominator_ - o.numerator_ * denominator_;
    denominator_ = denominator_ * o.denominator_;
    normalize();
    return *this;
  }
  modular_function& operator*=(const modular_function& o) {
    numerator_ = numerator_ * o.numerator_;
    denominator_ = denominator_ * o.denominator_;
    normalize();
    return *this;
  }
  modular_function& operator/=(const modular_function& o) {
    VARIETAS_ASSERT(!o.is_zero());
    numerator_ = numerator_ * o.denominator_;
    denominator_ = denominator_ * o.numerator_;
    normalize();
    return *this;
  }

  friend modular_function operator+(modular_function a, const modular_function& b) { return a += b; }
  friend modular_function operator-(modular_function a, const modular_function& b) { return a -= b; }
  friend modular_function operator*(modular_function a, const modular_function& b) { return a *= b; }
  friend modular_function operator/(modular_function a, const modular_function& b) { return a /= b; }
  friend modular_function operator-(const modular_function& a) {
    return modular_function(-a.numerator_, a.denominator_);
  }
  friend bool operator==(const modular_function& a, const modular_function& b) {
    return a.numerator_ * b.denominator_ == b.numerator_ * a.denominator_;
  }
  friend bool operator!=(const modular_function& a, const modular_function& b) { return !(a == b); }

 private:
  // Instrumentation: the largest numerator seen, reported periodically.
  //
  // Over F_p a coefficient operation is a machine multiply, so anything the run
  // spends its time on is the number of terms in these polynomials rather than
  // the size of the numbers in them. Recording the maximum makes that visible.
  static void note_size(std::size_t terms) {
#ifndef VARIETAS_EXPERIMENT_TRACE
    (void)terms;
#else
    static std::size_t largest = 0;
    static auto started = std::chrono::steady_clock::now();
    static double next = 5.0;
    if (terms > largest) {
      largest = terms;
    }
    const double elapsed =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - started).count();
    if (elapsed > next) {
      next += 15.0;
      std::fprintf(stderr, "%7.0f s   largest parameter polynomial %8zu terms\n", elapsed,
                   largest);
      std::fflush(stderr);
    }
#endif
  }

  void normalize() {
    note_size(numerator_.size() > denominator_.size() ? numerator_.size() : denominator_.size());
    if (numerator_.is_zero()) {
      denominator_ = parameter_polynomial::constant(modular(1));
      return;
    }
    remove_monomial_content(numerator_, denominator_);
    remove_common_factor(numerator_, denominator_);

    const modular lead = denominator_.leading_coefficient();
    if (lead != modular(1)) {
      const modular scale = modular::inverted(lead);
      numerator_ = numerator_ * parameter_polynomial::constant(scale);
      denominator_ = denominator_ * parameter_polynomial::constant(scale);
    }
  }

  static void remove_common_factor(parameter_polynomial& n, parameter_polynomial& d) {
    if (n.degree() == 0 || d.degree() == 0) {
      return;
    }
    const parameter_polynomial g = polynomial_gcd(n, d);
    if (g.degree() == 0) {
      return;
    }
    n = exact_quotient(n, g);
    d = exact_quotient(d, g);
  }

  static parameter_polynomial exact_quotient(const parameter_polynomial& p,
                                             const parameter_polynomial& g) {
    const auto result = divide(p, std::vector<parameter_polynomial>{g});
    VARIETAS_ASSERT(result.remainder.is_zero());
    return result.quotients.front();
  }

  static void remove_monomial_content(parameter_polynomial& n, parameter_polynomial& d) {
    monomial<P> common = n.terms().front().mon;
    for (const auto& t : n.terms()) {
      common = monomial<P>::gcd(common, t.mon);
    }
    for (const auto& t : d.terms()) {
      common = monomial<P>::gcd(common, t.mon);
    }
    if (common.is_one()) {
      return;
    }
    n = divide_by_monomial(n, common);
    d = divide_by_monomial(d, common);
  }

  static parameter_polynomial divide_by_monomial(const parameter_polynomial& p,
                                                 const monomial<P>& m) {
    std::vector<typename parameter_polynomial::term> terms;
    terms.reserve(p.size());
    for (const auto& t : p.terms()) {
      terms.push_back({monomial<P>::divide(t.mon, m), t.coeff});
    }
    return parameter_polynomial(std::move(terms));
  }

  parameter_polynomial numerator_;
  parameter_polynomial denominator_;
};

template <std::size_t P>
struct coefficient_traits<modular_function<P>> {
  using value_type = modular_function<P>;

  static constexpr bool is_exact = true;

  static value_type zero() { return value_type(modular(0)); }
  static value_type one() { return value_type(modular(1)); }
  static bool is_zero(const value_type& c) { return c.is_zero(); }
  static value_type inverse(const value_type& c) {
    VARIETAS_ASSERT(!c.is_zero());
    return value_type(c.denominator(), c.numerator());
  }
  static value_type negate(const value_type& c) { return -c; }

  // A rational function has no value until the parameters are supplied, so the
  // bridge to floating point is not available and is not needed here.
  static double to_double(const value_type&) {
    VARIETAS_ASSERT(false);
    return 0.0;
  }
  static value_type from_double(double c) {
    return value_type(modular(static_cast<long long>(c)));
  }
};

}  // namespace varietas

#endif
