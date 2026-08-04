#ifndef VARIETAS_CORE_POLYNOMIAL_HPP
#define VARIETAS_CORE_POLYNOMIAL_HPP

#include <algorithm>
#include <array>
#include <cstddef>
#include <initializer_list>
#include <vector>

#include "varietas/core/config.hpp"
#include "varietas/core/monomial.hpp"

namespace varietas {

// Sparse multivariate polynomial in N variables over Coeff, terms kept sorted
// strictly decreasing under Order and free of zero coefficients. The order is
// part of the type: polynomials carrying different orders do not mix.
//
// Coeff is a template parameter so that varietas_codegen can instantiate the
// same type over exact rationals. varietas_core only ever uses double.
template <class Coeff, std::size_t N, class Order>
class polynomial {
 public:
  using coefficient_type = Coeff;
  using monomial_type = monomial<N>;
  using order_type = Order;
  using traits = coefficient_traits<Coeff>;

  static constexpr std::size_t num_vars = N;

  struct term {
    monomial_type mon;
    Coeff coeff;
  };

  polynomial() = default;

  polynomial(std::initializer_list<term> terms) : terms_(terms) { normalize(); }

  explicit polynomial(std::vector<term> terms) : terms_(std::move(terms)) { normalize(); }

  static polynomial constant(const Coeff& c) {
    polynomial p;
    if (!traits::is_zero(c)) {
      p.terms_.push_back(term{monomial_type::one(), c});
    }
    return p;
  }

  static polynomial variable(std::size_t i, typename monomial_type::exponent_type power = 1) {
    polynomial p;
    p.terms_.push_back(term{monomial_type::variable(i, power), traits::one()});
    return p;
  }

  static polynomial from_monomial(const monomial_type& m, const Coeff& c) {
    polynomial p;
    if (!traits::is_zero(c)) {
      p.terms_.push_back(term{m, c});
    }
    return p;
  }

  bool is_zero() const noexcept { return terms_.empty(); }
  std::size_t size() const noexcept { return terms_.size(); }
  const std::vector<term>& terms() const noexcept { return terms_; }

  typename monomial_type::degree_type degree() const noexcept {
    typename monomial_type::degree_type d = 0;
    for (const term& t : terms_) {
      d = std::max(d, t.mon.degree());
    }
    return d;
  }

  const term& leading_term() const {
    VARIETAS_ASSERT(!terms_.empty());
    return terms_.front();
  }

  const monomial_type& leading_monomial() const { return leading_term().mon; }
  const Coeff& leading_coefficient() const { return leading_term().coeff; }

  Coeff coefficient_of(const monomial_type& m) const {
    const auto it = std::lower_bound(
        terms_.begin(), terms_.end(), m,
        [](const term& t, const monomial_type& key) { return Order::compare(t.mon, key) > 0; });
    if (it != terms_.end() && it->mon == m) {
      return it->coeff;
    }
    return traits::zero();
  }

  polynomial& operator+=(const polynomial& other) {
    terms_.insert(terms_.end(), other.terms_.begin(), other.terms_.end());
    normalize();
    return *this;
  }

  polynomial& operator-=(const polynomial& other) {
    terms_.reserve(terms_.size() + other.terms_.size());
    for (const term& t : other.terms_) {
      terms_.push_back(term{t.mon, -t.coeff});
    }
    normalize();
    return *this;
  }

  polynomial& operator*=(const Coeff& c) {
    if (traits::is_zero(c)) {
      terms_.clear();
      return *this;
    }
    for (term& t : terms_) {
      t.coeff = t.coeff * c;
    }
    return *this;
  }

  polynomial& operator*=(const polynomial& other) {
    *this = *this * other;
    return *this;
  }

  friend polynomial operator+(polynomial a, const polynomial& b) {
    a += b;
    return a;
  }

  friend polynomial operator-(polynomial a, const polynomial& b) {
    a -= b;
    return a;
  }

  friend polynomial operator-(const polynomial& a) {
    polynomial r;
    r.terms_.reserve(a.terms_.size());
    for (const term& t : a.terms_) {
      r.terms_.push_back(term{t.mon, -t.coeff});
    }
    return r;
  }

  friend polynomial operator*(const polynomial& a, const polynomial& b) {
    polynomial r;
    if (a.is_zero() || b.is_zero()) {
      return r;
    }
    r.terms_.reserve(a.terms_.size() * b.terms_.size());
    for (const term& s : a.terms_) {
      for (const term& t : b.terms_) {
        r.terms_.push_back(term{s.mon * t.mon, s.coeff * t.coeff});
      }
    }
    r.normalize();
    return r;
  }

  friend polynomial operator*(polynomial a, const Coeff& c) {
    a *= c;
    return a;
  }

  friend polynomial operator*(const Coeff& c, polynomial a) {
    a *= c;
    return a;
  }

  friend bool operator==(const polynomial& a, const polynomial& b) {
    if (a.terms_.size() != b.terms_.size()) {
      return false;
    }
    for (std::size_t i = 0; i < a.terms_.size(); ++i) {
      if (a.terms_[i].mon != b.terms_[i].mon || !(a.terms_[i].coeff == b.terms_[i].coeff)) {
        return false;
      }
    }
    return true;
  }

  friend bool operator!=(const polynomial& a, const polynomial& b) { return !(a == b); }

  polynomial multiply_by_monomial(const monomial_type& m, const Coeff& c) const {
    polynomial r;
    if (traits::is_zero(c)) {
      return r;
    }
    r.terms_.reserve(terms_.size());
    for (const term& t : terms_) {
      r.terms_.push_back(term{t.mon * m, t.coeff * c});
    }
    return r;
  }

  Coeff evaluate(const std::array<Coeff, N>& point) const {
    Coeff value = traits::zero();
    for (const term& t : terms_) {
      value = value + t.coeff * t.mon.template evaluate<Coeff>(point);
    }
    return value;
  }

  // Drops terms whose coefficient is negligible. Only meaningful for inexact
  // coefficients, and never applied implicitly by the arithmetic above.
  void prune(const Coeff& tolerance) {
    terms_.erase(std::remove_if(terms_.begin(), terms_.end(),
                                [&](const term& t) {
                                  return !(t.coeff > tolerance || t.coeff < -tolerance);
                                }),
                 terms_.end());
  }

 private:
  void normalize() {
    std::stable_sort(terms_.begin(), terms_.end(), [](const term& a, const term& b) {
      return Order::compare(a.mon, b.mon) > 0;
    });

    std::size_t out = 0;
    for (std::size_t in = 0; in < terms_.size();) {
      Coeff accumulated = terms_[in].coeff;
      const monomial_type mon = terms_[in].mon;
      std::size_t next = in + 1;
      while (next < terms_.size() && terms_[next].mon == mon) {
        accumulated = accumulated + terms_[next].coeff;
        ++next;
      }
      if (!traits::is_zero(accumulated)) {
        terms_[out].mon = mon;
        terms_[out].coeff = accumulated;
        ++out;
      }
      in = next;
    }
    terms_.resize(out);
  }

  std::vector<term> terms_;
};

}  // namespace varietas

#endif
