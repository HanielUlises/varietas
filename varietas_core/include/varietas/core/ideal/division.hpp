#ifndef VARIETAS_CORE_IDEAL_DIVISION_HPP
#define VARIETAS_CORE_IDEAL_DIVISION_HPP

#include <cstddef>
#include <vector>

#include "varietas/core/config.hpp"
#include "varietas/core/monomial.hpp"
#include "varietas/core/polynomial.hpp"

namespace varietas {

// Output of the multivariate division algorithm: f = sum_i q_i g_i + r, where
// no monomial of r is divisible by any leading monomial of the divisors.
template <class Poly>
struct division_result {
  std::vector<Poly> quotients;
  Poly remainder;
};

namespace detail {

// Index of the first divisor whose leading monomial divides m, or the size of
// the family if there is none. First-fit is the classical selection rule and
// is what makes the quotients depend on the order of the divisors.
template <class Poly, std::size_t N>
std::size_t find_divisor(const std::vector<Poly>& divisors, const monomial<N>& m) noexcept {
  for (std::size_t i = 0; i < divisors.size(); ++i) {
    if (!divisors[i].is_zero() && monomial<N>::divides(divisors[i].leading_monomial(), m)) {
      return i;
    }
  }
  return divisors.size();
}

}  // namespace detail

// Division of f by an ordered family of divisors (Cox, Little and O'Shea,
// Theorem 2.3.3). The quotients depend on the order of the divisors unless the
// family is a Gröbner basis, in which case the remainder is independent of it
// and is the normal form of f modulo the ideal.
template <class Coeff, std::size_t N, class Order>
division_result<polynomial<Coeff, N, Order>> divide(
    const polynomial<Coeff, N, Order>& f,
    const std::vector<polynomial<Coeff, N, Order>>& divisors) {
  using poly = polynomial<Coeff, N, Order>;
  using mono = monomial<N>;
  using traits = coefficient_traits<Coeff>;

  division_result<poly> result;
  result.quotients.assign(divisors.size(), poly());

  poly p = f;
  std::vector<typename poly::term> remainder_terms;

  while (!p.is_zero()) {
    const typename poly::term head = p.leading_term();
    const std::size_t divisor = detail::find_divisor(divisors, head.mon);

    if (divisor == divisors.size()) {
      // The head is irreducible, so it belongs to the remainder. Successive
      // heads strictly decrease, hence remainder_terms stays sorted.
      remainder_terms.push_back(head);
      p.drop_leading_term();
      continue;
    }

    const poly& g = divisors[divisor];
    const mono shift = mono::divide(head.mon, g.leading_monomial());
    const Coeff factor = head.coeff * traits::inverse(g.leading_coefficient());

    result.quotients[divisor] += poly::from_monomial(shift, factor);
    p.subtract_monomial_multiple(shift, factor, g);
  }

  result.remainder = poly(std::move(remainder_terms));
  return result;
}

// The remainder alone, without accumulating quotients. Modulo a Gröbner basis
// this is the canonical representative of f in the quotient ring, and it
// vanishes exactly when f lies in the ideal. This is the inner loop of
// Buchberger's algorithm, so it avoids the quotient bookkeeping of divide.
template <class Coeff, std::size_t N, class Order>
polynomial<Coeff, N, Order> normal_form(
    const polynomial<Coeff, N, Order>& f,
    const std::vector<polynomial<Coeff, N, Order>>& divisors) {
  using poly = polynomial<Coeff, N, Order>;
  using mono = monomial<N>;
  using traits = coefficient_traits<Coeff>;

  poly p = f;
  std::vector<typename poly::term> remainder_terms;

  while (!p.is_zero()) {
    const typename poly::term head = p.leading_term();
    const std::size_t divisor = detail::find_divisor(divisors, head.mon);

    if (divisor == divisors.size()) {
      remainder_terms.push_back(head);
      p.drop_leading_term();
      continue;
    }

    const poly& g = divisors[divisor];
    p.subtract_monomial_multiple(mono::divide(head.mon, g.leading_monomial()),
                                 head.coeff * traits::inverse(g.leading_coefficient()), g);
  }

  return poly(std::move(remainder_terms));
}

// Reduction of f against every divisor except the one at index skip. Used to
// pass from a minimal Gröbner basis to the reduced one.
template <class Coeff, std::size_t N, class Order>
polynomial<Coeff, N, Order> normal_form_excluding(
    const polynomial<Coeff, N, Order>& f,
    const std::vector<polynomial<Coeff, N, Order>>& divisors, std::size_t skip) {
  std::vector<polynomial<Coeff, N, Order>> others;
  others.reserve(divisors.size());
  for (std::size_t i = 0; i < divisors.size(); ++i) {
    if (i != skip) {
      others.push_back(divisors[i]);
    }
  }
  return normal_form(f, others);
}

}  // namespace varietas

#endif
