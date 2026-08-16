#ifndef VARIETAS_CORE_IDEAL_DIVISION_HPP
#define VARIETAS_CORE_IDEAL_DIVISION_HPP

#include <cstddef>
#include <vector>

#include "varietas/core/config.hpp"
#include "varietas/core/monomial.hpp"
#include "varietas/core/polynomial.hpp"

namespace varietas {

// Output of the multivariate division algorithm: f = sum_i q_i g_i + r, with
// no monomial of r divisible by any leading monomial of the divisors.
template <class Poly>
struct division_result {
  std::vector<Poly> quotients;
  Poly remainder;
};

// Division of f by an ordered family of divisors (Cox, Little, O'Shea,
// Theorem 2.3.3). The result depends on the order of the divisors unless the
// family is a Gröbner basis, in which case the remainder is the unique normal
// form of f modulo the ideal.
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

    std::size_t divisor = divisors.size();
    for (std::size_t i = 0; i < divisors.size(); ++i) {
      if (!divisors[i].is_zero() && mono::divides(divisors[i].leading_monomial(), head.mon)) {
        divisor = i;
        break;
      }
    }

    if (divisor == divisors.size()) {
      // The leading term is irreducible; it belongs to the remainder and is
      // removed so that division proceeds on the tail. Successive heads are
      // strictly decreasing, so remainder_terms stays sorted.
      remainder_terms.push_back(head);
      p.subtract_monomial_multiple(mono::one(), traits::one(), poly::from_monomial(head.mon, head.coeff));
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

// The remainder alone. Modulo a Gröbner basis this is the canonical
// representative of f in the quotient ring, and it vanishes exactly when f
// lies in the ideal.
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

    std::size_t divisor = divisors.size();
    for (std::size_t i = 0; i < divisors.size(); ++i) {
      if (!divisors[i].is_zero() && mono::divides(divisors[i].leading_monomial(), head.mon)) {
        divisor = i;
        break;
      }
    }

    if (divisor == divisors.size()) {
      remainder_terms.push_back(head);
      p.subtract_monomial_multiple(mono::one(), traits::one(), poly::from_monomial(head.mon, head.coeff));
      continue;
    }

    const poly& g = divisors[divisor];
    p.subtract_monomial_multiple(mono::divide(head.mon, g.leading_monomial()),
                                 head.coeff * traits::inverse(g.leading_coefficient()), g);
  }

  return poly(std::move(remainder_terms));
}

}  // namespace varietas

#endif
