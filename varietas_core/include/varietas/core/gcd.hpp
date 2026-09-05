#ifndef VARIETAS_CORE_GCD_HPP
#define VARIETAS_CORE_GCD_HPP

#include <cstddef>
#include <vector>

#include "varietas/core/config.hpp"
#include "varietas/core/ideal/division.hpp"
#include "varietas/core/monomial.hpp"
#include "varietas/core/polynomial.hpp"

namespace varietas {

// Greatest common divisor in k[x_1, ..., x_N], by the subresultant polynomial
// remainder sequence.
//
// The algorithm is the classical recursion on the number of variables: a
// polynomial in N variables is read as a univariate polynomial in one main
// variable whose coefficients are polynomials in the others, the content, that
// is the gcd of those coefficients, is taken out by a recursive call, and the
// primitive parts are carried through a Euclidean remainder sequence. Ordinary
// remainders do not exist over a ring, so the step is a pseudo-remainder, and
// the growth that pseudo-division causes is removed at each step by the exact
// division that makes the sequence subresultant rather than merely Euclidean.
//
// Everything is expressed in the flat polynomial type rather than in a nested
// univariate one. The coefficient of v^k is a polynomial in the same N
// variables that happens not to involve v, so the recursion needs no second
// representation and no conversion: it is the same type throughout, with one
// fewer variable in play at each level.
//
// The result is monic, which is the canonical representative over a field.

namespace detail {

// The degree of p in variable v alone.
template <class Coeff, std::size_t N, class Order>
int degree_in(const polynomial<Coeff, N, Order>& p, std::size_t v) {
  int d = -1;
  for (const auto& t : p.terms()) {
    d = std::max(d, static_cast<int>(t.mon[v]));
  }
  return d;
}

// The coefficient of v^k: a polynomial in the same variables, with the
// exponent of v struck out.
template <class Coeff, std::size_t N, class Order>
polynomial<Coeff, N, Order> coefficient_in(const polynomial<Coeff, N, Order>& p,
                                           std::size_t v, int k) {
  using poly = polynomial<Coeff, N, Order>;
  std::vector<typename poly::term> terms;
  for (const auto& t : p.terms()) {
    if (static_cast<int>(t.mon[v]) != k) {
      continue;
    }
    auto exponents = t.mon.exponents();
    exponents[v] = 0;
    terms.push_back({monomial<N>(exponents), t.coeff});
  }
  return poly(std::move(terms));
}

// The leading coefficient with respect to v.
template <class Coeff, std::size_t N, class Order>
polynomial<Coeff, N, Order> leading_coefficient_in(const polynomial<Coeff, N, Order>& p,
                                                   std::size_t v) {
  return coefficient_in(p, v, degree_in(p, v));
}

template <class Coeff, std::size_t N, class Order>
polynomial<Coeff, N, Order> monomial_power(std::size_t v, int k) {
  using poly = polynomial<Coeff, N, Order>;
  return poly::from_monomial(monomial<N>::variable(v, static_cast<typename monomial<N>::exponent_type>(k)),
                             coefficient_traits<Coeff>::one());
}

// Exact division, where the divisor is known to divide the dividend. Division
// by a single polynomial leaves no remainder exactly when it divides, so a
// nonzero remainder here is a bug in the caller rather than a case to handle.
template <class Coeff, std::size_t N, class Order>
polynomial<Coeff, N, Order> divide_exact(const polynomial<Coeff, N, Order>& a,
                                         const polynomial<Coeff, N, Order>& b) {
  VARIETAS_ASSERT(!b.is_zero());
  if (a.is_zero()) {
    return polynomial<Coeff, N, Order>();
  }
  const auto result = divide(a, std::vector<polynomial<Coeff, N, Order>>{b});
  VARIETAS_ASSERT(result.remainder.is_zero());
  return result.quotients.front();
}

// Pseudo-remainder in v: the remainder of lc_v(b)^e a divided by b, for the
// least e that makes the division go through over the coefficient ring.
template <class Coeff, std::size_t N, class Order>
polynomial<Coeff, N, Order> pseudo_remainder_in(polynomial<Coeff, N, Order> a,
                                                const polynomial<Coeff, N, Order>& b,
                                                std::size_t v) {
  using poly = polynomial<Coeff, N, Order>;

  const int db = degree_in(b, v);
  VARIETAS_ASSERT(db >= 0);
  const poly lb = leading_coefficient_in(b, v);

  int e = degree_in(a, v) - db + 1;
  if (e < 0) {
    e = 0;
  }
  while (!a.is_zero() && degree_in(a, v) >= db) {
    const int da = degree_in(a, v);
    const poly la = leading_coefficient_in(a, v);
    a = lb * a - la * b * monomial_power<Coeff, N, Order>(v, da - db);
    --e;
  }
  for (int i = 0; i < e; ++i) {
    a = a * lb;
  }
  return a;
}

// The first variable that either polynomial actually involves, or N when
// neither involves any and both are constants.
template <class Coeff, std::size_t N, class Order>
std::size_t main_variable(const polynomial<Coeff, N, Order>& a,
                          const polynomial<Coeff, N, Order>& b) {
  for (std::size_t v = 0; v < N; ++v) {
    if (degree_in(a, v) > 0 || degree_in(b, v) > 0) {
      return v;
    }
  }
  return N;
}

}  // namespace detail

template <class Coeff, std::size_t N, class Order>
polynomial<Coeff, N, Order> polynomial_gcd(const polynomial<Coeff, N, Order>& a,
                                           const polynomial<Coeff, N, Order>& b);

namespace detail {

// The gcd of the coefficients of p with respect to v, which is the content of p
// read as a univariate polynomial in v. The recursion is here: every
// coefficient is free of v, so the call below runs with one fewer variable in
// play.
template <class Coeff, std::size_t N, class Order>
polynomial<Coeff, N, Order> content_in(const polynomial<Coeff, N, Order>& p, std::size_t v) {
  using poly = polynomial<Coeff, N, Order>;
  poly result;
  for (int k = degree_in(p, v); k >= 0; --k) {
    const poly c = coefficient_in(p, v, k);
    if (c.is_zero()) {
      continue;
    }
    result = result.is_zero() ? c : polynomial_gcd(result, c);
    if (result.degree() == 0) {
      break;  // a unit: nothing further can divide it
    }
  }
  return result;
}

}  // namespace detail

// The greatest common divisor, monic.
//
// Only meaningful over an exact field: over floating point a coefficient that
// should cancel to zero leaves a residue instead, the remainder sequence never
// terminates where it should, and the answer is a statement about a different
// pair of polynomials. The refusal is at compile time, as it is for
// polynomial::prune in the opposite direction.
template <class Coeff, std::size_t N, class Order>
polynomial<Coeff, N, Order> polynomial_gcd(const polynomial<Coeff, N, Order>& a,
                                           const polynomial<Coeff, N, Order>& b) {
  static_assert(coefficient_traits<Coeff>::is_exact,
                "polynomial_gcd is only sound over an exact coefficient field: over floating "
                "point the remainder sequence does not terminate where it should");

  using poly = polynomial<Coeff, N, Order>;

  if (a.is_zero()) {
    return b.is_zero() ? poly() : b.monic();
  }
  if (b.is_zero()) {
    return a.monic();
  }

  const std::size_t v = detail::main_variable(a, b);
  if (v == N) {
    // Both are nonzero constants, and every nonzero constant is a unit.
    return poly::constant(coefficient_traits<Coeff>::one());
  }

  poly A = a;
  poly B = b;
  if (detail::degree_in(A, v) < detail::degree_in(B, v)) {
    std::swap(A, B);
  }

  // Take the contents out, so the remainder sequence runs on primitive parts
  // and the content of the answer is restored at the end.
  const poly content_a = detail::content_in(A, v);
  const poly content_b = detail::content_in(B, v);
  const poly common_content = polynomial_gcd(content_a, content_b);

  A = detail::divide_exact(A, content_a);
  B = detail::divide_exact(B, content_b);

  // The subresultant remainder sequence. The division by g h^delta is what
  // keeps the coefficients from growing the way a bare pseudo-remainder
  // sequence makes them grow, and it is exact at every step.
  poly g = poly::constant(coefficient_traits<Coeff>::one());
  poly h = poly::constant(coefficient_traits<Coeff>::one());

  while (true) {
    const int delta = detail::degree_in(A, v) - detail::degree_in(B, v);
    const poly remainder = detail::pseudo_remainder_in(A, B, v);

    if (remainder.is_zero()) {
      break;
    }
    if (detail::degree_in(remainder, v) == 0) {
      B = poly::constant(coefficient_traits<Coeff>::one());
      break;
    }

    A = B;

    poly divisor = g;
    for (int i = 0; i < delta; ++i) {
      divisor = divisor * h;
    }
    B = detail::divide_exact(remainder, divisor);

    g = detail::leading_coefficient_in(A, v);

    // h <- g^delta h^(1-delta), computed as an exact division so that it stays
    // inside the ring for every delta.
    poly numerator = poly::constant(coefficient_traits<Coeff>::one());
    for (int i = 0; i < delta; ++i) {
      numerator = numerator * g;
    }
    if (delta == 0) {
      // h is unchanged
    } else if (delta == 1) {
      h = g;
    } else {
      poly denominator = poly::constant(coefficient_traits<Coeff>::one());
      for (int i = 0; i < delta - 1; ++i) {
        denominator = denominator * h;
      }
      h = detail::divide_exact(numerator, denominator);
    }
  }

  const poly primitive = detail::divide_exact(B, detail::content_in(B, v));
  return (common_content * primitive).monic();
}

}  // namespace varietas

#endif
