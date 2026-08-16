#ifndef VARIETAS_CORE_QUOTIENT_QUOTIENT_BASIS_HPP
#define VARIETAS_CORE_QUOTIENT_QUOTIENT_BASIS_HPP

#include <algorithm>
#include <array>
#include <cstddef>
#include <vector>

#include "varietas/core/config.hpp"
#include "varietas/core/monomial.hpp"
#include "varietas/core/polynomial.hpp"

namespace varietas {

// Monomial basis of the quotient ring k[x]/I, together with the verdict on
// whether that quotient is finite dimensional at all.
//
// The standard monomials are those lying outside the ideal of leading terms;
// by Macaulay's theorem their classes form a basis of the quotient as a vector
// space. The quotient is finite dimensional exactly when, for every variable,
// some element of the basis has a pure power of that variable as its leading
// monomial (the finiteness criterion of Cox, Little and O'Shea, Chapter 5,
// Theorem 5.3.6). For a kinematic ideal this is the algebraic statement that
// the pose admits finitely many joint solutions.
template <std::size_t N>
struct quotient_basis {
  bool is_zero_dimensional = false;
  // Ascending under the order, so index 0 is the class of 1 whenever the
  // quotient is nontrivial.
  std::vector<monomial<N>> monomials;

  std::size_t dimension() const noexcept { return monomials.size(); }

  // Position of m in the basis, or dimension() if m is not standard.
  std::size_t index_of(const monomial<N>& m) const noexcept {
    for (std::size_t i = 0; i < monomials.size(); ++i) {
      if (monomials[i] == m) {
        return i;
      }
    }
    return monomials.size();
  }

  bool contains(const monomial<N>& m) const noexcept { return index_of(m) < monomials.size(); }
};

// Computes the standard monomials of the ideal whose Gröbner basis is given.
// If the quotient is not finite dimensional the verdict is recorded and the
// monomial list is left empty, since no finite basis exists to report.
template <class Coeff, std::size_t N, class Order>
quotient_basis<N> standard_monomials(const std::vector<polynomial<Coeff, N, Order>>& basis) {
  using mono = monomial<N>;
  using exponent_type = typename mono::exponent_type;

  quotient_basis<N> result;

  // The unit ideal: the quotient is the zero ring, finite dimensional of
  // dimension zero, with no standard monomials.
  for (const auto& g : basis) {
    if (!g.is_zero() && g.leading_monomial().is_one()) {
      result.is_zero_dimensional = true;
      return result;
    }
  }

  // Smallest pure power x_i^d appearing among the leading monomials; absent
  // for some variable, the quotient is infinite dimensional.
  std::array<std::size_t, N> bound{};
  for (std::size_t i = 0; i < N; ++i) {
    bound[i] = 0;
    for (const auto& g : basis) {
      if (g.is_zero()) {
        continue;
      }
      const mono& lm = g.leading_monomial();
      if (lm.degree() == lm[i] && lm[i] > 0) {
        if (bound[i] == 0 || lm[i] < bound[i]) {
          bound[i] = lm[i];
        }
      }
    }
    if (bound[i] == 0) {
      result.is_zero_dimensional = false;
      return result;
    }
  }

  result.is_zero_dimensional = true;

  // Every standard monomial has exponent below the pure-power bound in each
  // variable, so the search is confined to that box.
  std::array<exponent_type, N> exponents{};
  for (;;) {
    const mono candidate(exponents);
    const bool reducible = std::any_of(basis.begin(), basis.end(), [&](const auto& g) {
      return !g.is_zero() && mono::divides(g.leading_monomial(), candidate);
    });
    if (!reducible) {
      result.monomials.push_back(candidate);
    }

    std::size_t carry = 0;
    while (carry < N) {
      exponents[carry] = static_cast<exponent_type>(exponents[carry] + 1);
      if (exponents[carry] < bound[carry]) {
        break;
      }
      exponents[carry] = 0;
      ++carry;
    }
    if (carry == N) {
      break;
    }
  }

  std::sort(result.monomials.begin(), result.monomials.end(),
            [](const mono& a, const mono& b) { return Order::compare(a, b) < 0; });

  return result;
}

}  // namespace varietas

#endif
