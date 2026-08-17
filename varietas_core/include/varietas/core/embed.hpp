#ifndef VARIETAS_CORE_EMBED_HPP
#define VARIETAS_CORE_EMBED_HPP

#include <cstddef>
#include <vector>

#include "varietas/core/config.hpp"
#include "varietas/core/monomial.hpp"
#include "varietas/core/polynomial.hpp"

namespace varietas {

// Moving polynomials between rings of different numbers of variables.
//
// Several constructions work by adjoining variables to the ring, computing
// there, and reading the answer back in the original one: saturation adjoins
// the Rabinowitsch variable, and the parametric formulation will adjoin the
// pose coordinates. Both need the same two operations, and both need them to
// be honest about the order — a polynomial carries its monomial order in its
// type precisely so that terms written under different orders cannot be mixed,
// and an embedding lands in a ring whose order is a different one, so the terms
// must be re-sorted rather than copied across.

// The monomial of M variables whose exponents are those of m placed at offset,
// the remaining variables not appearing.
template <std::size_t M, std::size_t N>
monomial<M> embed_at(const monomial<N>& m, std::size_t offset) {
  static_assert(M >= N, "embedding must not lose variables");
  VARIETAS_ASSERT(offset + N <= M);

  std::array<typename monomial<M>::exponent_type, M> exponents{};
  for (std::size_t i = 0; i < N; ++i) {
    exponents[offset + i] = m[i];
  }
  return monomial<M>(exponents);
}

// The inverse, defined only on monomials that do not involve the variables
// outside the window. Callers establish that by construction — after an
// elimination, the surviving basis elements are free of the eliminated block —
// so a violation is a programming error rather than a runtime condition.
template <std::size_t N, std::size_t M>
monomial<N> project_from(const monomial<M>& m, std::size_t offset) {
  static_assert(N <= M, "projection must not invent variables");
  VARIETAS_ASSERT(offset + N <= M);

  std::array<typename monomial<N>::exponent_type, N> exponents{};
  for (std::size_t i = 0; i < M; ++i) {
    if (i < offset || i >= offset + N) {
      VARIETAS_ASSERT(m[i] == 0);
    } else {
      exponents[i - offset] = m[i];
    }
  }
  return monomial<N>(exponents);
}

// The same polynomial read in a larger ring, under that ring's order. The
// target order is explicit rather than deduced: the whole purpose of enlarging
// the ring is usually to work under an elimination order, which is not the
// order the input was written under.
template <class TargetOrder, std::size_t M, class Coeff, std::size_t N, class Order>
polynomial<Coeff, M, TargetOrder> embed_polynomial(const polynomial<Coeff, N, Order>& f,
                                                   std::size_t offset) {
  using target = polynomial<Coeff, M, TargetOrder>;

  std::vector<typename target::term> terms;
  terms.reserve(f.size());
  for (const auto& t : f.terms()) {
    terms.push_back(typename target::term{embed_at<M>(t.mon, offset), t.coeff});
  }
  // The constructor sorts under TargetOrder, which is the point: the term list
  // of the input is sorted under a different order and is not reusable.
  return target(std::move(terms));
}

template <class TargetOrder, std::size_t N, class Coeff, std::size_t M, class Order>
polynomial<Coeff, N, TargetOrder> project_polynomial(const polynomial<Coeff, M, Order>& f,
                                                     std::size_t offset) {
  using target = polynomial<Coeff, N, TargetOrder>;

  std::vector<typename target::term> terms;
  terms.reserve(f.size());
  for (const auto& t : f.terms()) {
    terms.push_back(typename target::term{project_from<N>(t.mon, offset), t.coeff});
  }
  return target(std::move(terms));
}

}  // namespace varietas

#endif
