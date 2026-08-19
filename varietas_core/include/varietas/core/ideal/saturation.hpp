#ifndef VARIETAS_CORE_IDEAL_SATURATION_HPP
#define VARIETAS_CORE_IDEAL_SATURATION_HPP

#include <cstddef>
#include <utility>
#include <vector>

#include "varietas/core/config.hpp"
#include "varietas/core/embed.hpp"
#include "varietas/core/ideal/buchberger.hpp"
#include "varietas/core/monomial.hpp"
#include "varietas/core/order/block.hpp"
#include "varietas/core/order/lex.hpp"
#include "varietas/core/polynomial.hpp"

namespace varietas {

// The saturation I : h^∞ = { f : f h^k ∈ I for some k }, computed by
// Rabinowitsch's trick.
//
// Geometrically, saturating removes from V(I) every component contained in the
// hypersurface h = 0: over an algebraically closed field V(I : h^∞) is the
// Zariski closure of V(I) \ V(h). That is exactly the operation the half-angle
// substitution requires. Writing t = tan(q/2) and clearing the denominators
// 1 + t_i^2 introduces components on which some 1 + t_i^2 vanishes, that is,
// t_i = ±i — the images of the configurations q_i = π that the substitution
// sent to infinity. They are not configurations of the robot, and they are
// invisible over the reals, where 1 + t^2 never vanishes. Left in the ideal
// they inflate dim_k A, and the finiteness verdict then counts solutions the
// arm does not have, which would make the completeness certificate a statement
// about the wrong variety.
//
// The computation adjoins one variable y and the generator 1 - y·h. In the
// larger ring, V(I + (1 - y h)) is the part of V(I) where h ≠ 0, with y
// recording 1/h; eliminating y projects that back down, and the Closure Theorem
// identifies the result as the saturation. So the whole construction is an
// elimination, and it uses the block order machinery exactly as it stands.
//
// The order matters twice over. Adjoined variable first, so that the block
// order eliminates it; and the trailing block ordered by the caller's own
// Order, so that the surviving elements are a Gröbner basis of the saturation
// under that same order, ready to use without a second Buchberger run.

namespace detail {

// y is variable 0 of the enlarged ring; the original variables follow it.
template <std::size_t N, class Order>
using saturation_order = block_order<1, lex, Order>;

constexpr std::size_t saturation_offset = 1;

}  // namespace detail

// Generators of I : h^∞, as a reduced Gröbner basis under Order.
//
// Requires h to be nonzero. Saturating by zero is the whole ring by the
// definition above, which is a degenerate answer that no caller wants and that
// almost always indicates the multiplier was computed rather than intended.
template <class Coeff, std::size_t N, class Order>
std::vector<polynomial<Coeff, N, Order>> saturate(
    const std::vector<polynomial<Coeff, N, Order>>& generators,
    const polynomial<Coeff, N, Order>& h, buchberger_statistics* statistics = nullptr) {
  using poly = polynomial<Coeff, N, Order>;
  using enlarged_order = detail::saturation_order<N, Order>;
  using enlarged = polynomial<Coeff, N + 1, enlarged_order>;

  VARIETAS_ASSERT(!h.is_zero());

  // The zero ideal saturates to itself: f·h^k = 0 with h ≠ 0 forces f = 0 in an
  // integral domain, and there is nothing for Buchberger to do.
  bool all_zero = true;
  for (const poly& f : generators) {
    all_zero = all_zero && f.is_zero();
  }
  if (all_zero) {
    return {};
  }

  std::vector<enlarged> lifted;
  lifted.reserve(generators.size() + 1);
  for (const poly& f : generators) {
    if (!f.is_zero()) {
      lifted.push_back(
          embed_polynomial<enlarged_order, N + 1>(f, detail::saturation_offset));
    }
  }

  // 1 - y·h, the generator that inverts h.
  const enlarged y = enlarged::variable(0);
  const enlarged lifted_h =
      embed_polynomial<enlarged_order, N + 1>(h, detail::saturation_offset);
  lifted.push_back(enlarged::constant(coefficient_traits<Coeff>::one()) - y * lifted_h);

  const std::vector<enlarged> basis = groebner_basis(lifted, statistics);

  // Elimination: the elements free of y generate I : h^∞, and because the
  // trailing block is ordered by Order they are already a Gröbner basis of it
  // under Order.
  std::vector<poly> saturated;
  for (const enlarged& g : basis) {
    bool involves_y = false;
    for (const auto& t : g.terms()) {
      if (t.mon[0] != 0) {
        involves_y = true;
        break;
      }
    }
    if (!involves_y) {
      saturated.push_back(
          project_polynomial<Order, N>(g, detail::saturation_offset));
    }
  }

  // Already a Gröbner basis; this only puts it in the canonical reduced form,
  // so that two saturations of the same ideal compare equal.
  reduce_groebner_basis(saturated);
  return saturated;
}

// Saturation by a product, which is the form the kinematic ideal calls for:
// the denominators 1 + t_i^2 cleared by the half-angle substitution, one for
// each revolute joint. Saturating by the product removes every component on
// which any one of them vanishes, since V(h1···hm) is the union of the V(hi).
template <class Coeff, std::size_t N, class Order>
std::vector<polynomial<Coeff, N, Order>> saturate_by_product(
    const std::vector<polynomial<Coeff, N, Order>>& generators,
    const std::vector<polynomial<Coeff, N, Order>>& multipliers,
    buchberger_statistics* statistics = nullptr) {
  using poly = polynomial<Coeff, N, Order>;

  poly product = poly::constant(coefficient_traits<Coeff>::one());
  for (const poly& h : multipliers) {
    VARIETAS_ASSERT(!h.is_zero());
    product = product * h;
  }
  return saturate(generators, product, statistics);
}

// Splitting a variety along a hypersurface.
//
// V(I) = V(I : h^∞) ∪ V(I + (h)), exactly, for any h: a point of V(I) either
// lies on h = 0, and then it is in the second piece, or it does not, and then
// it lies in the closure of V(I) \ V(h), which is the first. The identity is
// the one honest step towards a decomposition that this library can take with
// what it has, and it is worth being precise about what it is not.
//
// It is not a primary decomposition. The pieces are not irreducible, they are
// not primary, and nothing here chooses h — a real decomposition algorithm
// finds its own splittings, and finding them needs factorisation of
// multivariate polynomials over Q, which varietas does not have and which is a
// project of its own. What the caller gets is the ability to say "separate the
// configurations where this vanishes from the ones where it does not" and have
// the two answers be exact and exhaustive.
//
// It is nonetheless what a decomposition is made of. Applied to the singular
// locus of an arm with a generator the geometry suggests — the sine of an elbow
// angle, say — it separates the elbow-straight branch from the elbow-folded one
// and each piece can then be measured, eliminated to the workspace, or split
// again.
template <class Coeff, std::size_t N, class Order>
struct ideal_splitting {
  // I : h^∞, the closure of the part of V(I) away from h = 0.
  std::vector<polynomial<Coeff, N, Order>> away;
  // I + (h), the part of V(I) lying on h = 0.
  std::vector<polynomial<Coeff, N, Order>> on;
};

template <class Coeff, std::size_t N, class Order>
ideal_splitting<Coeff, N, Order> split_along(
    const std::vector<polynomial<Coeff, N, Order>>& generators,
    const polynomial<Coeff, N, Order>& h, buchberger_statistics* statistics = nullptr) {
  VARIETAS_ASSERT(!h.is_zero());

  ideal_splitting<Coeff, N, Order> result;
  result.away = saturate(generators, h, statistics);

  std::vector<polynomial<Coeff, N, Order>> on_generators = generators;
  on_generators.push_back(h);
  result.on = groebner_basis(on_generators);
  return result;
}

// 1 + t_i^2, the denominator the tangent half-angle substitution clears at the
// i-th variable. Named here because the saturation it calls for is the reason
// this header exists.
template <class Coeff, std::size_t N, class Order>
polynomial<Coeff, N, Order> half_angle_denominator(std::size_t i) {
  using poly = polynomial<Coeff, N, Order>;
  VARIETAS_ASSERT(i < N);
  return poly::constant(coefficient_traits<Coeff>::one()) +
         poly::variable(i, 2);
}

}  // namespace varietas

#endif
