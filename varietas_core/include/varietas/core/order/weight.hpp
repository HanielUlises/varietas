#ifndef VARIETAS_CORE_ORDER_WEIGHT_HPP
#define VARIETAS_CORE_ORDER_WEIGHT_HPP

#include <cstddef>
#include <cstdint>

#include "varietas/core/monomial.hpp"
#include "varietas/core/order/grevlex.hpp"
#include "varietas/core/order/order_id.hpp"

namespace varietas {

// Weighted order: compare the inner product of the exponent vector with a
// weight vector, and break ties with TieBreak.
//
// Weights is any type exposing
//
//     static constexpr std::int64_t weight(std::size_t variable) noexcept;
//
// The weights must be nonnegative for the relation to be a monomial order;
// with a monomial order as TieBreak, nonnegativity is also sufficient. Weights
// are how the half-angle variables of a kinematic ideal are made cheaper than
// the pose variables without leaving the graded world entirely.
template <class Weights, class TieBreak = grevlex>
struct weight_order {
  static constexpr order_id id = order_id::weight;

  template <std::size_t N>
  static constexpr int compare(const monomial<N>& a, const monomial<N>& b) noexcept {
    std::int64_t wa = 0;
    std::int64_t wb = 0;
    for (std::size_t i = 0; i < N; ++i) {
      const std::int64_t w = Weights::weight(i);
      VARIETAS_ASSERT(w >= 0);
      wa += w * static_cast<std::int64_t>(a[i]);
      wb += w * static_cast<std::int64_t>(b[i]);
    }
    if (wa != wb) {
      return wa > wb ? 1 : -1;
    }
    return TieBreak::compare(a, b);
  }
};

}  // namespace varietas

#endif
