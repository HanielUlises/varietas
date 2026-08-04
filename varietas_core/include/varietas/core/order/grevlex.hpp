#ifndef VARIETAS_CORE_ORDER_GREVLEX_HPP
#define VARIETAS_CORE_ORDER_GREVLEX_HPP

#include <cstddef>

#include "varietas/core/monomial.hpp"
#include "varietas/core/order/order_id.hpp"

namespace varietas {

// Graded reverse lexicographic order: total degree first, then the monomial
// with the smaller exponent in the last differing variable is the larger one.
struct grevlex {
  static constexpr order_id id = order_id::grevlex;

  template <std::size_t N>
  static constexpr int compare(const monomial<N>& a, const monomial<N>& b) noexcept {
    if (a.degree() != b.degree()) {
      return a.degree() > b.degree() ? 1 : -1;
    }
    for (std::size_t i = N; i > 0; --i) {
      const std::size_t k = i - 1;
      if (a[k] != b[k]) {
        return a[k] < b[k] ? 1 : -1;
      }
    }
    return 0;
  }
};

}  // namespace varietas

#endif
