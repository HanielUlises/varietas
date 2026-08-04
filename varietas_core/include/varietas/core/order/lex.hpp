#ifndef VARIETAS_CORE_ORDER_LEX_HPP
#define VARIETAS_CORE_ORDER_LEX_HPP

#include <cstddef>

#include "varietas/core/monomial.hpp"
#include "varietas/core/order/order_id.hpp"

namespace varietas {

// Lexicographic order. compare returns -1 if a < b, 0 if a == b, 1 if a > b.
struct lex {
  static constexpr order_id id = order_id::lex;

  template <std::size_t N>
  static constexpr int compare(const monomial<N>& a, const monomial<N>& b) noexcept {
    for (std::size_t i = 0; i < N; ++i) {
      if (a[i] != b[i]) {
        return a[i] > b[i] ? 1 : -1;
      }
    }
    return 0;
  }
};

}  // namespace varietas

#endif
