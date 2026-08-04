#ifndef VARIETAS_CORE_ORDER_GRLEX_HPP
#define VARIETAS_CORE_ORDER_GRLEX_HPP

#include <cstddef>

#include "varietas/core/monomial.hpp"
#include "varietas/core/order/lex.hpp"
#include "varietas/core/order/order_id.hpp"

namespace varietas {

// Graded lexicographic order: total degree first, lex as tie break.
struct grlex {
  static constexpr order_id id = order_id::grlex;

  template <std::size_t N>
  static constexpr int compare(const monomial<N>& a, const monomial<N>& b) noexcept {
    if (a.degree() != b.degree()) {
      return a.degree() > b.degree() ? 1 : -1;
    }
    return lex::compare(a, b);
  }
};

}  // namespace varietas

#endif
