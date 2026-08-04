#ifndef VARIETAS_CORE_ORDER_ORDER_ID_HPP
#define VARIETAS_CORE_ORDER_ORDER_ID_HPP

#include <cstdint>

namespace varietas {

// Recorded in every generated header so that a mismatch between the order used
// offline and the order assumed at runtime is a compile error rather than a
// wrong solution set.
enum class order_id : std::uint8_t {
  lex = 0,
  grlex = 1,
  grevlex = 2,
  block = 3,
  weight = 4,
};

constexpr const char* to_string(order_id id) noexcept {
  return id == order_id::lex       ? "lex"
         : id == order_id::grlex   ? "grlex"
         : id == order_id::grevlex ? "grevlex"
         : id == order_id::block   ? "block"
         : id == order_id::weight  ? "weight"
                                   : "unknown";
}

}  // namespace varietas

#endif
