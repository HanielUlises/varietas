#ifndef VARIETAS_CORE_CONFIG_HPP
#define VARIETAS_CORE_CONFIG_HPP

#include <cassert>
#include <cstddef>
#include <cstdint>

#define VARIETAS_ASSERT(condition) assert(condition)

namespace varietas {

using index_type = std::size_t;
using scalar = double;

// Specialized by varietas_codegen for its exact rational type. varietas_core
// instantiates it only with double.
template <class Coeff>
struct coefficient_traits;

template <>
struct coefficient_traits<double> {
  // False for double, true for the rational type used offline. Algorithms that
  // are only sound over an exact field consult this flag rather than assuming.
  static constexpr bool is_exact = false;

  static constexpr double zero() noexcept { return 0.0; }
  static constexpr double one() noexcept { return 1.0; }
  static constexpr bool is_zero(double c) noexcept { return c == 0.0; }
  static constexpr double inverse(double c) noexcept { return 1.0 / c; }
  static constexpr double negate(double c) noexcept { return -c; }
  static constexpr double to_double(double c) noexcept { return c; }
};

}  // namespace varietas

#endif
