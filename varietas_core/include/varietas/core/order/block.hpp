#ifndef VARIETAS_CORE_ORDER_BLOCK_HPP
#define VARIETAS_CORE_ORDER_BLOCK_HPP

#include <array>
#include <cstddef>

#include "varietas/core/monomial.hpp"
#include "varietas/core/order/order_id.hpp"

namespace varietas {

// Block (product) order on N variables split after the first Split of them:
// compare the leading block under First, and only on a tie compare the
// trailing block under Second.
//
// This is the elimination order of the theory. Because First is a monomial
// order on the leading block, the monomial 1 is minimal there, hence every
// monomial that involves a variable of the leading block dominates every
// monomial that does not. Consequently, for a Gröbner basis G of an ideal I
// under this order, the elements of G free of the first Split variables
// generate the elimination ideal I intersected with k[x_Split, ..., x_{N-1}],
// which is the mechanism used to project a kinematic ideal onto the joint
// variables.
template <std::size_t Split, class First, class Second>
struct block_order {
  static constexpr order_id id = order_id::block;
  static constexpr std::size_t split = Split;

  template <std::size_t N>
  static constexpr int compare(const monomial<N>& a, const monomial<N>& b) noexcept {
    static_assert(Split <= N, "block split exceeds the number of variables");

    const int leading = First::compare(project<N, 0, Split>(a), project<N, 0, Split>(b));
    if (leading != 0) {
      return leading;
    }
    return Second::compare(project<N, Split, N - Split>(a), project<N, Split, N - Split>(b));
  }

 private:
  template <std::size_t N, std::size_t Offset, std::size_t Count>
  static constexpr monomial<Count> project(const monomial<N>& m) noexcept {
    std::array<typename monomial<N>::exponent_type, Count> exponents{};
    for (std::size_t i = 0; i < Count; ++i) {
      exponents[i] = m[Offset + i];
    }
    return monomial<Count>(exponents);
  }
};

}  // namespace varietas

#endif
