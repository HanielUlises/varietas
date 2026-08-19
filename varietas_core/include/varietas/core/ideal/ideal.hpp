#ifndef VARIETAS_CORE_IDEAL_IDEAL_HPP
#define VARIETAS_CORE_IDEAL_IDEAL_HPP

#include <cstddef>
#include <utility>
#include <vector>

#include "varietas/core/config.hpp"
#include "varietas/core/ideal/buchberger.hpp"
#include "varietas/core/ideal/dimension.hpp"
#include "varietas/core/ideal/division.hpp"
#include "varietas/core/monomial.hpp"
#include "varietas/core/order/order_id.hpp"
#include "varietas/core/polynomial.hpp"
#include "varietas/core/quotient/quotient_basis.hpp"

namespace varietas {

// The elements of a Gröbner basis that do not involve the first Split
// variables. By the Elimination Theorem these generate the elimination ideal
// I ∩ k[x_Split, ..., x_{N-1}], and are a Gröbner basis of it under the order
// induced on the remaining variables, provided Order is an elimination order
// for that split.
//
// Free of the ideal class because a basis is often already in hand — saturation
// returns one, and so does any earlier elimination — and rebuilding an ideal
// around it only to have Buchberger walk the critical pairs of a set that is
// already a Gröbner basis is pure waste.
template <class Coeff, std::size_t N, class Order>
std::vector<polynomial<Coeff, N, Order>> eliminated_generators(
    const std::vector<polynomial<Coeff, N, Order>>& basis, std::size_t split) {
  VARIETAS_ASSERT(split <= N);

  std::vector<polynomial<Coeff, N, Order>> eliminated;
  for (const polynomial<Coeff, N, Order>& g : basis) {
    bool free_of_leading_block = true;
    for (const auto& t : g.terms()) {
      for (std::size_t i = 0; i < split && free_of_leading_block; ++i) {
        free_of_leading_block = t.mon[i] == 0;
      }
      if (!free_of_leading_block) {
        break;
      }
    }
    if (free_of_leading_block) {
      eliminated.push_back(g);
    }
  }
  return eliminated;
}

// An ideal given by generators, with its reduced Gröbner basis computed on
// first demand and cached afterwards. Every question that needs the basis goes
// through basis(), so no caller can accidentally ask a membership question
// against the raw generators.
template <class Coeff, std::size_t N, class Order>
class ideal {
 public:
  using polynomial_type = polynomial<Coeff, N, Order>;
  using monomial_type = monomial<N>;
  using order_type = Order;

  static constexpr std::size_t num_vars = N;

  ideal() = default;

  explicit ideal(std::vector<polynomial_type> generators)
      : generators_(std::move(generators)) {}

  const std::vector<polynomial_type>& generators() const noexcept { return generators_; }

  void add_generator(polynomial_type f) {
    generators_.push_back(std::move(f));
    computed_ = false;
    basis_.clear();
    statistics_ = buchberger_statistics{};
  }

  const std::vector<polynomial_type>& basis() const {
    if (!computed_) {
      basis_ = groebner_basis(generators_, &statistics_);
      computed_ = true;
    }
    return basis_;
  }

  const buchberger_statistics& statistics() const {
    basis();
    return statistics_;
  }

  bool contains(const polynomial_type& f) const { return is_member(f, basis()); }

  // The whole ring, equivalently an empty variety over the algebraic closure.
  bool is_unit() const { return is_unit_ideal(basis()); }

  quotient_basis<N> quotient() const { return standard_monomials(basis()); }

  bool is_zero_dimensional() const { return quotient().is_zero_dimensional; }

  // dim V(I), with the empty variety reported as such. The finiteness verdict
  // above is the same statement at dimension zero — the two are checked against
  // each other in the tests — but it is the one the quotient construction needs
  // and it comes with the standard monomials, so both are kept.
  affine_dimension<N> dimension() const { return ideal_dimension(basis()); }

  // Generators of the elimination ideal I ∩ k[x_Split, ..., x_{N-1}], namely
  // the basis elements that do not involve any of the first Split variables.
  // By the Elimination Theorem these generate the elimination ideal, provided
  // Order is an elimination order for that split; the assertion below is the
  // only guard the runtime can offer, and varietas_codegen is responsible for
  // choosing the order.
  std::vector<polynomial_type> eliminate(std::size_t split) const {
    VARIETAS_ASSERT(Order::id == order_id::lex || Order::id == order_id::block);
    return eliminated_generators(basis(), split);
  }

 private:
  std::vector<polynomial_type> generators_;
  mutable std::vector<polynomial_type> basis_;
  mutable buchberger_statistics statistics_;
  mutable bool computed_ = false;
};

}  // namespace varietas

#endif
