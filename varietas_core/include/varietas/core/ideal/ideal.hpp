#ifndef VARIETAS_CORE_IDEAL_IDEAL_HPP
#define VARIETAS_CORE_IDEAL_IDEAL_HPP

#include <cstddef>
#include <utility>
#include <vector>

#include "varietas/core/config.hpp"
#include "varietas/core/ideal/buchberger.hpp"
#include "varietas/core/ideal/division.hpp"
#include "varietas/core/monomial.hpp"
#include "varietas/core/order/order_id.hpp"
#include "varietas/core/polynomial.hpp"
#include "varietas/core/quotient/quotient_basis.hpp"

namespace varietas {

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

  // Generators of the elimination ideal I ∩ k[x_Split, ..., x_{N-1}], namely
  // the basis elements that do not involve any of the first Split variables.
  // By the Elimination Theorem these generate the elimination ideal, provided
  // Order is an elimination order for that split; the assertion below is the
  // only guard the runtime can offer, and varietas_codegen is responsible for
  // choosing the order.
  std::vector<polynomial_type> eliminate(std::size_t split) const {
    VARIETAS_ASSERT(split <= N);
    VARIETAS_ASSERT(Order::id == order_id::lex || Order::id == order_id::block);

    std::vector<polynomial_type> eliminated;
    for (const polynomial_type& g : basis()) {
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

 private:
  std::vector<polynomial_type> generators_;
  mutable std::vector<polynomial_type> basis_;
  mutable buchberger_statistics statistics_;
  mutable bool computed_ = false;
};

}  // namespace varietas

#endif
