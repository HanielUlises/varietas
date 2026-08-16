#ifndef VARIETAS_CORE_QUOTIENT_ACTION_MATRIX_HPP
#define VARIETAS_CORE_QUOTIENT_ACTION_MATRIX_HPP

#include <cstddef>
#include <vector>

#include <Eigen/Core>

#include "varietas/core/config.hpp"
#include "varietas/core/ideal/division.hpp"
#include "varietas/core/monomial.hpp"
#include "varietas/core/polynomial.hpp"
#include "varietas/core/quotient/quotient_basis.hpp"

namespace varietas {

using action_matrix_type = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;

// Matrix of the multiplication operator m_f : A -> A, a |-> f a, on the
// quotient algebra A = k[x]/I written in the basis of standard monomials.
//
// Column j holds the coordinates of the normal form of f b_j, so that for a
// coordinate vector v one has coords(f a) = M coords(a). This is the object
// the generated header stores: at runtime the inverse kinematics is the
// spectral decomposition of a matrix of this fixed size, with no iteration and
// no initial guess.
//
// The basis must come from a reduced Gröbner basis of I under the same order,
// otherwise the normal forms are not canonical and the matrix is meaningless.
template <class Coeff, std::size_t N, class Order>
action_matrix_type action_matrix(const polynomial<Coeff, N, Order>& multiplier,
                                 const std::vector<polynomial<Coeff, N, Order>>& basis,
                                 const quotient_basis<N>& quotient) {
  using poly = polynomial<Coeff, N, Order>;
  using traits = coefficient_traits<Coeff>;

  VARIETAS_ASSERT(quotient.is_zero_dimensional);

  const std::size_t d = quotient.dimension();
  action_matrix_type matrix = action_matrix_type::Zero(static_cast<Eigen::Index>(d),
                                                       static_cast<Eigen::Index>(d));

  for (std::size_t j = 0; j < d; ++j) {
    const poly image = normal_form(multiplier * poly::from_monomial(quotient.monomials[j],
                                                                    traits::one()),
                                   basis);
    for (const auto& t : image.terms()) {
      const std::size_t i = quotient.index_of(t.mon);
      // A normal form is supported on standard monomials, so the lookup can
      // only fail if the basis is not a Gröbner basis of the ideal.
      VARIETAS_ASSERT(i < d);
      matrix(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(j)) =
          traits::to_double(t.coeff);
    }
  }

  return matrix;
}

// Convenience overload for the multiplication operator of a single variable,
// which is the case the solver actually uses.
template <class Coeff, std::size_t N, class Order>
action_matrix_type variable_action_matrix(std::size_t variable,
                                          const std::vector<polynomial<Coeff, N, Order>>& basis,
                                          const quotient_basis<N>& quotient) {
  VARIETAS_ASSERT(variable < N);
  return action_matrix(polynomial<Coeff, N, Order>::variable(variable), basis, quotient);
}

}  // namespace varietas

#endif
