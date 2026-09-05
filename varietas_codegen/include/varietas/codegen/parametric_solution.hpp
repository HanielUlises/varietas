#ifndef VARIETAS_CODEGEN_PARAMETRIC_SOLUTION_HPP
#define VARIETAS_CODEGEN_PARAMETRIC_SOLUTION_HPP

#include <cstddef>
#include <string>
#include <vector>

#include "varietas/codegen/rational_function.hpp"
#include "varietas/core/config.hpp"
#include "varietas/core/ideal/division.hpp"
#include "varietas/core/order/order_id.hpp"
#include "varietas/core/polynomial.hpp"
#include "varietas/core/monomial.hpp"
#include "varietas/core/quotient/quotient_basis.hpp"

namespace varietas {
namespace codegen {

// The action matrix with the pose left symbolic.
//
// varietas_core already builds action matrices, but its version ends each entry
// with traits::to_double, which is the right thing when the basis was computed
// over Q and the eigendecomposition is the next step. Over Q(x, ..) it is not
// available at all: a rational function has no value until a pose is supplied,
// and coefficient_traits<rational_function>::to_double says so with an assert.
//
// So the construction is repeated here with the conversion left out. It is the
// same construction, where the j-th column is the normal form of x_v times the
// j-th standard monomial read off in the quotient basis, and the only
// difference is what the entry is stored as. That is the whole of the offline half: the
// matrix whose entries are functions of the pose is the object the emitter
// writes out, and evaluating it at a pose is what the generated header does.
template <std::size_t P>
struct parametric_matrix {
  std::size_t dimension = 0;

  // Column-major, dimension * dimension entries, to match the layout the
  // emitted header writes and Eigen's default.
  std::vector<rational_function<P>> entries;

  const rational_function<P>& operator()(std::size_t row, std::size_t column) const {
    VARIETAS_ASSERT(row < dimension && column < dimension);
    return entries[column * dimension + row];
  }

  rational_function<P>& operator()(std::size_t row, std::size_t column) {
    VARIETAS_ASSERT(row < dimension && column < dimension);
    return entries[column * dimension + row];
  }
};

template <std::size_t N, std::size_t P, class Order>
parametric_matrix<P> parametric_action_matrix(
    std::size_t variable,
    const std::vector<polynomial<rational_function<P>, N, Order>>& basis,
    const quotient_basis<N>& quotient) {
  using field = rational_function<P>;
  using poly = polynomial<field, N, Order>;
  using traits = coefficient_traits<field>;

  VARIETAS_ASSERT(variable < N);
  VARIETAS_ASSERT(quotient.is_zero_dimensional);

  const std::size_t d = quotient.dimension();
  parametric_matrix<P> matrix;
  matrix.dimension = d;
  matrix.entries.assign(d * d, traits::zero());

  const poly multiplier = poly::variable(variable);
  for (std::size_t j = 0; j < d; ++j) {
    const poly image =
        normal_form(multiplier * poly::from_monomial(quotient.monomials[j], traits::one()), basis);
    for (const auto& t : image.terms()) {
      const std::size_t i = quotient.index_of(t.mon);
      // The normal form of an element of the quotient is supported on the
      // standard monomials by construction; anything else means the basis
      // handed in was not a Gröbner basis of the ideal the quotient came from.
      VARIETAS_ASSERT(i < d);
      matrix(i, j) = t.coeff;
    }
  }
  return matrix;
}

// The coordinates of the normal form of each variable in the standard basis,
// with the pose left symbolic.
//
// A variable is usually not itself a standard monomial: an ideal that forces
// one joint to be a function of the others reduces it away, which is exactly
// what happened to u in the fixture these tests use. So recovering x_i at a
// point is not a matter of reading off a coordinate, it is evaluating the
// normal form of x_i there, and over Q(p) those coefficients are functions of
// the pose, which is why they have to be emitted alongside the matrices rather
// than baked in as numbers.
template <std::size_t N, std::size_t P, class Order>
std::vector<std::vector<rational_function<P>>> parametric_variable_coordinates(
    const std::vector<polynomial<rational_function<P>, N, Order>>& basis,
    const quotient_basis<N>& quotient) {
  using field = rational_function<P>;
  using poly = polynomial<field, N, Order>;
  using traits = coefficient_traits<field>;

  const std::size_t d = quotient.dimension();
  std::vector<std::vector<field>> coordinates(N, std::vector<field>(d, traits::zero()));
  for (std::size_t i = 0; i < N; ++i) {
    const poly reduced = normal_form(poly::variable(i), basis);
    for (const auto& t : reduced.terms()) {
      const std::size_t k = quotient.index_of(t.mon);
      VARIETAS_ASSERT(k < d);
      coordinates[i][k] = t.coeff;
    }
  }
  return coordinates;
}

// Everything the emitter needs about a solved system, and nothing about how it
// was posed. Keeping it to this means the emitter is testable without a robot:
// a parametric_solution can be built by hand in a unit test.
template <std::size_t N, std::size_t P>
struct parametric_solution {
  // The order the basis was computed under, stored so that the generated header
  // can be checked against the runtime's assumption rather than trusted.
  order_id order = order_id::grevlex;

  // Names used for the unknowns in comments in the generated header. Empty is
  // allowed; the emitter falls back to x0, x1, ...
  std::vector<std::string> unknown_names;

  // Likewise for the pose parameters, which appear in the generated signature.
  std::vector<std::string> parameter_names;

  quotient_basis<N> quotient;

  // One per unknown, indexed as the unknowns are.
  std::vector<parametric_matrix<P>> action;

  // Index of the monomial 1 in the quotient basis. 1 is a standard monomial for
  // every proper ideal, since no leading monomial of a nonunit basis is
  // constant, and the eigenvalue method divides by the eigenvector's component
  // there. It is the value the evaluation functional gives to 1, and a
  // functional that sends 1 to zero is not an evaluation at a point.
  std::size_t one_index = 0;

  // N rows of dimension() coordinates each: row i is the normal form of x_i.
  std::vector<std::vector<rational_function<P>>> variable_coordinates;

  static constexpr std::size_t num_unknowns = N;
  static constexpr std::size_t num_parameters = P;

  std::size_t dimension() const noexcept { return quotient.dimension(); }

  // The invariants the emitter would otherwise have to trust. Checked at the
  // top of emit(), because a header generated from an inconsistent solution
  // compiles perfectly and answers wrongly, which is the worst failure this
  // code can have.
  bool is_well_formed() const {
    if (!quotient.is_zero_dimensional || quotient.dimension() == 0) {
      return false;
    }
    if (action.size() != N) {
      return false;
    }
    for (const auto& m : action) {
      if (m.dimension != quotient.dimension()) {
        return false;
      }
      if (m.entries.size() != m.dimension * m.dimension) {
        return false;
      }
    }
    if (one_index >= quotient.dimension()) {
      return false;
    }
    if (quotient.monomials[one_index].degree() != 0) {
      return false;  // one_index must point at 1, not at some other monomial
    }
    if (variable_coordinates.size() != N) {
      return false;
    }
    for (const auto& row : variable_coordinates) {
      if (row.size() != quotient.dimension()) {
        return false;
      }
    }
    return true;
  }
};

}  // namespace codegen
}  // namespace varietas

#endif
