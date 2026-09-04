#ifndef VARIETAS_CODEGEN_TEST_PLANAR_FIXTURE_HPP
#define VARIETAS_CODEGEN_TEST_PLANAR_FIXTURE_HPP

#include <cstddef>
#include <cstdint>
#include <vector>

#include "varietas/codegen/parametric_solution.hpp"
#include "varietas/codegen/rational_function.hpp"
#include "varietas/core/ideal/ideal.hpp"
#include "varietas/core/monomial.hpp"
#include "varietas/core/order/grevlex.hpp"

namespace varietas_test {

// A parametric system small enough to solve by hand and awkward enough to be
// worth generating code for.
//
//   u^2 - x
//   v   - y u
//
// Under grevlex the leading term of v - y u is u, not v, so it is u that gets
// eliminated: the reduced basis is {u - v/y, v^2 - x y^2} and the standard
// monomials are {1, v}. dim_k A = 2, and since u = v/y and v^2 = x y^2,
//
//   u.1 = v/y      u.v = v^2/y = x y
//   v.1 = v        v.v = v^2   = x y^2
//
// which read off in the basis {1, v} as
//
//   M_u = [ 0    xy ]      M_v = [ 0  xy^2 ]
//         [ 1/y  0  ]            [ 1  0    ]
//
// Both carry the pose, so a generated header for this system has to evaluate
// something rather than look a constant up, and the closed form above is what
// the generated code is checked against.
//
// The 1/y is not an accident of this example and is the reason emit() writes a
// guard at all. Solving for u in terms of v divides by y, so the parametric
// basis describes the system only where y != 0; y = 0 is a locus the one basis
// cannot cover, and a generated header has to say so rather than return
// infinities.
inline varietas::codegen::parametric_solution<2, 2> planar_solution() {
  using varietas::grevlex;
  using field = varietas::rational_function<2>;
  using parametric = varietas::polynomial<field, 2, grevlex>;

  const auto u = parametric::variable(0);
  const auto v = parametric::variable(1);
  const auto x = field::parameter(0);
  const auto y = field::parameter(1);

  const std::vector<parametric> generators{
      u * u - parametric::constant(x),
      v - parametric::constant(y) * u,
  };

  const varietas::ideal<field, 2, grevlex> ideal(generators);

  varietas::codegen::parametric_solution<2, 2> solution;
  solution.order = varietas::order_id::grevlex;
  solution.unknown_names = {"u", "v"};
  solution.parameter_names = {"x", "y"};
  solution.quotient = ideal.quotient();
  solution.one_index = solution.quotient.index_of(varietas::monomial<2>::one());
  for (std::size_t variable = 0; variable < 2; ++variable) {
    solution.action.push_back(varietas::codegen::parametric_action_matrix<2, 2, grevlex>(
        variable, ideal.basis(), solution.quotient));
  }
  solution.variable_coordinates =
      varietas::codegen::parametric_variable_coordinates<2, 2, grevlex>(ideal.basis(),
                                                                        solution.quotient);
  return solution;
}

}  // namespace varietas_test

#endif
