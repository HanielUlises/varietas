#ifndef VARIETAS_CORE_SOLVE_SPECTRAL_HPP
#define VARIETAS_CORE_SOLVE_SPECTRAL_HPP

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <vector>

#include <Eigen/Core>
#include <Eigen/Eigenvalues>

#include "varietas/core/config.hpp"
#include "varietas/core/monomial.hpp"
#include "varietas/core/polynomial.hpp"
#include "varietas/core/quotient/action_matrix.hpp"
#include "varietas/core/quotient/quotient_basis.hpp"

namespace varietas {

// Why a solve may fail to produce points. Failure is always reported as a
// hypothesis that was violated, never as a silently truncated solution set.
enum class solve_status {
  ok = 0,
  // The variety is empty over the algebraic closure: the Gröbner basis is {1}.
  empty_variety,
  // The quotient is infinite dimensional, so the system has infinitely many
  // solutions and no action matrix exists.
  positive_dimensional,
  // The standard monomials do not contain 1 and every variable, so the
  // coordinates of a point cannot be read off an eigenvector directly.
  basis_lacks_linear_monomials,
  // The eigendecomposition of the action matrix did not converge.
  eigensolver_failed,
};

constexpr const char* to_string(solve_status status) noexcept {
  return status == solve_status::ok                    ? "ok"
         : status == solve_status::empty_variety       ? "empty_variety"
         : status == solve_status::positive_dimensional ? "positive_dimensional"
         : status == solve_status::basis_lacks_linear_monomials
             ? "basis_lacks_linear_monomials"
         : status == solve_status::eigensolver_failed ? "eigensolver_failed"
                                                      : "unknown";
}

template <std::size_t N>
struct solution_set {
  solve_status status = solve_status::ok;
  // All points of the variety over the complex numbers, counted once each.
  std::vector<std::array<std::complex<double>, N>> points;

  bool ok() const noexcept { return status == solve_status::ok; }

  // The points whose imaginary part is negligible, which are the only ones a
  // manipulator can be commanded to.
  std::vector<std::array<double, N>> real_points(double tolerance = 1e-8) const {
    std::vector<std::array<double, N>> real;
    for (const auto& p : points) {
      bool is_real = true;
      for (std::size_t i = 0; i < N && is_real; ++i) {
        is_real = std::abs(p[i].imag()) <= tolerance;
      }
      if (is_real) {
        std::array<double, N> q{};
        for (std::size_t i = 0; i < N; ++i) {
          q[i] = p[i].real();
        }
        real.push_back(q);
      }
    }
    return real;
  }
};

// Solves a zero-dimensional system by eigendecomposition of an action matrix.
//
// The construction is the classical one of Stetter and Möller. For a generic
// linear form u, the eigenvalues of the multiplication operator m_u are the
// values of u at the points of the variety, and the left eigenvectors are the
// evaluation functionals at those points. Since a left eigenvector v satisfies
// v_m = m(p) v_1 for every standard monomial m, the coordinates of p are read
// off the entries indexed by 1 and by the variables, after normalising by the
// entry indexed by 1.
//
// The linear form is taken with fixed pseudo-random coefficients rather than a
// single variable, because a variable that takes the same value at two distinct
// points produces a repeated eigenvalue and an ambiguous eigenvector. Genericity
// fails only on a proper algebraic subset of the coefficient space; the residual
// check that callers can run on the returned points is what detects it.
//
// The basis must be a reduced Gröbner basis of the ideal under Order.
template <class Coeff, std::size_t N, class Order>
solution_set<N> solve_zero_dimensional(const std::vector<polynomial<Coeff, N, Order>>& basis,
                                       double tolerance = 1e-9) {
  using poly = polynomial<Coeff, N, Order>;

  solution_set<N> result;

  const quotient_basis<N> quotient = standard_monomials(basis);
  if (!quotient.is_zero_dimensional) {
    result.status = solve_status::positive_dimensional;
    return result;
  }
  if (quotient.dimension() == 0) {
    result.status = solve_status::empty_variety;
    return result;
  }

  const std::size_t one_index = quotient.index_of(monomial<N>::one());
  std::array<std::size_t, N> variable_index{};
  for (std::size_t i = 0; i < N; ++i) {
    variable_index[i] = quotient.index_of(monomial<N>::variable(i));
    if (variable_index[i] >= quotient.dimension()) {
      result.status = solve_status::basis_lacks_linear_monomials;
      return result;
    }
  }
  if (one_index >= quotient.dimension()) {
    result.status = solve_status::basis_lacks_linear_monomials;
    return result;
  }

  // A fixed generic linear form. The coefficients are deterministic so that
  // two runs on the same ideal return the same points in the same arrangement.
  poly separating;
  for (std::size_t i = 0; i < N; ++i) {
    const double c = 1.0 + 0.37 * static_cast<double>(i) + 0.11 * static_cast<double>(i * i);
    separating += poly::from_monomial(monomial<N>::variable(i),
                                      coefficient_traits<Coeff>::from_double(c));
  }

  const action_matrix_type matrix = action_matrix(separating, basis, quotient);

  // Left eigenvectors of M are the right eigenvectors of its transpose.
  Eigen::EigenSolver<action_matrix_type> solver(matrix.transpose(), true);
  if (solver.info() != Eigen::Success) {
    result.status = solve_status::eigensolver_failed;
    return result;
  }

  const auto vectors = solver.eigenvectors();
  result.points.reserve(static_cast<std::size_t>(vectors.cols()));

  for (Eigen::Index k = 0; k < vectors.cols(); ++k) {
    const auto v = vectors.col(k);
    const std::complex<double> scale = v(static_cast<Eigen::Index>(one_index));
    if (std::abs(scale) <= tolerance) {
      // The functional does not evaluate 1 to a nonzero value, which happens
      // when the eigenvalue is repeated and the eigenvector is an arbitrary
      // element of an eigenspace of dimension greater than one. Such a vector
      // carries no point, so it is dropped rather than reported as a solution.
      continue;
    }

    std::array<std::complex<double>, N> point{};
    for (std::size_t i = 0; i < N; ++i) {
      point[i] = v(static_cast<Eigen::Index>(variable_index[i])) / scale;
    }
    result.points.push_back(point);
  }

  return result;
}

// Largest modulus of the generators evaluated at a candidate point. Callers
// use it to certify the output of the spectral solve, which is the only step
// of the pipeline where floating point enters.
template <class Coeff, std::size_t N, class Order>
double residual(const std::vector<polynomial<Coeff, N, Order>>& generators,
                const std::array<std::complex<double>, N>& point) {
  using traits = coefficient_traits<Coeff>;

  double worst = 0.0;
  for (const auto& g : generators) {
    std::complex<double> value{0.0, 0.0};
    for (const auto& t : g.terms()) {
      std::complex<double> monomial_value{1.0, 0.0};
      for (std::size_t i = 0; i < N; ++i) {
        for (std::size_t k = 0; k < t.mon[i]; ++k) {
          monomial_value *= point[i];
        }
      }
      value += traits::to_double(t.coeff) * monomial_value;
    }
    worst = std::max(worst, std::abs(value));
  }
  return worst;
}

}  // namespace varietas

#endif
