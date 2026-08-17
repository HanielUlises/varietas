#ifndef VARIETAS_CODEGEN_RATIONAL_HPP
#define VARIETAS_CODEGEN_RATIONAL_HPP

#include <cstdint>
#include <string>

#include <gmpxx.h>

#include "varietas/core/config.hpp"

namespace varietas {

// Exact coefficient field for the offline pipeline.
//
// Buchberger's algorithm decides termination by asking whether an S-polynomial
// reduces to zero. Over a floating point field that question is undecidable in
// practice: cancellation leaves a residue of the order of the rounding error,
// the pair is never discarded, the basis accretes spurious elements, and the
// finiteness verdict read off the leading terms is computed from a basis of a
// different ideal. The coefficient blowup that makes this unavoidable is also
// what rules out fixed-width integers, so the numerator and denominator are
// arbitrary precision.
//
// mpq_class keeps its value in lowest terms with a positive denominator after
// every operation, so equality with zero is a structural test and no
// normalisation step is needed here.
using rational = mpq_class;

// Exact constructors. Prefer these to conversion from double: from_double is
// exact but converts the binary value actually stored, so from_double(0.1) is
// 3602879701896397/36028797018963968 and not 1/10.
inline rational make_rational(std::int64_t numerator, std::int64_t denominator = 1) {
  VARIETAS_ASSERT(denominator != 0);
  rational q(numerator, denominator);
  q.canonicalize();
  return q;
}

// Parses "3", "-7/2" and the like. GMP requires the value be canonicalised
// after construction from a string, which it does not do itself.
inline rational rational_from_string(const std::string& text, int base = 10) {
  rational q(text, base);
  q.canonicalize();
  return q;
}

template <>
struct coefficient_traits<rational> {
  // Consulted by algorithms that are only sound over an exact field, and by
  // polynomial::prune, which refuses to compile against an exact coefficient.
  static constexpr bool is_exact = true;

  static rational zero() { return rational(0, 1); }
  static rational one() { return rational(1, 1); }

  // mpq_class is canonical, so the numerator alone decides.
  static bool is_zero(const rational& c) { return sgn(c) == 0; }

  static rational inverse(const rational& c) {
    VARIETAS_ASSERT(!is_zero(c));
    rational r;
    mpq_inv(r.get_mpq_t(), c.get_mpq_t());
    return r;
  }

  static rational negate(const rational& c) { return rational(-c); }

  // The bridge to the numerical side. Action matrices and the spectral solver
  // cross it deliberately: the Gröbner basis is computed exactly, and only the
  // eigendecomposition that recovers the points works in floating point.
  static double to_double(const rational& c) { return c.get_d(); }

  static rational from_double(double c) { return rational(c); }
};

}  // namespace varietas

#endif
