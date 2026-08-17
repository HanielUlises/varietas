#ifndef VARIETAS_URDF_RATIONAL_APPROXIMATION_HPP
#define VARIETAS_URDF_RATIONAL_APPROXIMATION_HPP

#include <algorithm>
#include <cmath>
#include <cstddef>

#include <gmpxx.h>

#include "varietas/codegen/rational.hpp"
#include "varietas/kinematics/rigid_transform.hpp"

namespace varietas {
namespace urdf_import {

// Recovering exact geometry from a file that does not contain any.
//
// A URDF stores its angles as decimal literals, and the angles that matter are
// almost always multiples of a right angle: the KUKA iiwa writes pi/2 as
// 1.57079632679, which is neither pi/2 nor a number with any exact meaning of
// its own. Reading that literal as a rational and handing it to Buchberger
// produces an exact answer about a robot whose axes are misaligned by 1e-12
// radians — a robot nobody built. What is wanted is the number the author
// meant, together with a bound on the distance to the number they wrote.
//
// The recovery works because of a fact about quaternions rather than about
// angles. A rotation by a right angle has quaternion (cos 45, sin 45, 0, 0),
// whose entries are irrational; but a quaternion is homogeneous, and dividing
// through by any nonzero entry leaves (1, 1, 0, 0), which is integral, and
// which rotation_from_quaternion turns into an exactly orthogonal rational
// matrix without ever normalising or taking a square root. So the search is
// projective: divide by the largest entry, then approximate each of the four
// by a rational of bounded denominator. Every multiple of a right angle is
// recovered exactly this way; a genuinely oblique angle is not, and is reported
// with the deviation it would introduce rather than silently accepted.

// The closest rational to value whose denominator does not exceed the bound, by
// the standard continued fraction argument: the convergents of the expansion
// are the best approximations of their denominator, and the best approximation
// under a bound is either the last convergent within it or the semiconvergent
// that just reaches it.
inline rational limit_denominator(const rational& value, const mpz_class& bound) {
  VARIETAS_ASSERT(bound >= 1);

  mpz_class n = value.get_num();
  mpz_class d = value.get_den();
  if (d <= bound) {
    return value;
  }

  const bool negative = n < 0;
  if (negative) {
    n = -n;
  }

  mpz_class p0 = 0, q0 = 1, p1 = 1, q1 = 0;
  while (true) {
    const mpz_class a = n / d;
    const mpz_class q2 = q0 + a * q1;
    if (q2 > bound) {
      break;
    }
    mpz_class p2 = p0 + a * p1;
    p0 = p1;
    q0 = q1;
    p1 = p2;
    q1 = q2;

    const mpz_class remainder = n - a * d;
    n = d;
    d = remainder;
    if (d == 0) {
      break;
    }
  }

  // The two candidates: the semiconvergent that just fits under the bound, and
  // the last convergent itself.
  const mpz_class k = (bound - q0) / q1;
  rational lower(p0 + k * p1, q0 + k * q1);
  rational upper(p1, q1);
  lower.canonicalize();
  upper.canonicalize();

  rational target(negative ? -value : value);
  rational best = abs(upper - target) <= abs(lower - target) ? upper : lower;
  return negative ? rational(-best) : best;
}

inline rational rationalize(double value, long max_denominator = 1000000) {
  VARIETAS_ASSERT(max_denominator >= 1);
  // rational(value) is the binary value exactly; the approximation then walks
  // back from it to the decimal, or the fraction, the author is likely to have
  // meant.
  return limit_denominator(rational(value), mpz_class(max_denominator));
}

// The double nearest a rational.
//
// mpq_class::get_d truncates towards zero rather than rounding to nearest, so
// it can land a whole unit in the last place away from the value and always
// errs in the same direction. That is too coarse to decide whether a recovered
// number is the one the file held: 409/2000 truncates to a double one ulp below
// the double that parsing "0.2045" produces, and the recovery would be reported
// as having moved a length that it did not touch. Truncating and then examining
// the neighbour costs one exact comparison and removes both the bias and the
// ambiguity.
inline double nearest_double(const rational& value) {
  const double truncated = value.get_d();
  if (!std::isfinite(truncated)) {
    return truncated;
  }
  const double neighbour =
      std::nextafter(truncated, value >= rational(0) ? HUGE_VAL : -HUGE_VAL);
  if (!std::isfinite(neighbour)) {
    return truncated;
  }
  const rational to_truncated = abs(rational(truncated) - value);
  const rational to_neighbour = abs(rational(neighbour) - value);
  return to_neighbour < to_truncated ? neighbour : truncated;
}

// A scalar recovered exactly enough, with the distance to what the file said.
//
// The two fields answer different questions. The deviation is how far the
// recovered number is from the double in the file; round_trips says whether
// that gap is below the resolution of a double at all, so that the recovery is
// not merely close to the file but indistinguishable from it. A length written
// as 0.1575 recovers as 63/400 and round trips — the residue of 3e-17 is the
// file's own decimal-to-binary rounding, which the recovery undoes rather than
// commits. A truncated pi does not round trip, and should not: there the
// recovery deliberately moves the robot onto the right angle that was meant.
struct scalar_snap {
  rational value;
  double deviation = 0.0;
  bool round_trips = false;
};

inline scalar_snap snap_scalar(double value, long max_denominator = 1000000) {
  scalar_snap snap;
  snap.value = rationalize(value, max_denominator);
  const double recovered = nearest_double(snap.value);
  snap.deviation = std::fabs(recovered - value);
  snap.round_trips = recovered == value;
  return snap;
}

// The angle, in radians, of the rotation carrying a to b, which is the
// geometrically meaningful measure of how far a snapped rotation has moved.
//
// Not by the trace: recovering the angle as acos((tr - 1) / 2) is useless
// exactly where this quantity matters. Near the identity the trace differs from
// three by the square of the angle, so a deviation of 1e-12 radians moves the
// trace by 1e-24, which vanishes into the rounding of a trace of order one, and
// the answer comes back a flat zero — a rotation reported as recovered exactly
// when it was not. The chordal form has no such cancellation: for the relative
// rotation of angle t, the Frobenius distance ||a - b|| is 2 sqrt(2) sin(t/2)
// exactly, and inverting that is accurate all the way down and still correct
// out to a half turn.
inline double rotation_distance(const matrix3<double>& a, const matrix3<double>& b) {
  double sum_of_squares = 0.0;
  for (std::size_t i = 0; i < 3; ++i) {
    for (std::size_t j = 0; j < 3; ++j) {
      const double difference = a(i, j) - b(i, j);
      sum_of_squares += difference * difference;
    }
  }
  const double half_sine =
      std::min(1.0, std::sqrt(sum_of_squares) / (2.0 * std::sqrt(2.0)));
  return 2.0 * std::asin(half_sine);
}

// A rotation recovered as an exact rational quaternion, with the angle between
// it and the one the file described.
struct rotation_snap {
  rational w, x, y, z;
  double deviation_radians = 0.0;
  // True when the recovered rotation is the stated one to the last bit of a
  // double, so that nothing about the robot has moved. False for a truncated
  // right angle, where moving it is the whole point.
  bool round_trips = false;

  matrix3<rational> exact() const { return rotation_from_quaternion<rational>(w, x, y, z); }
};

// Projective rationalisation of a quaternion, in the (x, y, z, w) order urdfdom
// reports. The result is exactly orthogonal of determinant one whatever the
// input, since that holds for every nonzero rational quaternion; the quantity
// that has to be checked is therefore not orthogonality but the deviation, and
// that is what is returned alongside.
inline rotation_snap snap_rotation(double qx, double qy, double qz, double qw,
                                   long max_denominator = 1000000) {
  const double components[4] = {qw, qx, qy, qz};

  std::size_t pivot = 0;
  for (std::size_t i = 1; i < 4; ++i) {
    if (std::fabs(components[i]) > std::fabs(components[pivot])) {
      pivot = i;
    }
  }
  // A quaternion of norm one always has an entry of magnitude at least 1/2, so
  // the pivot cannot be small unless the input was not a rotation at all.
  VARIETAS_ASSERT(std::fabs(components[pivot]) > 1e-6);

  rotation_snap snap;
  const double scale = components[pivot];
  rational* out[4] = {&snap.w, &snap.x, &snap.y, &snap.z};
  for (std::size_t i = 0; i < 4; ++i) {
    *out[i] = (i == pivot) ? rational(1)
                           : rationalize(components[i] / scale, max_denominator);
  }

  // The deviation is measured between rotations, not between quaternions, so
  // that the sign and the scale of the representative cannot affect it.
  matrix3<double> recovered;
  const matrix3<rational> exact = snap.exact();
  for (std::size_t i = 0; i < 3; ++i) {
    for (std::size_t j = 0; j < 3; ++j) {
      recovered(i, j) = nearest_double(exact(i, j));
    }
  }
  const matrix3<double> stated =
      rotation_from_quaternion<double>(qw, qx, qy, qz);
  snap.deviation_radians = rotation_distance(stated, recovered);
  snap.round_trips = recovered == stated;
  return snap;
}

}  // namespace urdf_import
}  // namespace varietas

#endif
