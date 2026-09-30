#ifndef VARIETAS_CODEGEN_PRIME_FIELD_HPP
#define VARIETAS_CODEGEN_PRIME_FIELD_HPP

#include <cstdint>

#include <gmpxx.h>

#include "varietas/codegen/rational.hpp"
#include "varietas/core/config.hpp"

namespace varietas {

// Z/p for a word-sized prime chosen at run time, as a coefficient field the
// rest of the library can compute over.
//
// This is what lets the parametric solve be reconstructed rather than carried
// out symbolically (varietas/ik/reconstructed_ik.hpp). There, the kinematic
// ideal is solved at many poses, each a Grobner basis with no parameters at
// all, and doing those over Z/p rather than over Q keeps every coefficient one
// machine word, whatever the arm's rationals look like.
//
// The prime is a property of the computation rather than of the type, since
// the reconstruction walks through several of them, and making it a template
// parameter would instantiate the whole of Buchberger once per prime. So it
// is selected per thread, and every residue in play belongs to the prime that
// was selected when it was made. Mixing residues across a change of prime is a
// bug the type cannot catch, and the only code that changes it is the
// reconstruction, which does so between samples and never inside one.
class residue {
 public:
  using word = std::uint64_t;

  // Below 2^31, so that a product of two residues fits in 63 bits.
  static void select(word p) noexcept {
    VARIETAS_ASSERT(p > 2 && p < (word{1} << 31));
    prime_ = p;
  }
  static word prime() noexcept { return prime_; }

  // Set when a rational was mapped here whose denominator the prime divides.
  // The image is then undefined, the residue made from it is meaningless, and
  // whatever used it has to be discarded. The flag is how the caller learns
  // that it happened, and it stays set until cleared.
  static bool undefined_image() noexcept { return undefined_; }
  static void clear_undefined_image() noexcept { undefined_ = false; }

  residue() noexcept : value_(0) {}

  // NOLINTNEXTLINE(google-explicit-constructor): integers are field elements,
  // and the generic code writes Coeff(0) and Coeff(1).
  residue(long long v) noexcept
      : value_(static_cast<word>(((v % static_cast<long long>(prime_)) +
                                  static_cast<long long>(prime_)) %
                                 static_cast<long long>(prime_))) {}

  // NOLINTNEXTLINE(google-explicit-constructor): the embedding of an exact
  // chain, whose coefficients are rationals, is by conversion.
  residue(const rational& q) : value_(0) {
    const unsigned long p = static_cast<unsigned long>(prime_);
    const word numerator = mpz_fdiv_ui(q.get_num().get_mpz_t(), p);
    const word denominator = mpz_fdiv_ui(q.get_den().get_mpz_t(), p);
    if (denominator == 0) {
      undefined_ = true;
      return;
    }
    value_ = numerator * inverse_of(denominator) % prime_;
  }

  static residue raw(word v) noexcept {
    residue r;
    r.value_ = v;
    return r;
  }

  word value() const noexcept { return value_; }

  friend residue operator+(residue a, residue b) noexcept {
    const word s = a.value_ + b.value_;
    return raw(s >= prime_ ? s - prime_ : s);
  }
  friend residue operator-(residue a, residue b) noexcept {
    return raw(a.value_ >= b.value_ ? a.value_ - b.value_ : a.value_ + prime_ - b.value_);
  }
  friend residue operator*(residue a, residue b) noexcept {
    return raw(a.value_ * b.value_ % prime_);
  }
  friend residue operator/(residue a, residue b) noexcept {
    return raw(a.value_ * inverse_of(b.value_) % prime_);
  }
  friend residue operator-(residue a) noexcept { return raw(a.value_ == 0 ? 0 : prime_ - a.value_); }

  residue& operator+=(residue o) noexcept { return *this = *this + o; }
  residue& operator-=(residue o) noexcept { return *this = *this - o; }
  residue& operator*=(residue o) noexcept { return *this = *this * o; }
  residue& operator/=(residue o) noexcept { return *this = *this / o; }

  friend bool operator==(residue a, residue b) noexcept { return a.value_ == b.value_; }
  friend bool operator!=(residue a, residue b) noexcept { return a.value_ != b.value_; }

  // Extended Euclid.
  static word inverse_of(word a) noexcept {
    VARIETAS_ASSERT(a % prime_ != 0);
    std::int64_t r0 = static_cast<std::int64_t>(prime_);
    std::int64_t r1 = static_cast<std::int64_t>(a % prime_);
    std::int64_t t0 = 0;
    std::int64_t t1 = 1;
    while (r1 != 0) {
      const std::int64_t q = r0 / r1;
      const std::int64_t r2 = r0 - q * r1;
      r0 = r1;
      r1 = r2;
      const std::int64_t t2 = t0 - q * t1;
      t0 = t1;
      t1 = t2;
    }
    return static_cast<word>(t0 < 0 ? t0 + static_cast<std::int64_t>(prime_) : t0);
  }

 private:
  static inline thread_local word prime_ = 2147483647;
  static inline thread_local bool undefined_ = false;

  word value_;
};

template <>
struct coefficient_traits<residue> {
  static constexpr bool is_exact = true;

  static residue zero() noexcept { return residue(); }
  static residue one() noexcept { return residue::raw(1); }
  static bool is_zero(residue c) noexcept { return c.value() == 0; }
  static residue inverse(residue c) noexcept {
    VARIETAS_ASSERT(c.value() != 0);
    return residue::raw(residue::inverse_of(c.value()));
  }
  static residue negate(residue c) noexcept { return -c; }

  // There is no numerical side to cross to from a prime field, but the chain's
  // validation asks for a double to compare against zero. The representative
  // in [0, p) answers the only question it asks, whether a value is zero.
  static double to_double(residue c) noexcept { return static_cast<double>(c.value()); }
  static residue from_double(double c) noexcept { return residue(static_cast<long long>(c)); }
};

}  // namespace varietas

#endif
