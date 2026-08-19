#ifndef VARIETAS_CODEGEN_RATIONAL_FUNCTION_HPP
#define VARIETAS_CODEGEN_RATIONAL_FUNCTION_HPP

#include <array>
#include <cstddef>
#include <string>
#include <vector>

#include "varietas/codegen/rational.hpp"
#include "varietas/core/config.hpp"
#include "varietas/core/order/grevlex.hpp"
#include "varietas/core/polynomial.hpp"

namespace varietas {

// The coefficient field the emitter computes over: rational functions of the
// pose, with rational coefficients.
//
// The point is a fact about how the pose enters the problem. A residual is
// built as translation(i) - denominator * target[i], so the target appears in
// the generators only as a coefficient, and affinely at that. Fixing it to a
// number therefore produces a Gröbner basis that answers exactly one pose,
// which is no use to a robot. Leaving it symbolic — adjoining x, y, z to the
// coefficient field rather than to the polynomial ring — produces a basis whose
// coefficients are rational functions of the pose, and the action matrix
// assembled from it has entries that are rational functions too. That matrix is
// what the generated header stores, and evaluating it at a pose is what the
// runtime does.
//
// Nothing in varietas_core changes to accommodate this. The ideal, quotient and
// solving layers are templated on the coefficient type and ask it only for the
// field operations, so the same Buchberger implementation that runs over Q runs
// over Q(x, y, z) — the second time this parameterisation has paid for itself.
//
// P is the number of pose parameters: 3 for a position target, 6 for a pose.
template <std::size_t P>
class rational_function {
 public:
  using parameter_order = grevlex;
  using parameter_polynomial = polynomial<rational, P, parameter_order>;

  static constexpr std::size_t num_parameters = P;

  rational_function() : numerator_(), denominator_(parameter_polynomial::constant(rational(1))) {}

  // NOLINTNEXTLINE(google-explicit-constructor) — a rational is a constant
  // function, and the implicit conversion is what lets the existing generic
  // code build coefficients without knowing the field.
  rational_function(const rational& c)
      : numerator_(parameter_polynomial::constant(c)),
        denominator_(parameter_polynomial::constant(rational(1))) {}

  explicit rational_function(const parameter_polynomial& p)
      : numerator_(p), denominator_(parameter_polynomial::constant(rational(1))) {}

  rational_function(const parameter_polynomial& n, const parameter_polynomial& d)
      : numerator_(n), denominator_(d) {
    VARIETAS_ASSERT(!d.is_zero());
    normalize();
  }

  // The parameter itself: x, y, z as elements of the coefficient field.
  static rational_function parameter(std::size_t i) {
    VARIETAS_ASSERT(i < P);
    return rational_function(parameter_polynomial::variable(i));
  }

  const parameter_polynomial& numerator() const noexcept { return numerator_; }
  const parameter_polynomial& denominator() const noexcept { return denominator_; }

  bool is_zero() const noexcept { return numerator_.is_zero(); }

  // True when the function does not depend on the pose at all, which is the
  // only case in which it has a value without one being supplied.
  bool is_constant() const noexcept {
    return numerator_.degree() == 0 && denominator_.degree() == 0;
  }

  rational constant_value() const {
    VARIETAS_ASSERT(is_constant());
    if (numerator_.is_zero()) {
      return rational(0);
    }
    return numerator_.leading_coefficient() / denominator_.leading_coefficient();
  }

  // The value at a pose. This is what the runtime does, and what the tests use
  // to check a parametric basis against the numeric one.
  rational evaluate(const std::array<rational, P>& at) const {
    const rational d = denominator_.evaluate(at);
    VARIETAS_ASSERT(!coefficient_traits<rational>::is_zero(d));
    return numerator_.evaluate(at) / d;
  }

  rational_function& operator+=(const rational_function& other) {
    numerator_ = numerator_ * other.denominator_ + other.numerator_ * denominator_;
    denominator_ = denominator_ * other.denominator_;
    normalize();
    return *this;
  }

  rational_function& operator-=(const rational_function& other) {
    numerator_ = numerator_ * other.denominator_ - other.numerator_ * denominator_;
    denominator_ = denominator_ * other.denominator_;
    normalize();
    return *this;
  }

  rational_function& operator*=(const rational_function& other) {
    numerator_ = numerator_ * other.numerator_;
    denominator_ = denominator_ * other.denominator_;
    normalize();
    return *this;
  }

  rational_function& operator/=(const rational_function& other) {
    VARIETAS_ASSERT(!other.is_zero());
    numerator_ = numerator_ * other.denominator_;
    denominator_ = denominator_ * other.numerator_;
    normalize();
    return *this;
  }

  friend rational_function operator+(rational_function a, const rational_function& b) {
    return a += b;
  }
  friend rational_function operator-(rational_function a, const rational_function& b) {
    return a -= b;
  }
  friend rational_function operator*(rational_function a, const rational_function& b) {
    return a *= b;
  }
  friend rational_function operator/(rational_function a, const rational_function& b) {
    return a /= b;
  }

  friend rational_function operator-(const rational_function& a) {
    rational_function r;
    r.numerator_ = -a.numerator_;
    r.denominator_ = a.denominator_;
    return r;
  }

  // Equality is of the function, not of the written form: a/b == c/d exactly
  // when ad == cb. Normalisation makes the representation canonical up to the
  // common factors it fails to cancel, so the cross product is what decides.
  friend bool operator==(const rational_function& a, const rational_function& b) {
    return a.numerator_ * b.denominator_ == b.numerator_ * a.denominator_;
  }
  friend bool operator!=(const rational_function& a, const rational_function& b) {
    return !(a == b);
  }

 private:
  // Reduce the representation as far as cheap operations allow: the rational
  // content, the monomial content, and a sign convention on the denominator.
  //
  // What this deliberately does not do is a full multivariate GCD. The
  // consequence is expression swell rather than incorrectness — a common factor
  // that is a genuine polynomial survives — and whether it is tolerable is a
  // question about the arms actually put through the emitter, not one that can
  // be settled in the abstract. It is measured rather than assumed.
  void normalize() {
    if (numerator_.is_zero()) {
      denominator_ = parameter_polynomial::constant(rational(1));
      return;
    }

    remove_monomial_content(numerator_, denominator_);
    remove_rational_content(numerator_, denominator_);

    // A negative leading coefficient on the denominator is moved to the
    // numerator, so that equal functions are written the same way.
    if (sgn(denominator_.leading_coefficient()) < 0) {
      numerator_ = -numerator_;
      denominator_ = -denominator_;
    }
  }

  // The largest monomial dividing every term of both, divided out of each.
  static void remove_monomial_content(parameter_polynomial& n, parameter_polynomial& d) {
    monomial<P> common = n.terms().front().mon;
    for (const auto& t : n.terms()) {
      common = monomial<P>::gcd(common, t.mon);
    }
    for (const auto& t : d.terms()) {
      common = monomial<P>::gcd(common, t.mon);
    }
    if (common.is_one()) {
      return;
    }
    n = divide_by_monomial(n, common);
    d = divide_by_monomial(d, common);
  }

  static parameter_polynomial divide_by_monomial(const parameter_polynomial& p,
                                                 const monomial<P>& m) {
    std::vector<typename parameter_polynomial::term> terms;
    terms.reserve(p.size());
    for (const auto& t : p.terms()) {
      VARIETAS_ASSERT(monomial<P>::divides(m, t.mon));
      terms.push_back({monomial<P>::divide(t.mon, m), t.coeff});
    }
    return parameter_polynomial(std::move(terms));
  }

  // Clear denominators and divide out the integer content of both, which keeps
  // the coefficients from growing without bound under repeated arithmetic.
  static void remove_rational_content(parameter_polynomial& n, parameter_polynomial& d) {
    const rational scale = content(n) / content(d);
    if (coefficient_traits<rational>::is_zero(scale)) {
      return;
    }
    n = n / content(n);
    d = d / content(d);
    n = n * scale;
    // n and d are now integral and primitive up to the shared scale; move the
    // scale's own content back out so neither side carries it twice.
    const rational shared = content(n);
    if (!coefficient_traits<rational>::is_zero(shared)) {
      n = n / shared;
      d = d / shared;
    }
  }

  // The rational content: the gcd of the numerators over the lcm of the
  // denominators, so that dividing by it makes the polynomial integral and
  // primitive.
  static rational content(const parameter_polynomial& p) {
    if (p.is_zero()) {
      return rational(1);
    }
    mpz_class numerator_gcd = 0;
    mpz_class denominator_lcm = 1;
    for (const auto& t : p.terms()) {
      mpz_gcd(numerator_gcd.get_mpz_t(), numerator_gcd.get_mpz_t(),
              t.coeff.get_num().get_mpz_t());
      mpz_lcm(denominator_lcm.get_mpz_t(), denominator_lcm.get_mpz_t(),
              t.coeff.get_den().get_mpz_t());
    }
    rational c(numerator_gcd, denominator_lcm);
    c.canonicalize();
    return c;
  }

  parameter_polynomial numerator_;
  parameter_polynomial denominator_;
};

// The field operations the generic algorithms consult. Everything varietas_core
// asks of a coefficient is here, which is why nothing there has to change.
template <std::size_t P>
struct coefficient_traits<rational_function<P>> {
  using value_type = rational_function<P>;

  static constexpr bool is_exact = true;

  static value_type zero() { return value_type(rational(0)); }
  static value_type one() { return value_type(rational(1)); }

  static bool is_zero(const value_type& c) { return c.is_zero(); }

  static value_type inverse(const value_type& c) {
    VARIETAS_ASSERT(!is_zero(c));
    return value_type(c.denominator(), c.numerator());
  }

  static value_type negate(const value_type& c) { return -c; }

  // A rational function has no value until a pose is supplied, so the bridge to
  // the numerical side exists only for the constants. Crossing it on anything
  // else is a bug in the caller: the emitter writes the expression out and lets
  // the runtime evaluate it, which is the entire point of computing here.
  static double to_double(const value_type& c) {
    VARIETAS_ASSERT(c.is_constant());
    return c.constant_value().get_d();
  }

  static value_type from_double(double c) { return value_type(rational(c)); }
};

}  // namespace varietas

#endif
