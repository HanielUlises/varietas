#ifndef VARIETAS_CODEGEN_RATIONAL_FUNCTION_HPP
#define VARIETAS_CODEGEN_RATIONAL_FUNCTION_HPP

#include <array>
#include <cstddef>
#include <string>
#include <vector>

#include "varietas/codegen/rational.hpp"
#include "varietas/core/config.hpp"
#include "varietas/core/gcd.hpp"
#include "varietas/core/ideal/division.hpp"
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
  // Reduce the representation to lowest terms.
  //
  // The cheap reductions come first because they are cheap — the monomial
  // content and the rational content are read straight off the terms — and the
  // polynomial gcd, which is the expensive one, then runs on smaller inputs.
  //
  // Without the gcd this arithmetic is still correct but unusable: the entries
  // of a parametric action matrix run to degree seventy and beyond, when the
  // functions they represent are of degree one or two. Cancellation is not a
  // tidying step here, it is what makes the emitted expressions finite.
  void normalize() {
    if (numerator_.is_zero()) {
      denominator_ = parameter_polynomial::constant(rational(1));
      return;
    }

    remove_monomial_content(numerator_, denominator_);
    remove_rational_content(numerator_, denominator_);
    remove_common_factor(numerator_, denominator_);
    remove_rational_content(numerator_, denominator_);

    // A negative leading coefficient on the denominator is moved to the
    // numerator, so that equal functions are written the same way.
    if (sgn(denominator_.leading_coefficient()) < 0) {
      numerator_ = -numerator_;
      denominator_ = -denominator_;
    }
  }

  // A cheap test that answers "no common factor" or "cannot tell".
  //
  // The subresultant gcd is where a parametric solve spends its time — around
  // eighty-six per cent of it on a three-joint arm, with the cost per call
  // climbing as the parameter polynomials grow — and the great majority of
  // those calls return 1. Deciding coprimality cheaply, and only paying for the
  // exact gcd when the cheap test cannot rule it out, is what keeps that
  // majority cheap.
  //
  // The test specialises every parameter but one at fixed values in a small
  // prime field, leaving two univariate polynomials whose gcd is a Euclidean
  // algorithm on machine integers. A common factor of positive degree in the
  // remaining variable almost always survives that specialisation, so a
  // constant gcd downstairs is strong evidence of coprimality upstairs.
  //
  // Strong evidence and not proof, and that is the point worth being precise
  // about: this test is consulted only to decide whether to *skip* a
  // cancellation. Skipping one that was available leaves the fraction in
  // higher terms than necessary, which costs size and nothing else — the value
  // the fraction stands for is unchanged. So the failure mode is a larger
  // expression, never a wrong one, and the test is allowed to be a heuristic
  // where the arithmetic around it is not.
  //
  // The values are fixed rather than random so that the same inputs always
  // reach the same verdict: a generated header must not depend on when it was
  // generated.
  static bool certainly_coprime(const parameter_polynomial& n, const parameter_polynomial& d) {
    constexpr long long prime = 1000003;  // fits a long long product comfortably

    // Fixed, nonzero and distinct, so that the specialisation is not the zero
    // map on some factor by accident, and so that the same inputs always reach
    // the same verdict — a generated header must not depend on when it was
    // generated.
    std::array<long long, P> at{};
    for (std::size_t v = 0; v < P; ++v) {
      at[v] = static_cast<long long>(2 + 3 * v);
    }

    // Every variable has to be tried, not just the most promising one.
    //
    // Specialising a variable collapses any factor that involves only that
    // variable into a constant, where it is invisible: keeping x and setting
    // y = 5 turns a shared factor of (y + 2) into 7, and the univariate gcd in
    // x then reports coprimality that is not there. But a nonconstant factor
    // has positive degree in some variable, and both polynomials are divisible
    // by it, so that variable is one where both have positive degree — keeping
    // it is what makes the factor visible. Testing every such variable
    // therefore cannot miss a factor for this reason.
    for (std::size_t main = 0; main < P; ++main) {
      const unsigned dn = degree_in(n, main);
      const unsigned dd = degree_in(d, main);
      if (dn == 0 || dd == 0) {
        continue;  // no common factor can have positive degree here
      }

      std::vector<long long> un;
      std::vector<long long> ud;
      if (!specialise(n, main, at, prime, un) || !specialise(d, main, at, prime, ud)) {
        return false;
      }
      // A degree that dropped under specialisation means a leading coefficient
      // vanished there, and the conclusion would not follow. Fall back.
      if (un.size() != dn + 1u || ud.size() != dd + 1u) {
        return false;
      }
      if (univariate_gcd_degree(un, ud, prime) != 0) {
        return false;  // a shared factor survives here; let the exact gcd decide
      }
    }

    // Either every variable that could have carried a shared factor reported
    // none, or no variable has positive degree in both — in which case no
    // nonconstant polynomial divides both and the conclusion needs no test.
    return true;
  }

  static unsigned degree_in(const parameter_polynomial& p, std::size_t v) {
    unsigned d = 0;
    for (const auto& t : p.terms()) {
      const unsigned e = static_cast<unsigned>(t.mon[v]);
      if (e > d) {
        d = e;
      }
    }
    return d;
  }

  static long long power_mod(long long base, unsigned exponent, long long prime) {
    long long result = 1;
    long long factor = ((base % prime) + prime) % prime;
    while (exponent != 0) {
      if ((exponent & 1u) != 0u) {
        result = result * factor % prime;
      }
      factor = factor * factor % prime;
      exponent >>= 1;
    }
    return result;
  }

  // The polynomial with every variable but `main` evaluated, reduced mod prime
  // and returned as coefficients of increasing powers of `main`. False if a
  // coefficient's denominator vanishes mod prime, where the reduction is not
  // defined.
  static bool specialise(const parameter_polynomial& p, std::size_t main,
                         const std::array<long long, P>& at, long long prime,
                         std::vector<long long>& out) {
    out.assign(degree_in(p, main) + 1u, 0);
    for (const auto& t : p.terms()) {
      const mpz_class& numerator = t.coeff.get_num();
      const mpz_class& denominator = t.coeff.get_den();
      // The cast pins down which of GMP's mixed-mode overloads applies; a
      // long long is ambiguous against them.
      const long modulus = static_cast<long>(prime);
      const long long a = mpz_class(numerator % modulus).get_si();
      const long long b = mpz_class(denominator % modulus).get_si();
      const long long bb = ((b % prime) + prime) % prime;
      if (bb == 0) {
        return false;
      }
      long long value = ((a % prime) + prime) % prime;
      value = value * power_mod(bb, static_cast<unsigned>(prime - 2), prime) % prime;

      for (std::size_t v = 0; v < P; ++v) {
        if (v == main) {
          continue;
        }
        value = value * power_mod(at[v], static_cast<unsigned>(t.mon[v]), prime) % prime;
      }
      const std::size_t k = static_cast<std::size_t>(t.mon[main]);
      out[k] = (out[k] + value) % prime;
    }
    while (out.size() > 1 && out.back() == 0) {
      out.pop_back();
    }
    return true;
  }

  // Degree of the gcd of two univariate polynomials over Z/prime, by Euclid.
  static std::size_t univariate_gcd_degree(std::vector<long long> a, std::vector<long long> b,
                                           long long prime) {
    const auto trim = [](std::vector<long long>& p) {
      while (p.size() > 1 && p.back() == 0) {
        p.pop_back();
      }
    };
    const auto vanishes = [](const std::vector<long long>& p) {
      return p.size() == 1 && p[0] == 0;
    };

    trim(a);
    trim(b);
    while (!vanishes(b)) {
      // a %= b. Each step cancels the leading term, so the degree falls and the
      // loop ends.
      const long long inverse = power_mod(b.back(), static_cast<unsigned>(prime - 2), prime);
      while (!vanishes(a) && a.size() >= b.size()) {
        const long long factor = a.back() % prime * inverse % prime;
        const std::size_t shift = a.size() - b.size();
        for (std::size_t i = 0; i < b.size(); ++i) {
          const long long updated = (a[i + shift] - factor * b[i] % prime) % prime;
          a[i + shift] = (updated + prime) % prime;
        }
        trim(a);
      }
      a.swap(b);
    }
    trim(a);
    return a.size() - 1;
  }

  // The polynomial common factor, divided out of each. The gcd is monic and
  // divides both, so each division is exact.
  static void remove_common_factor(parameter_polynomial& n, parameter_polynomial& d) {
    if (n.degree() == 0 || d.degree() == 0) {
      return;  // one side is a unit; nothing of positive degree can be shared
    }
    if (certainly_coprime(n, d)) {
      return;  // nothing to cancel, and no exact gcd computed to find that out
    }
    const parameter_polynomial g = polynomial_gcd(n, d);
    if (g.degree() == 0) {
      return;
    }
    n = exact_quotient(n, g);
    d = exact_quotient(d, g);
  }

  static parameter_polynomial exact_quotient(const parameter_polynomial& p,
                                             const parameter_polynomial& g) {
    const auto result = divide(p, std::vector<parameter_polynomial>{g});
    VARIETAS_ASSERT(result.remainder.is_zero());
    return result.quotients.front();
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
