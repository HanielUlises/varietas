#ifndef VARIETAS_CODEGEN_MODULAR_GCD_HPP
#define VARIETAS_CODEGEN_MODULAR_GCD_HPP

#include <algorithm>
#include <array>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <map>
#include <utility>
#include <vector>

#include <gmpxx.h>

#include "varietas/codegen/rational.hpp"
#include "varietas/core/config.hpp"
#include "varietas/core/monomial.hpp"
#include "varietas/core/polynomial.hpp"

namespace varietas {

// Greatest common divisor in Q[x_1, ..., x_N] by evaluation and interpolation.
//
// This is Brown's dense modular algorithm, and it exists because of a
// measurement rather than on principle. The subresultant remainder sequence in
// varietas/core/gcd.hpp costs about the fourth power of the number of terms in
// its operands, over Q and over a prime field alike, so the size of the
// coefficients was never the problem and computing the same sequence modulo
// several primes would not have been the answer. What has to go is the
// sequence itself.
//
// Brown's construction never forms a multivariate remainder. The integer
// polynomials are reduced modulo a word-sized prime; there, one variable at a
// time is specialised at points of the field until only the main variable is
// left, where the gcd is a Euclidean algorithm on machine words; and the
// multivariate gcd is rebuilt from those univariate images by Newton
// interpolation, one variable at a time on the way back up. The images modulo
// several primes are then combined by the Chinese remainder theorem, and the
// rational coefficients read back by rational reconstruction.
//
// Nothing in that chain is a proof by itself, so the answer is checked at the
// end by exact division of both operands, over Z, which by Gauss's lemma
// decides divisibility over Q. A candidate that fails the check is thrown away
// and the computation continues with more primes. The
// algorithm is therefore allowed to be wrong in the middle, where it is fast,
// and is never wrong at the end.
//
// Two facts carry the rest, and both concern what an image can look like.
// Specialising a variable, or reducing modulo a prime, maps the true gcd G onto
// a common divisor of the two images, so an image gcd is always a multiple of
// the image of G. When the leading coefficient of G survives the map, the
// image of G keeps its leading monomial, and an image gcd can then only have
// that leading monomial or a larger one. Larger means unlucky, and the image is
// discarded; the smallest leading monomial seen so far is the one believed.
// Keeping the leading coefficient alive is what the choice of points and
// primes is for.
//
// The result is monic under Order, the same normalisation polynomial_gcd uses,
// so the two are interchangeable.

// Counters for the cost measurements in doc/experiments. They are not read by
// the library.
struct modular_gcd_statistics {
  std::size_t calls = 0;
  std::size_t primes = 0;
  std::size_t trial_divisions = 0;
  std::size_t failed_trials = 0;
  std::size_t evaluation_points = 0;       // univariate images, all levels
  std::size_t interpolation_checks = 0;    // divisions modulo p inside the recursion
  std::size_t interpolation_rejections = 0;
  std::size_t largest_operand = 0;  // terms
  double seconds = 0.0;             // in modular_gcd, all of it
  double trial_seconds = 0.0;       // of which in the exact trial division
};

inline modular_gcd_statistics& modular_gcd_counters() {
  static modular_gcd_statistics counters;
  return counters;
}

namespace modular_detail {

using word = std::uint64_t;

// Arithmetic in Z/p for a prime below 2^31, so that a product of two residues
// fits in sixty-three bits.
struct prime_field {
  word p;

  word add(word a, word b) const noexcept {
    const word s = a + b;
    return s >= p ? s - p : s;
  }
  word sub(word a, word b) const noexcept { return a >= b ? a - b : a + p - b; }
  word mul(word a, word b) const noexcept { return a * b % p; }
  word neg(word a) const noexcept { return a == 0 ? 0 : p - a; }

  word pow(word base, word exponent) const noexcept {
    word result = 1;
    base %= p;
    while (exponent != 0) {
      if ((exponent & 1u) != 0u) {
        result = mul(result, base);
      }
      base = mul(base, base);
      exponent >>= 1;
    }
    return result;
  }

  // Extended Euclid rather than Fermat, a few dozen divisions against the
  // thirty-one squarings of a^(p-2), on the path every univariate division
  // takes, which made it the most frequent single call in a parametric solve.
  word inv(word a) const noexcept {
    VARIETAS_ASSERT(a % p != 0);
    std::int64_t r0 = static_cast<std::int64_t>(p);
    std::int64_t r1 = static_cast<std::int64_t>(a % p);
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
    return static_cast<word>(t0 < 0 ? t0 + static_cast<std::int64_t>(p) : t0);
  }
};

// A univariate polynomial over Z/p, coefficients in increasing degree, with no
// trailing zeros. The zero polynomial is the empty vector.
using upoly = std::vector<word>;

inline void trim(upoly& a) {
  while (!a.empty() && a.back() == 0) {
    a.pop_back();
  }
}

inline int degree(const upoly& a) { return static_cast<int>(a.size()) - 1; }

inline word evaluate(const prime_field& F, const upoly& a, word x) {
  word value = 0;
  for (auto it = a.rbegin(); it != a.rend(); ++it) {
    value = F.add(F.mul(value, x), *it);
  }
  return value;
}

inline upoly scale(const prime_field& F, upoly a, word c) {
  if (c == 0) {
    return {};
  }
  for (word& coefficient : a) {
    coefficient = F.mul(coefficient, c);
  }
  return a;
}

inline upoly add(const prime_field& F, upoly a, const upoly& b) {
  if (a.size() < b.size()) {
    a.resize(b.size(), 0);
  }
  for (std::size_t i = 0; i < b.size(); ++i) {
    a[i] = F.add(a[i], b[i]);
  }
  trim(a);
  return a;
}

inline upoly multiply(const prime_field& F, const upoly& a, const upoly& b) {
  if (a.empty() || b.empty()) {
    return {};
  }
  upoly r(a.size() + b.size() - 1, 0);
  for (std::size_t i = 0; i < a.size(); ++i) {
    if (a[i] == 0) {
      continue;
    }
    for (std::size_t j = 0; j < b.size(); ++j) {
      r[i + j] = F.add(r[i + j], F.mul(a[i], b[j]));
    }
  }
  trim(r);
  return r;
}

inline upoly monic(const prime_field& F, upoly a) {
  if (a.empty()) {
    return a;
  }
  return scale(F, std::move(a), F.inv(a.back()));
}

// Quotient and remainder. The divisor is nonzero.
inline std::pair<upoly, upoly> divide(const prime_field& F, upoly a, const upoly& b) {
  VARIETAS_ASSERT(!b.empty());
  if (a.size() < b.size()) {
    return {upoly{}, std::move(a)};
  }
  upoly q(a.size() - b.size() + 1, 0);
  const word inverse = F.inv(b.back());
  for (std::size_t shift = a.size() - b.size() + 1; shift-- > 0;) {
    const word factor = F.mul(a[shift + b.size() - 1], inverse);
    q[shift] = factor;
    if (factor != 0) {
      for (std::size_t i = 0; i < b.size(); ++i) {
        a[i + shift] = F.sub(a[i + shift], F.mul(factor, b[i]));
      }
    }
  }
  a.resize(b.size() - 1);
  trim(a);
  trim(q);
  return {std::move(q), std::move(a)};
}

// Euclid, monic. The gcd of zero and zero is zero.
inline upoly gcd(const prime_field& F, upoly a, upoly b) {
  while (!b.empty()) {
    upoly r = divide(F, std::move(a), b).second;
    a = std::move(b);
    b = std::move(r);
  }
  return monic(F, std::move(a));
}

// A multivariate polynomial over Z/p. Terms are sorted strictly decreasing in
// the lexicographic order with x_0 most significant, and carry no zero
// coefficient.
//
// That one choice is what makes the recursion cheap. At level k only the
// variables x_0, ..., x_k are in play, and the terms that agree in x_0, ...,
// x_{k-1} are contiguous and ordered by their exponent of x_k, so reading the
// polynomial as one in x_0, ..., x_{k-1} with coefficients in Z/p[x_k] is a
// single pass that cuts the term list into runs, and nothing is sorted.
template <std::size_t N>
struct mterm {
  std::array<std::uint16_t, N> e;
  word c;
};

template <std::size_t N>
using mpoly = std::vector<mterm<N>>;

template <std::size_t N>
bool lex_greater(const std::array<std::uint16_t, N>& a, const std::array<std::uint16_t, N>& b) {
  for (std::size_t i = 0; i < N; ++i) {
    if (a[i] != b[i]) {
      return a[i] > b[i];
    }
  }
  return false;
}

// Compares the exponents of x_0, ..., x_{k-1}, ignoring x_k and beyond.
template <std::size_t N>
int compare_head(const std::array<std::uint16_t, N>& a, const std::array<std::uint16_t, N>& b,
                 std::size_t k) {
  for (std::size_t i = 0; i < k; ++i) {
    if (a[i] != b[i]) {
      return a[i] > b[i] ? 1 : -1;
    }
  }
  return 0;
}

template <std::size_t N>
void sort_lex(mpoly<N>& a) {
  std::sort(a.begin(), a.end(),
            [](const mterm<N>& s, const mterm<N>& t) { return lex_greater(s.e, t.e); });
}

// The polynomial at level k read with coefficients in Z/p[x_k]: one entry per
// monomial in x_0, ..., x_{k-1}, in decreasing lexicographic order, each with
// its coefficient as a dense univariate polynomial in x_k.
template <std::size_t N>
struct run {
  std::array<std::uint16_t, N> head;  // x_k and beyond set to zero
  upoly coefficient;
};

template <std::size_t N>
std::vector<run<N>> split(const mpoly<N>& a, std::size_t k) {
  std::vector<run<N>> runs;
  for (const auto& t : a) {
    if (runs.empty() || compare_head(runs.back().head, t.e, k) != 0) {
      run<N> r;
      r.head = t.e;
      for (std::size_t i = k; i < N; ++i) {
        r.head[i] = 0;
      }
      runs.push_back(std::move(r));
    }
    upoly& u = runs.back().coefficient;
    const std::size_t d = t.e[k];
    if (u.size() <= d) {
      u.resize(d + 1, 0);
    }
    u[d] = t.c;
  }
  return runs;
}

template <std::size_t N>
mpoly<N> join(const std::vector<run<N>>& runs, std::size_t k) {
  mpoly<N> a;
  for (const auto& r : runs) {
    for (std::size_t d = r.coefficient.size(); d-- > 0;) {
      if (r.coefficient[d] == 0) {
        continue;
      }
      mterm<N> t{r.head, r.coefficient[d]};
      t.e[k] = static_cast<std::uint16_t>(d);
      a.push_back(t);
    }
  }
  return a;
}

// x_k set to alpha. Contiguity again makes every run collapse to one term.
template <std::size_t N>
mpoly<N> specialise(const prime_field& F, const std::vector<run<N>>& runs, word alpha) {
  mpoly<N> a;
  for (const auto& r : runs) {
    const word value = evaluate(F, r.coefficient, alpha);
    if (value != 0) {
      a.push_back(mterm<N>{r.head, value});
    }
  }
  return a;
}

template <std::size_t N>
mpoly<N> multiply_scalar(const prime_field& F, mpoly<N> a, word c) {
  if (c == 0) {
    return {};
  }
  for (auto& t : a) {
    t.c = F.mul(t.c, c);
  }
  return a;
}

// a - c * m * b, as one merge, both sorted and the shift by a monomial keeping
// b sorted.
template <std::size_t N>
mpoly<N> subtract_multiple(const prime_field& F, const mpoly<N>& a,
                           const std::array<std::uint16_t, N>& m, word c, const mpoly<N>& b) {
  mpoly<N> r;
  r.reserve(a.size() + b.size());
  std::size_t i = 0;
  std::size_t j = 0;
  const auto shifted = [&](std::size_t idx) {
    auto e = b[idx].e;
    for (std::size_t v = 0; v < N; ++v) {
      e[v] = static_cast<std::uint16_t>(e[v] + m[v]);
    }
    return e;
  };
  while (i < a.size() || j < b.size()) {
    if (j == b.size()) {
      r.push_back(a[i++]);
      continue;
    }
    const auto e = shifted(j);
    if (i == a.size() || lex_greater(e, a[i].e)) {
      r.push_back(mterm<N>{e, F.neg(F.mul(c, b[j].c))});
      ++j;
    } else if (lex_greater(a[i].e, e)) {
      r.push_back(a[i++]);
    } else {
      const word v = F.sub(a[i].c, F.mul(c, b[j].c));
      if (v != 0) {
        r.push_back(mterm<N>{e, v});
      }
      ++i;
      ++j;
    }
  }
  return r;
}

// True when b divides a in Z/p[x_0, ..., x_{N-1}].
template <std::size_t N>
bool divides(const prime_field& F, const mpoly<N>& b, mpoly<N> a) {
  VARIETAS_ASSERT(!b.empty());
  const word inverse = F.inv(b.front().c);
  while (!a.empty()) {
    std::array<std::uint16_t, N> m{};
    for (std::size_t v = 0; v < N; ++v) {
      if (a.front().e[v] < b.front().e[v]) {
        return false;
      }
      m[v] = static_cast<std::uint16_t>(a.front().e[v] - b.front().e[v]);
    }
    a = subtract_multiple(F, a, m, F.mul(a.front().c, inverse), b);
  }
  return true;
}

// The gcd in Z/p[x_0, ..., x_k], with every variable past x_k absent from both
// operands. Normalised so that its lexicographic leading coefficient is one.
//
// `seed` shifts the evaluation points, so that a call that has to be repeated
// with a different prime also meets different points.
template <std::size_t N>
mpoly<N> brown(const prime_field& F, const mpoly<N>& a, const mpoly<N>& b, std::size_t k,
               word seed) {
  if (a.empty() || b.empty()) {
    mpoly<N> r = a.empty() ? b : a;
    return r.empty() ? r : multiply_scalar(F, r, F.inv(r.front().c));
  }

  if (k == 0) {
    // With one variable left, Euclid on dense vectors.
    upoly ua = split(a, 0).front().coefficient;
    upoly ub = split(b, 0).front().coefficient;
    std::vector<run<N>> g(1);
    g[0].head = {};
    g[0].coefficient = gcd(F, std::move(ua), std::move(ub));
    return join(g, 0);
  }

  auto ra = split(a, k);
  auto rb = split(b, k);

  // Contents in Z/p[x_k], and the primitive parts.
  upoly content_a;
  for (const auto& r : ra) {
    content_a = gcd(F, std::move(content_a), r.coefficient);
  }
  upoly content_b;
  for (const auto& r : rb) {
    content_b = gcd(F, std::move(content_b), r.coefficient);
  }
  const upoly content = gcd(F, content_a, content_b);
  for (auto& r : ra) {
    r.coefficient = divide(F, r.coefficient, content_a).first;
  }
  for (auto& r : rb) {
    r.coefficient = divide(F, r.coefficient, content_b).first;
  }

  // The leading coefficients, in x_0, ..., x_{k-1}, are polynomials in x_k.
  // Their gcd is a multiple of the leading coefficient of the answer, and is
  // what every interpolated image is scaled to lead with, so that the images
  // agree on a normalisation before they are combined.
  const upoly gamma = gcd(F, ra.front().coefficient, rb.front().coefficient);

  int degree_a = 0;
  for (const auto& r : ra) {
    degree_a = std::max(degree_a, degree(r.coefficient));
  }
  int degree_b = 0;
  for (const auto& r : rb) {
    degree_b = std::max(degree_b, degree(r.coefficient));
  }
  // The answer times gamma / lc has at most this degree in x_k, so this many
  // points plus one determine it.
  const int bound = degree(gamma) + std::min(degree_a, degree_b);

  // The candidate the interpolant stands for, which is its primitive part in x_0, ...,
  // x_{k-1}, the content of the answer put back, and the lexicographic leading
  // coefficient made one. Accepted only if it divides both operands.
  const auto accept = [&](const std::vector<run<N>>& h, mpoly<N>& out) {
    upoly content_h;
    for (const auto& r : h) {
      content_h = gcd(F, std::move(content_h), r.coefficient);
    }
    std::vector<run<N>> primitive = h;
    for (auto& r : primitive) {
      r.coefficient = divide(F, r.coefficient, content_h).first;
    }
    mpoly<N> candidate = join(primitive, k);
    candidate = multiply_scalar(F, std::move(candidate), F.inv(candidate.front().c));
    ++modular_gcd_counters().interpolation_checks;
    if (!divides(F, candidate, join(ra, k)) || !divides(F, candidate, join(rb, k))) {
      ++modular_gcd_counters().interpolation_rejections;
      return false;
    }
    for (auto& r : primitive) {
      r.coefficient = multiply(F, r.coefficient, content);
    }
    out = join(primitive, k);
    out = multiply_scalar(F, std::move(out), F.inv(out.front().c));
    return true;
  };

  std::vector<run<N>> h;  // the interpolant, one univariate per monomial
  upoly q{1};             // prod (x_k - alpha) over the points used
  std::array<std::uint16_t, N> leading{};
  bool have = false;
  int points = 0;

  for (word alpha = 1 + seed % (F.p - 1);; alpha = alpha % (F.p - 1) + 1) {
    const word g_alpha = evaluate(F, gamma, alpha);
    if (g_alpha == 0) {
      continue;  // the leading coefficient of the answer may vanish here
    }
    ++modular_gcd_counters().evaluation_points;
    mpoly<N> image = brown(F, specialise(F, ra, alpha), specialise(F, rb, alpha), k - 1, seed);
    image = multiply_scalar(F, std::move(image), g_alpha);  // image was monic

    const auto& lm = image.front().e;
    if (have && lex_greater(lm, leading)) {
      continue;  // unlucky point: this image carries a spurious factor
    }
    bool stable = false;
    if (!have || lex_greater(leading, lm)) {
      // First image, or every previous one was unlucky. Start again.
      h = split(image, k);
      q = upoly{F.neg(alpha), 1};
      leading = lm;
      have = true;
      points = 1;
    } else {
      // Newton interpolation, h + (image - h(alpha)) q / q(alpha), monomial by monomial,
      // over the union of the two supports.
      const word q_inverse = F.inv(evaluate(F, q, alpha));
      stable = true;
      std::vector<run<N>> updated;
      std::size_t i = 0;
      std::size_t j = 0;
      while (i < h.size() || j < image.size()) {
        run<N> r;
        word target = 0;
        const int cmp = i == h.size()       ? -1
                        : j == image.size() ? 1
                                            : compare_head(h[i].head, image[j].e, k);
        if (cmp > 0) {
          r = h[i++];
        } else if (cmp < 0) {
          r.head = image[j].e;
          target = image[j++].c;
        } else {
          r = h[i++];
          target = image[j++].c;
        }
        const word correction = F.sub(target, evaluate(F, r.coefficient, alpha));
        if (correction != 0) {
          stable = false;
          r.coefficient = add(F, r.coefficient, scale(F, q, F.mul(correction, q_inverse)));
        }
        if (!r.coefficient.empty()) {
          updated.push_back(std::move(r));
        }
      }
      h = std::move(updated);
      q = multiply(F, q, upoly{F.neg(alpha), 1});
      ++points;
    }

    // Two reasons to try the candidate. Past the degree bound the interpolant
    // is determined, so a failure there means some image was wrong, which the
    // division catches and more points repair. Before the bound, a point that
    // changed nothing suggests the interpolant has already stopped moving,
    // which is the common case, since the bound is for the worst case and sparse
    // answers of low degree settle long before it.
    if (points > bound || stable) {
      mpoly<N> result;
      if (accept(h, result)) {
        return result;
      }
    }
  }
}

// Primes below 2^31, from the top down. Trial division is enough at this size
// and the list is only walked a few entries deep.
inline word nth_prime_below_2_31(std::size_t n) {
  static std::vector<word> found;
  word candidate = found.empty() ? 2147483647u : found.back() - 2;
  while (found.size() <= n) {
    bool prime = true;
    for (word d = 3; d * d <= candidate; d += 2) {
      if (candidate % d == 0) {
        prime = false;
        break;
      }
    }
    if (prime) {
      found.push_back(candidate);
    }
    candidate -= 2;
  }
  return found[n];
}

inline word reduce(const mpz_class& z, word p) {
  return static_cast<word>(mpz_fdiv_ui(z.get_mpz_t(), static_cast<unsigned long>(p)));
}

// n/d with |n|, d below sqrt(m/2) and n/d = u mod m, if there is one.
inline bool rational_reconstruction(const mpz_class& u, const mpz_class& m, rational& out) {
  mpz_class bound;
  mpz_class half = m / 2;
  mpz_sqrt(bound.get_mpz_t(), half.get_mpz_t());

  mpz_class r0 = m;
  mpz_class r1 = u;
  mpz_class t0 = 0;
  mpz_class t1 = 1;
  while (r1 > bound) {
    const mpz_class quotient = r0 / r1;
    mpz_class r2 = r0 - quotient * r1;
    mpz_class t2 = t0 - quotient * t1;
    r0 = std::move(r1);
    r1 = std::move(r2);
    t0 = std::move(t1);
    t1 = std::move(t2);
  }
  if (t1 == 0 || abs(t1) > bound) {
    return false;
  }
  mpz_class g;
  mpz_gcd(g.get_mpz_t(), r1.get_mpz_t(), t1.get_mpz_t());
  if (g != 1) {
    return false;
  }
  out = rational(r1, t1);
  out.canonicalize();
  return true;
}

// p written as content * (an integer polynomial with no common factor), the
// content rational. The terms keep the order of p.
template <std::size_t N>
struct integer_form {
  rational content;
  std::vector<std::pair<monomial<N>, mpz_class>> terms;
};

template <std::size_t N, class Order>
integer_form<N> primitive_integer(const polynomial<rational, N, Order>& p) {
  mpz_class denominator_lcm = 1;
  mpz_class numerator_gcd = 0;
  for (const auto& t : p.terms()) {
    mpz_lcm(denominator_lcm.get_mpz_t(), denominator_lcm.get_mpz_t(),
            t.coeff.get_den().get_mpz_t());
    mpz_gcd(numerator_gcd.get_mpz_t(), numerator_gcd.get_mpz_t(),
            t.coeff.get_num().get_mpz_t());
  }
  integer_form<N> out;
  out.content = rational(numerator_gcd, denominator_lcm);
  out.content.canonicalize();
  out.terms.reserve(p.size());
  for (const auto& t : p.terms()) {
    mpz_class z = t.coeff.get_num() * (denominator_lcm / t.coeff.get_den()) / numerator_gcd;
    out.terms.emplace_back(t.mon, std::move(z));
  }
  return out;
}

// Exact division over Z, or false as soon as it is clear there is none.
//
// This is the certificate the whole algorithm rests on, and it used to be the
// general division routine over Q, which was most of the cost of a gcd because it
// rewrites the entire remainder at every step, rational coefficients and all.
// Over Z, with the remainder held in an ordered map and updated in place, a
// step costs the size of the divisor rather than the size of the remainder,
// and the arithmetic is on integers that need no canonicalising.
//
// Gauss's lemma is what makes Z the right ring. The divisor is primitive, so
// if it divides a primitive integer polynomial over Q it divides it over Z,
// and a quotient coefficient that fails to be an integer is already a proof
// that it does not divide.
template <std::size_t N, class Order>
bool divide_integer(const std::vector<std::pair<monomial<N>, mpz_class>>& a,
                    const std::vector<std::pair<monomial<N>, mpz_class>>& g,
                    std::vector<std::pair<monomial<N>, mpz_class>>& quotient) {
  VARIETAS_ASSERT(!g.empty());
  const auto descending = [](const monomial<N>& s, const monomial<N>& t) {
    return Order::compare(s, t) > 0;
  };
  std::map<monomial<N>, mpz_class, decltype(descending)> remainder(descending);
  for (const auto& [m, c] : a) {
    remainder.emplace(m, c);
  }
  const monomial<N>& lead = g.front().first;
  const mpz_class& lead_coefficient = g.front().second;

  quotient.clear();
  mpz_class q;
  mpz_class r;
  while (!remainder.empty()) {
    const auto head = remainder.begin();
    if (!monomial<N>::divides(lead, head->first)) {
      return false;
    }
    mpz_tdiv_qr(q.get_mpz_t(), r.get_mpz_t(), head->second.get_mpz_t(),
                lead_coefficient.get_mpz_t());
    if (r != 0) {
      return false;
    }
    const monomial<N> shift = monomial<N>::divide(head->first, lead);
    remainder.erase(head);
    for (std::size_t i = 1; i < g.size(); ++i) {
      auto [it, inserted] = remainder.try_emplace(shift * g[i].first);
      mpz_submul(it->second.get_mpz_t(), q.get_mpz_t(), g[i].second.get_mpz_t());
      if (it->second == 0) {
        remainder.erase(it);
      }
    }
    quotient.emplace_back(shift, q);
  }
  return true;
}

template <std::size_t N>
mpoly<N> image_mod(const std::vector<std::pair<monomial<N>, mpz_class>>& a, word p) {
  mpoly<N> r;
  r.reserve(a.size());
  for (const auto& [mon, z] : a) {
    const word c = reduce(z, p);
    if (c != 0) {
      r.push_back(mterm<N>{mon.exponents(), c});
    }
  }
  sort_lex(r);
  return r;
}

}  // namespace modular_detail


// The gcd together with the cofactors a / gcd and b / gcd.
//
// The exact trial division that certifies the gcd computes the cofactors as its
// quotients, and a caller cancelling a common factor, which is every caller in
// this library, wants exactly those. Returning them saves doing the same
// division twice, and in a parametric solve that division is the largest part
// of the cost of the gcd.
template <class Poly>
struct gcd_and_cofactors {
  Poly gcd;
  Poly a_over_gcd;
  Poly b_over_gcd;
};

template <std::size_t N, class Order>
gcd_and_cofactors<polynomial<rational, N, Order>> modular_gcd_with_cofactors(
    const polynomial<rational, N, Order>& a, const polynomial<rational, N, Order>& b) {
  using poly = polynomial<rational, N, Order>;
  using result = gcd_and_cofactors<poly>;
  using namespace modular_detail;

  auto& counters = modular_gcd_counters();
  ++counters.calls;
  counters.largest_operand = std::max({counters.largest_operand, a.size(), b.size()});
  struct timer {
    double& total;
    std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
    ~timer() {
      total += std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    }
  } timed{counters.seconds};

  const poly one = poly::constant(rational(1));
  if (a.is_zero()) {
    if (b.is_zero()) {
      return result{poly(), poly(), poly()};
    }
    return result{b.monic(), poly(), poly::constant(b.leading_coefficient())};
  }
  if (b.is_zero()) {
    return result{a.monic(), poly::constant(a.leading_coefficient()), poly()};
  }
  if (a.degree() == 0 || b.degree() == 0) {
    return result{one, a, b};
  }

  const auto A = primitive_integer(a);
  const auto B = primitive_integer(b);
  // The leading coefficients under Order. A prime that divides neither keeps
  // the leading monomial of the true gcd in its image by Gauss's lemma, since the
  // primitive gcd divides the primitive operands over Z, so its leading
  // coefficient divides theirs.
  const mpz_class& lead_a = A.terms.front().second;
  const mpz_class& lead_b = B.terms.front().second;

  // The images combined so far, as coefficient residues modulo `modulus` keyed
  // by monomial, for the images that share the smallest leading monomial seen.
  std::map<std::array<std::uint16_t, N>, mpz_class> residues;
  mpz_class modulus = 1;
  monomial<N> leading;
  bool have = false;
  poly previous;

  for (std::size_t index = 0;; ++index) {
    const word p = nth_prime_below_2_31(index);
    if (reduce(lead_a, p) == 0 || reduce(lead_b, p) == 0) {
      continue;
    }
    ++counters.primes;
    const prime_field F{p};

    mpoly<N> image = brown(F, image_mod(A.terms, p), image_mod(B.terms, p), N - 1, index);

    // A constant image is conclusive at once. Every image is a multiple of the
    // image of the gcd, so the gcd itself is constant, and that is the answer
    // to most of the calls normalisation makes.
    if (image.size() == 1 && image.front().e == std::array<std::uint16_t, N>{}) {
      return result{one, a, b};
    }

    // Monic under Order rather than lex, which is the normalisation of the
    // answer and so the one under which the images agree coefficient by
    // coefficient.
    std::size_t top = 0;
    for (std::size_t i = 1; i < image.size(); ++i) {
      if (Order::compare(monomial<N>(image[i].e), monomial<N>(image[top].e)) > 0) {
        top = i;
      }
    }
    const monomial<N> lm(image[top].e);
    const word scale = F.inv(image[top].c);

    if (have && Order::compare(lm, leading) > 0) {
      continue;  // unlucky prime
    }
    if (!have || Order::compare(lm, leading) < 0) {
      residues.clear();
      modulus = 1;
      leading = lm;
      have = true;
    }

    // Chinese remaindering, monomial by monomial. A monomial absent from one
    // side has residue zero there.
    std::map<std::array<std::uint16_t, N>, word> current;
    for (const auto& t : image) {
      current[t.e] = F.mul(t.c, scale);
    }
    for (const auto& [e, c] : current) {
      residues.emplace(e, mpz_class(0));
    }
    const word m_mod_p = reduce(modulus, p);
    const word m_inverse = F.inv(m_mod_p);
    for (auto& [e, r] : residues) {
      const auto it = current.find(e);
      const word target = it == current.end() ? 0 : it->second;
      const word delta = F.mul(F.sub(target, reduce(r, p)), m_inverse);
      r += modulus * static_cast<unsigned long>(delta);
    }
    modulus *= static_cast<unsigned long>(p);

    // Rational reconstruction, and a trial division once two consecutive
    // reconstructions agree, because the check is exact and is the expensive part, so
    // it waits for the coefficients to stop moving.
    std::vector<typename poly::term> terms;
    bool reconstructed = true;
    for (const auto& [e, r] : residues) {
      rational c;
      if (!rational_reconstruction(r, modulus, c)) {
        reconstructed = false;
        break;
      }
      if (sgn(c) != 0) {
        terms.push_back({monomial<N>(e), c});
      }
    }
    if (!reconstructed) {
      continue;
    }
    poly candidate(std::move(terms));
    if (candidate != previous) {
      previous = std::move(candidate);
      continue;
    }

    // The candidate is monic already, since every image was scaled to lead with one
    // and one reconstructs to one. It is certified by dividing both primitive
    // operands by its primitive integer form, and the quotients are then
    // scaled back into the cofactors over Q: with a = c_a A and the candidate
    // g = G / lc(G), a / g = c_a lc(G) (A / G).
    ++counters.trial_divisions;
    timer trial{counters.trial_seconds};
    const auto G = primitive_integer(candidate);
    std::vector<std::pair<monomial<N>, mpz_class>> over_a;
    std::vector<std::pair<monomial<N>, mpz_class>> over_b;
    if (divide_integer<N, Order>(A.terms, G.terms, over_a) &&
        divide_integer<N, Order>(B.terms, G.terms, over_b)) {
      const auto cofactor = [&](const integer_form<N>& operand,
                                const std::vector<std::pair<monomial<N>, mpz_class>>& q) {
        const rational scale = operand.content * rational(G.terms.front().second);
        std::vector<typename poly::term> terms;
        terms.reserve(q.size());
        for (const auto& [m, c] : q) {
          terms.push_back({m, scale * rational(c)});
        }
        return poly(std::move(terms));
      };
      return result{std::move(candidate), cofactor(A, over_a), cofactor(B, over_b)};
    }
    ++counters.failed_trials;
  }
}

template <std::size_t N, class Order>
polynomial<rational, N, Order> modular_gcd(const polynomial<rational, N, Order>& a,
                                           const polynomial<rational, N, Order>& b) {
  return modular_gcd_with_cofactors(a, b).gcd;
}

}  // namespace varietas

#endif
