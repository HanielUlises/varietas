#ifndef VARIETAS_CODEGEN_FACTOR_HPP
#define VARIETAS_CODEGEN_FACTOR_HPP

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <map>
#include <numeric>
#include <random>
#include <utility>
#include <vector>

#include <gmpxx.h>

#include "varietas/codegen/modular_gcd.hpp"
#include "varietas/codegen/rational.hpp"
#include "varietas/core/config.hpp"
#include "varietas/core/gcd.hpp"
#include "varietas/core/monomial.hpp"
#include "varietas/core/polynomial.hpp"

namespace varietas {

// Factorisation over Q, the first part of it: squarefree decomposition in any
// number of variables, and complete factorisation into irreducibles in one.
//
// What split_along lacks is the ability to choose its own h, and the
// irreducible factors of a polynomial are exactly the hypersurfaces a
// decomposition would split along. Multivariate factorisation is reached
// through the univariate one, by specialising all but one variable, factoring
// what is left and lifting the factors back, so the univariate algorithm is
// where it has to start, and this is it.
//
// The squarefree decomposition is Yun's algorithm, which needs nothing but
// derivatives and gcds and so works in every number of variables at once: in
// characteristic zero a squarefree factor and its derivative share nothing,
// and a factor of multiplicity m divides the derivative m - 1 times. The gcds
// are modular_gcd's, and the cofactors it returns are the exact quotients the
// algorithm needs next.
//
// The univariate factorisation is Zassenhaus's. A squarefree integer
// polynomial is reduced modulo a prime that keeps it squarefree and of the
// same degree, factored there by Cantor and Zassenhaus, the factors are lifted
// by Hensel's lemma to a power of the prime past a bound on the coefficients
// of any true factor, and the true factors are recovered as products of the
// lifted ones. The last step is the expensive one and the only one that can
// be: a polynomial irreducible over Q can split into many factors modulo every
// prime, the polynomials of Swinnerton-Dyer being the standard example, and
// then every subset of the modular factors is a candidate. The lattice
// reduction of van Hoeij removes that exponential, and it is not here.
//
// The recombination is exact rather than tested. A candidate g and its
// cofactor h, both read symmetrically modulo p^l, are accepted when the
// product of their 1-norms is at most the bound B, and since p^l > 2B that
// inequality forces g h to equal lc(f) f over Z, not merely modulo p^l. Nothing
// in the answer is believed on the strength of an image. The construction
// throughout is the one in chapters 14 and 15 of von zur Gathen and Gerhard,
// Modern Computer Algebra.
//
// Factors are monic under Order, the normalisation the gcds use, and the unit
// carries the leading coefficient.

// Counters for the cost measurements in doc/experiments. They are not read by
// the library.
struct factor_statistics {
  std::size_t calls = 0;            // squarefree integer polynomials factored
  std::size_t primes = 0;           // primes at which a factorisation was counted
  std::size_t modular_factors = 0;  // at the prime chosen, all calls
  std::size_t true_factors = 0;
  std::size_t subsets = 0;          // candidates through the norm test
  std::size_t precision_bits = 0;   // the largest p^l lifted to
};

inline factor_statistics& factor_counters() {
  static factor_statistics counters;
  return counters;
}

// f = unit * product of base^multiplicity.
template <class Poly>
struct factor_power {
  Poly base;
  unsigned multiplicity;
};

template <class Poly>
struct factorisation {
  rational unit;
  std::vector<factor_power<Poly>> factors;
};

namespace factor_detail {

namespace md = modular_detail;
using md::word;

// Univariate arithmetic over Z/p beyond what the gcd needed.

inline md::upoly subtract(const md::prime_field& F, const md::upoly& a, const md::upoly& b) {
  return md::add(F, a, md::scale(F, b, F.p - 1));
}

inline md::upoly remainder(const md::prime_field& F, md::upoly a, const md::upoly& b) {
  return md::divide(F, std::move(a), b).second;
}

inline md::upoly derivative(const md::prime_field& F, const md::upoly& a) {
  md::upoly r;
  for (std::size_t i = 1; i < a.size(); ++i) {
    r.push_back(F.mul(a[i], static_cast<word>(i) % F.p));
  }
  md::trim(r);
  return r;
}

// base^e modulo m, by squaring along the bits of e.
inline md::upoly power_mod(const md::prime_field& F, const md::upoly& base, const mpz_class& e,
                           const md::upoly& m) {
  const md::upoly b = remainder(F, base, m);
  md::upoly result{1};
  for (std::size_t i = mpz_sizeinbase(e.get_mpz_t(), 2); i-- > 0;) {
    result = remainder(F, md::multiply(F, result, result), m);
    if (mpz_tstbit(e.get_mpz_t(), i) != 0) {
      result = remainder(F, md::multiply(F, result, b), m);
    }
  }
  return result;
}

// s a + t b = 1, with deg s < deg b and deg t < deg a, for a and b coprime.
inline std::pair<md::upoly, md::upoly> bezout(const md::prime_field& F, const md::upoly& a,
                                              const md::upoly& b) {
  md::upoly r0 = a;
  md::upoly r1 = b;
  md::upoly s0{1};
  md::upoly s1;
  md::upoly t0;
  md::upoly t1{1};
  while (!r1.empty()) {
    auto [q, r] = md::divide(F, r0, r1);
    r0 = std::move(r1);
    r1 = std::move(r);
    md::upoly s2 = subtract(F, s0, md::multiply(F, q, s1));
    md::upoly t2 = subtract(F, t0, md::multiply(F, q, t1));
    s0 = std::move(s1);
    s1 = std::move(s2);
    t0 = std::move(t1);
    t1 = std::move(t2);
  }
  VARIETAS_ASSERT(md::degree(r0) == 0);
  const word inverse = F.inv(r0.front());
  s0 = md::scale(F, std::move(s0), inverse);
  t0 = md::scale(F, std::move(t0), inverse);
  // Bring s below deg b; t absorbs the difference and drops below deg a with it.
  auto [q, s] = md::divide(F, std::move(s0), b);
  return {std::move(s), md::add(F, std::move(t0), md::multiply(F, q, a))};
}

// The distinct-degree factorisation of a monic squarefree f: for each d, the
// product of the irreducible factors of degree d, which are exactly the
// factors of x^(p^d) - x not already removed.
inline std::vector<std::pair<md::upoly, int>> distinct_degree(const md::prime_field& F,
                                                              md::upoly f) {
  std::vector<std::pair<md::upoly, int>> out;
  const md::upoly x{0, 1};
  const mpz_class p = static_cast<unsigned long>(F.p);
  md::upoly h = x;
  int d = 0;
  while (md::degree(f) >= 2 * (d + 1)) {
    ++d;
    h = power_mod(F, h, p, f);
    md::upoly g = md::gcd(F, subtract(F, h, x), f);
    if (md::degree(g) > 0) {
      f = md::divide(F, std::move(f), g).first;
      h = remainder(F, std::move(h), f);
      out.emplace_back(std::move(g), d);
    }
  }
  // Nothing of degree below half of what is left divides it, so it is irreducible.
  if (md::degree(f) > 0) {
    const int degree = md::degree(f);
    out.emplace_back(std::move(f), degree);
  }
  return out;
}

// Splits g, a product of irreducibles all of degree d, by Cantor and
// Zassenhaus: for a random a, a^((p^d - 1)/2) is 1 or -1 modulo each factor,
// independently and with equal odds, so its difference with 1 shares a proper
// factor with g at least half the time. The prime is odd.
inline void equal_degree(const md::prime_field& F, const md::upoly& g, int d, std::mt19937_64& rng,
                         std::vector<md::upoly>& out) {
  if (md::degree(g) == d) {
    out.push_back(g);
    return;
  }
  mpz_class e;
  mpz_ui_pow_ui(e.get_mpz_t(), static_cast<unsigned long>(F.p), static_cast<unsigned long>(d));
  e = (e - 1) / 2;
  std::uniform_int_distribution<word> coefficient(0, F.p - 1);
  while (true) {
    md::upoly a(static_cast<std::size_t>(md::degree(g)));
    for (word& c : a) {
      c = coefficient(rng);
    }
    md::trim(a);
    if (md::degree(a) < 1) {
      continue;
    }
    const md::upoly c = md::gcd(F, subtract(F, power_mod(F, a, e, g), md::upoly{1}), g);
    if (md::degree(c) > 0 && md::degree(c) < md::degree(g)) {
      equal_degree(F, c, d, rng, out);
      equal_degree(F, md::divide(F, g, c).first, d, rng, out);
      return;
    }
  }
}

// Integer polynomials, coefficients in increasing degree, no trailing zeros,
// used both over Z and modulo p^l with representatives in [0, p^l).
using zpoly = std::vector<mpz_class>;

inline void trim(zpoly& a) {
  while (!a.empty() && a.back() == 0) {
    a.pop_back();
  }
}

inline int degree(const zpoly& a) { return static_cast<int>(a.size()) - 1; }

inline zpoly reduce(zpoly a, const mpz_class& m) {
  for (mpz_class& c : a) {
    mpz_fdiv_r(c.get_mpz_t(), c.get_mpz_t(), m.get_mpz_t());
  }
  trim(a);
  return a;
}

inline zpoly add(zpoly a, const zpoly& b, const mpz_class& m) {
  if (a.size() < b.size()) {
    a.resize(b.size(), 0);
  }
  for (std::size_t i = 0; i < b.size(); ++i) {
    a[i] += b[i];
  }
  return reduce(std::move(a), m);
}

inline zpoly subtract(zpoly a, const zpoly& b, const mpz_class& m) {
  if (a.size() < b.size()) {
    a.resize(b.size(), 0);
  }
  for (std::size_t i = 0; i < b.size(); ++i) {
    a[i] -= b[i];
  }
  return reduce(std::move(a), m);
}

inline zpoly multiply(const zpoly& a, const zpoly& b, const mpz_class& m) {
  if (a.empty() || b.empty()) {
    return {};
  }
  zpoly r(a.size() + b.size() - 1, 0);
  for (std::size_t i = 0; i < a.size(); ++i) {
    for (std::size_t j = 0; j < b.size(); ++j) {
      mpz_addmul(r[i + j].get_mpz_t(), a[i].get_mpz_t(), b[j].get_mpz_t());
    }
  }
  return reduce(std::move(r), m);
}

// Quotient and remainder modulo m by a monic divisor.
inline std::pair<zpoly, zpoly> divide_monic(zpoly a, const zpoly& b, const mpz_class& m) {
  VARIETAS_ASSERT(!b.empty() && b.back() == 1);
  if (a.size() < b.size()) {
    return {zpoly{}, std::move(a)};
  }
  zpoly q(a.size() - b.size() + 1, 0);
  mpz_class factor;
  for (std::size_t shift = q.size(); shift-- > 0;) {
    mpz_fdiv_r(factor.get_mpz_t(), a[shift + b.size() - 1].get_mpz_t(), m.get_mpz_t());
    q[shift] = factor;
    for (std::size_t i = 0; i < b.size(); ++i) {
      mpz_submul(a[i + shift].get_mpz_t(), factor.get_mpz_t(), b[i].get_mpz_t());
    }
  }
  a.resize(b.size() - 1);
  return {reduce(std::move(q), m), reduce(std::move(a), m)};
}

inline zpoly lift(const md::upoly& a) {
  zpoly r;
  r.reserve(a.size());
  for (word c : a) {
    r.emplace_back(static_cast<unsigned long>(c));
  }
  return r;
}

inline md::upoly image(const zpoly& a, word p) {
  md::upoly r;
  r.reserve(a.size());
  for (const mpz_class& c : a) {
    r.push_back(md::reduce(c, p));
  }
  md::trim(r);
  return r;
}

// Representatives in (-m/2, m/2].
inline zpoly symmetric(zpoly a, const mpz_class& m) {
  const mpz_class half = m / 2;
  for (mpz_class& c : a) {
    if (c > half) {
      c -= m;
    }
  }
  return a;
}

inline zpoly primitive(zpoly a) {
  mpz_class g = 0;
  for (const mpz_class& c : a) {
    mpz_gcd(g.get_mpz_t(), g.get_mpz_t(), c.get_mpz_t());
  }
  VARIETAS_ASSERT(g != 0);
  const bool negate = a.back() < 0;
  for (mpz_class& c : a) {
    mpz_divexact(c.get_mpz_t(), c.get_mpz_t(), g.get_mpz_t());
    if (negate) {
      c = -c;
    }
  }
  return a;
}

inline mpz_class norm_1(const zpoly& a) {
  mpz_class n = 0;
  for (const mpz_class& c : a) {
    n += abs(c);
  }
  return n;
}

// One quadratic Hensel step: from f = g h and s g + t h = 1 modulo m, with h monic, to the same modulo M
// for any M dividing m^2. Every input is read as an integer lift of its class
// modulo m, which is all the step asks of it.
struct hensel_pair {
  zpoly g;
  zpoly h;
  zpoly s;
  zpoly t;
};

inline void hensel_step(const zpoly& f, hensel_pair& x, const mpz_class& M) {
  const zpoly e = subtract(f, multiply(x.g, x.h, M), M);
  auto [q, r] = divide_monic(multiply(x.s, e, M), x.h, M);
  zpoly g = add(x.g, add(multiply(x.t, e, M), multiply(q, x.g, M), M), M);
  zpoly h = add(x.h, r, M);

  const zpoly b = subtract(add(multiply(x.s, g, M), multiply(x.t, h, M), M), zpoly{1}, M);
  auto [c, d] = divide_monic(multiply(x.s, b, M), h, M);
  x.s = subtract(x.s, d, M);
  x.t = subtract(x.t, add(multiply(x.t, b, M), multiply(c, g, M), M), M);
  x.g = std::move(g);
  x.h = std::move(h);
}

// The factors [lo, hi) of f modulo p, monic and pairwise coprime there, with
// f = lc(f) times their product modulo p, lifted to monic factors modulo P, a
// power of p. The factors are split into two halves, the two-factor
// factorisation is lifted, and each half is lifted again inside its own
// product, which is a factor tree built as it is walked.
inline void hensel_lift(const zpoly& f, const std::vector<md::upoly>& factors, std::size_t lo,
                        std::size_t hi, const md::prime_field& F, const mpz_class& P,
                        std::vector<zpoly>& out) {
  if (hi - lo == 1) {
    mpz_class inverse;
    const int invertible = mpz_invert(inverse.get_mpz_t(), f.back().get_mpz_t(), P.get_mpz_t());
    VARIETAS_ASSERT(invertible != 0);
    (void)invertible;
    out.push_back(multiply(f, zpoly{inverse}, P));
    return;
  }
  const std::size_t mid = lo + (hi - lo) / 2;
  md::upoly g0{md::reduce(f.back(), F.p)};
  for (std::size_t i = lo; i < mid; ++i) {
    g0 = md::multiply(F, g0, factors[i]);
  }
  md::upoly h0{1};
  for (std::size_t i = mid; i < hi; ++i) {
    h0 = md::multiply(F, h0, factors[i]);
  }
  auto [s0, t0] = bezout(F, g0, h0);

  hensel_pair x{lift(g0), lift(h0), lift(s0), lift(t0)};
  mpz_class m = static_cast<unsigned long>(F.p);
  while (m < P) {
    const mpz_class M = std::min<mpz_class>(m * m, P);
    hensel_step(f, x, M);
    m = M;
  }
  hensel_lift(x.g, factors, lo, mid, F, P, out);
  hensel_lift(x.h, factors, mid, hi, F, P, out);
}

// Advances a combination of s indices from {0, ..., n - 1}, lexicographically.
inline bool next_combination(std::vector<std::size_t>& index, std::size_t n) {
  const std::size_t s = index.size();
  for (std::size_t i = s; i-- > 0;) {
    if (index[i] < n - s + i) {
      ++index[i];
      for (std::size_t j = i + 1; j < s; ++j) {
        index[j] = index[j - 1] + 1;
      }
      return true;
    }
  }
  return false;
}

// The irreducible factors over Z of a squarefree primitive f of positive
// degree and positive leading coefficient, each primitive with positive
// leading coefficient.
inline std::vector<zpoly> zassenhaus(const zpoly& f) {
  auto& counters = factor_counters();
  ++counters.calls;

  const int n = degree(f);
  VARIETAS_ASSERT(n >= 1 && f.back() > 0);
  if (n == 1) {
    ++counters.true_factors;
    return {f};
  }

  // B is Mignotte's bound on the coefficients of a factor of f, scaled by
  // lc(f): sqrt(n + 1) 2^n |f|_inf lc(f), with the square root rounded up.
  const mpz_class& b = f.back();
  mpz_class A = 0;
  for (const mpz_class& c : f) {
    A = std::max<mpz_class>(A, abs(c));
  }
  mpz_class B = n + 1;
  mpz_sqrt(B.get_mpz_t(), B.get_mpz_t());
  B += 1;
  mpz_mul_2exp(B.get_mpz_t(), B.get_mpz_t(), static_cast<mp_bitcnt_t>(n));
  B *= A * b;

  // A prime that divides neither the leading coefficient nor the discriminant
  // keeps the degree and the squarefreeness. Of the first few, the one with
  // the fewest factors is taken, since the recombination is what costs; the
  // count needs only the distinct-degree pass.
  constexpr std::size_t candidates = 5;
  word best_prime = 0;
  std::vector<std::pair<md::upoly, int>> best;
  std::size_t best_count = 0;
  std::size_t accepted = 0;
  for (std::size_t index = 0; accepted < candidates; ++index) {
    const word p = md::nth_prime_below_2_31(index);
    if (md::reduce(b, p) == 0) {
      continue;
    }
    const md::prime_field F{p};
    const md::upoly image_f = md::monic(F, image(f, p));
    if (md::degree(md::gcd(F, image_f, derivative(F, image_f))) != 0) {
      continue;
    }
    ++accepted;
    ++counters.primes;
    auto split = distinct_degree(F, image_f);
    std::size_t count = 0;
    for (const auto& [g, d] : split) {
      count += static_cast<std::size_t>(md::degree(g) / d);
    }
    if (best_prime == 0 || count < best_count) {
      best_prime = p;
      best = std::move(split);
      best_count = count;
    }
    if (count == 1) {
      break;
    }
  }
  counters.modular_factors += best_count;
  if (best_count == 1) {
    ++counters.true_factors;
    return {f};
  }

  const md::prime_field F{best_prime};
  std::mt19937_64 rng(best_prime);
  std::vector<md::upoly> modular;
  for (const auto& [g, d] : best) {
    equal_degree(F, g, d, rng, modular);
  }

  // p^l > 2B.
  const mpz_class p = static_cast<unsigned long>(best_prime);
  mpz_class P = p;
  while (P <= 2 * B) {
    P *= p;
  }
  counters.precision_bits =
      std::max<std::size_t>(counters.precision_bits, mpz_sizeinbase(P.get_mpz_t(), 2));

  std::vector<zpoly> lifted;
  hensel_lift(reduce(f, P), modular, 0, modular.size(), F, P, lifted);

  // Recombination. Subsets are tried in increasing size, and a factor found
  // removes its members; what is left when no subset of up to half of the
  // remaining factors works is irreducible. A candidate's constant term has to
  // divide lc(f*) f*(0), which settles most subsets before any product of
  // polynomials is formed.
  std::vector<zpoly> out;
  zpoly rest = f;
  mpz_class lead = b;
  std::size_t s = 1;
  while (2 * s <= lifted.size()) {
    bool found = false;
    std::vector<std::size_t> index(s);
    std::iota(index.begin(), index.end(), std::size_t{0});
    const mpz_class half = P / 2;
    const mpz_class trailing = lead * rest.front();
    do {
      if (trailing != 0) {
        mpz_class constant = lead;
        for (std::size_t i : index) {
          constant = constant * lifted[i].front() % P;
        }
        if (constant > half) {
          constant -= P;
        }
        if (constant == 0 || trailing % constant != 0) {
          continue;
        }
      }
      ++counters.subsets;
      zpoly g{lead};
      zpoly h{lead};
      std::size_t next = 0;
      for (std::size_t i = 0; i < lifted.size(); ++i) {
        if (next < s && index[next] == i) {
          g = multiply(g, lifted[i], P);
          ++next;
        } else {
          h = multiply(h, lifted[i], P);
        }
      }
      g = symmetric(std::move(g), P);
      h = symmetric(std::move(h), P);
      if (norm_1(g) * norm_1(h) <= B) {
        out.push_back(primitive(std::move(g)));
        rest = primitive(std::move(h));
        lead = rest.back();
        for (std::size_t k = s; k-- > 0;) {
          lifted.erase(lifted.begin() + static_cast<std::ptrdiff_t>(index[k]));
        }
        found = true;
        break;
      }
    } while (next_combination(index, lifted.size()));
    if (!found) {
      ++s;
    }
  }
  out.push_back(std::move(rest));
  counters.true_factors += out.size();
  return out;
}

template <std::size_t N, class Order>
polynomial<rational, N, Order> derivative(const polynomial<rational, N, Order>& p, std::size_t v) {
  using poly = polynomial<rational, N, Order>;
  std::vector<typename poly::term> terms;
  for (const auto& t : p.terms()) {
    if (t.mon[v] == 0) {
      continue;
    }
    auto exponents = t.mon.exponents();
    const unsigned long k = exponents[v];
    --exponents[v];
    terms.push_back({monomial<N>(exponents), t.coeff * rational(k)});
  }
  return poly(std::move(terms));
}

// The content with respect to v, monic, by the modular gcd.
template <std::size_t N, class Order>
polynomial<rational, N, Order> content_in(const polynomial<rational, N, Order>& p, std::size_t v) {
  using poly = polynomial<rational, N, Order>;
  poly c;
  for (int k = detail::degree_in(p, v); k >= 0; --k) {
    const poly coefficient = detail::coefficient_in(p, v, k);
    if (coefficient.is_zero()) {
      continue;
    }
    c = c.is_zero() ? coefficient.monic() : modular_gcd(c, coefficient);
    if (c.degree() == 0) {
      break;
    }
  }
  return c;
}

template <std::size_t N, class Order>
void accumulate(std::map<unsigned, polynomial<rational, N, Order>>& parts, unsigned multiplicity,
                const polynomial<rational, N, Order>& g) {
  auto [it, inserted] = parts.try_emplace(multiplicity, g);
  if (!inserted) {
    it->second = it->second * g;
  }
}

// Yun's algorithm on f, primitive in v. With f = a_1 a_2^2 ... a_k^k, the
// gcd with the derivative takes one off every multiplicity, b_1 = f / that is
// a_1 ... a_k, and each further gcd peels off the factor of the current
// least multiplicity. The cofactors of each gcd are the next b and the next
// numerator, so no division is done separately.
template <std::size_t N, class Order>
void yun(const polynomial<rational, N, Order>& f, std::size_t v,
         std::map<unsigned, polynomial<rational, N, Order>>& parts) {
  auto first = modular_gcd_with_cofactors(f, derivative(f, v));
  auto b = std::move(first.a_over_gcd);
  auto d = first.b_over_gcd - derivative(b, v);
  for (unsigned i = 1; b.degree() > 0; ++i) {
    auto step = modular_gcd_with_cofactors(b, d);
    if (step.gcd.degree() > 0) {
      accumulate(parts, i, step.gcd);
    }
    b = std::move(step.a_over_gcd);
    d = step.b_over_gcd - derivative(b, v);
  }
}

// The content with respect to the first variable f involves is a polynomial
// in the others, decomposed by the same recursion with one variable fewer;
// the primitive part is decomposed by Yun. A factor of the content does not
// involve v and every factor of the primitive part does, so the two
// decompositions are coprime to each other and multiply together.
template <std::size_t N, class Order>
void squarefree_parts(const polynomial<rational, N, Order>& f,
                      std::map<unsigned, polynomial<rational, N, Order>>& parts) {
  if (f.degree() == 0) {
    return;
  }
  std::size_t v = 0;
  while (detail::degree_in(f, v) <= 0) {
    ++v;
  }
  const auto c = content_in(f, v);
  yun(detail::divide_exact(f, c), v, parts);
  squarefree_parts(c, parts);
}

}  // namespace factor_detail

// The squarefree decomposition of a nonzero f in any number of variables:
// f = unit * g_1 g_2^2 ... g_k^k with the g_i squarefree, monic under Order
// and pairwise coprime, one factor for each multiplicity that occurs, in
// increasing multiplicity.
template <std::size_t N, class Order>
factorisation<polynomial<rational, N, Order>> squarefree_decomposition(
    const polynomial<rational, N, Order>& f) {
  VARIETAS_ASSERT(!f.is_zero());
  std::map<unsigned, polynomial<rational, N, Order>> parts;
  factor_detail::squarefree_parts(f.monic(), parts);

  factorisation<polynomial<rational, N, Order>> out;
  out.unit = f.leading_coefficient();
  for (auto& [multiplicity, g] : parts) {
    out.factors.push_back({g.monic(), multiplicity});
  }
  return out;
}

// The factorisation into irreducibles over Q of a nonzero f involving at most
// one of the N variables. The factors are monic and distinct, ordered by
// multiplicity and then by degree, and the unit is the leading coefficient of
// f. A polynomial in more than one variable is a precondition violation, not
// a case: that is the multivariate factorisation this one will be lifted into.
template <std::size_t N, class Order>
factorisation<polynomial<rational, N, Order>> factor_univariate(
    const polynomial<rational, N, Order>& f) {
  using poly = polynomial<rational, N, Order>;
  using factor_detail::zpoly;
  VARIETAS_ASSERT(!f.is_zero());

  std::size_t v = N;
  for (std::size_t i = 0; i < N; ++i) {
    if (detail::degree_in(f, i) > 0) {
      VARIETAS_ASSERT(v == N);
      v = i;
    }
  }

  const auto squarefree = squarefree_decomposition(f);
  factorisation<poly> out;
  out.unit = squarefree.unit;
  for (const auto& [g, multiplicity] : squarefree.factors) {
    const auto form = modular_detail::primitive_integer(g);
    zpoly z(static_cast<std::size_t>(detail::degree_in(g, v)) + 1, 0);
    for (const auto& [m, c] : form.terms) {
      z[m[v]] = c;
    }
    if (z.back() < 0) {
      for (mpz_class& c : z) {
        c = -c;
      }
    }
    for (const zpoly& q : factor_detail::zassenhaus(z)) {
      std::vector<typename poly::term> terms;
      for (std::size_t k = 0; k < q.size(); ++k) {
        if (q[k] != 0) {
          terms.push_back(
              {monomial<N>::variable(v, static_cast<typename monomial<N>::exponent_type>(k)),
               rational(q[k])});
        }
      }
      out.factors.push_back({poly(std::move(terms)).monic(), multiplicity});
    }
  }
  std::stable_sort(out.factors.begin(), out.factors.end(), [](const auto& a, const auto& b) {
    if (a.multiplicity != b.multiplicity) {
      return a.multiplicity < b.multiplicity;
    }
    return a.base.degree() < b.base.degree();
  });
  return out;
}

}  // namespace varietas

#endif
