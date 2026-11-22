#ifndef VARIETAS_CODEGEN_MODULAR_GROEBNER_HPP
#define VARIETAS_CODEGEN_MODULAR_GROEBNER_HPP

#include <array>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <map>
#include <thread>
#include <utility>
#include <vector>

#include <gmpxx.h>

#include "varietas/codegen/modular_gcd.hpp"
#include "varietas/codegen/prime_field.hpp"
#include "varietas/codegen/rational.hpp"
#include "varietas/core/ideal/buchberger.hpp"
#include "varietas/core/ideal/division.hpp"
#include "varietas/core/ideal/saturation.hpp"
#include "varietas/core/monomial.hpp"
#include "varietas/core/polynomial.hpp"
#include "varietas/core/quotient/quotient_basis.hpp"

namespace varietas {

// Gröbner bases over Q from their images modulo primes.
//
// Buchberger over Q pays for the coefficients it carries on the way, and those
// grow far beyond the ones it finally returns: the saturated ideal of a 3R arm
// at a dyadic pose takes four seconds, nearly all of it in arithmetic on
// rationals the answer does not contain. Over Z/p a coefficient is one word,
// and F4 computes the same basis there in about a millisecond. So the basis
// over Q is recovered from its images, as the modular gcd recovers a gcd and
// the reconstruction a rational function.
//
// For each prime p the reduced basis of the ideal modulo p is computed. Its
// leading monomials are those of the basis over Q for all but finitely many
// primes, the unlucky ones, and those are found by disagreement, since the
// images are grouped by their leading monomials and only the largest group is
// used. Within it the coefficients are combined by the Chinese remainder
// theorem and read back by rational reconstruction. A candidate is accepted
// only after the next prime of its group has agreed with it, and then only if
// it passes a check over Q.
//
// The check over Q is exact, and it is half of what one would want:
//
//   * every S-pair of the candidate G reduces to zero over Q, so G is a
//     Gröbner basis of the ideal it spans; and
//   * every generator of I reduces to zero modulo G, so I ⊆ <G>.
//
// For a saturation I : h^∞ the second is not enough, since the saturation is
// larger than I, and it is completed by checking that h is a unit modulo <G>,
// that is, that <G> is already saturated: then I : h^∞ ⊆ <G> : h^∞ = <G>.
// For a zero-dimensional <G> that is the nonsingularity of multiplication by
// h on the finite quotient, and nonsingularity modulo any prime proves it over
// Q, since a determinant that is nonzero modulo p is nonzero.
//
// What is not checked is the reverse containment, <G> ⊆ I, which would need
// each element of G written over the generators. If it failed, <G> would be an
// ideal strictly containing I with the same leading monomials as I has modulo
// every prime in the majority, and a fresh prime would have agreed with it, which
// requires every one of those primes to be unlucky in the same way. That is why
// the answer is right with high probability rather than certainly, and why a
// caller that needs certainty uses groebner_basis over Q, which is Buchberger,
// and is unchanged.
//
// Whatever cannot be checked, a saturation that is not zero-dimensional for
// instance, and whatever runs out of primes, is computed by Buchberger
// instead, so this never returns less than the exact computation would, only
// sooner, and the statistics say which route was taken.
namespace modular {

struct statistics {
  std::size_t primes = 0;           // images used
  std::size_t wasted_images = 0;    // computed in the last batch after the answer
  std::size_t unusable_primes = 0;  // a denominator or leading coefficient vanished
  std::size_t outvoted_primes = 0;  // leading monomials other than the majority's
  std::size_t reconstructions = 0;  // attempts at rational reconstruction
  std::size_t candidates_rejected = 0;
  std::size_t coefficient_bits = 0;  // the largest numerator or denominator
  bool lifted = false;               // false if Buchberger computed it instead
  double image_seconds = 0.0;
  double lifting_seconds = 0.0;
  double verification_seconds = 0.0;
};

struct options {
  std::size_t max_primes = 4096;
  // Images modulo different primes are independent, and the prime is selected
  // per thread, so they are computed this many at a time. The answer does not
  // depend on it: the images are folded in the order of their primes, whatever
  // order they finished in.
  //
  // Zero is half the hardware threads, at most four. More is not faster: F4's
  // matrices are large and sparse, the images compete for memory bandwidth, and
  // on an eight-core laptop one image of a six-joint arm took 2.6 s alone and
  // the nineteen it needed took 21 s four at a time, 24 s eight at a time and 33 s
  // sixteen at a time.
  std::size_t threads = 0;
};

namespace detail {

using word = residue::word;
using clock = std::chrono::steady_clock;

inline double since(clock::time_point start) {
  return std::chrono::duration<double>(clock::now() - start).count();
}

// The image of f under the prime currently selected. The caller checks
// residue::undefined_image() for a denominator the prime divides.
template <std::size_t N, class Order>
polynomial<residue, N, Order> image(const polynomial<rational, N, Order>& f) {
  std::vector<typename polynomial<residue, N, Order>::term> terms;
  terms.reserve(f.size());
  for (const auto& t : f.terms()) {
    terms.push_back({t.mon, residue(t.coeff)});
  }
  return polynomial<residue, N, Order>(std::move(terms));
}

// The images of a family, false if the prime is unusable for it: it divides
// a denominator, or it kills a leading coefficient, which changes the
// leading monomial of an input and with it, in general, the basis.
template <std::size_t N, class Order>
bool images(const std::vector<polynomial<rational, N, Order>>& family,
            std::vector<polynomial<residue, N, Order>>& out) {
  out.clear();
  out.reserve(family.size());
  for (const auto& f : family) {
    out.push_back(image(f));
    if (residue::undefined_image()) {
      return false;
    }
    if (!f.is_zero() && (out.back().is_zero() ||
                         out.back().leading_monomial() != f.leading_monomial())) {
      return false;
    }
  }
  return true;
}

template <std::size_t N>
using exponents = std::array<std::uint16_t, N>;

// The primes whose images share one set of leading monomials, and what they
// have accumulated.
template <std::size_t N, class Order>
struct group {
  std::vector<exponents<N>> leading;
  std::vector<std::map<exponents<N>, mpz_class>> coefficients;  // per element
  mpz_class modulus = 1;
  std::size_t primes = 0;
  std::size_t next_attempt = 1;
  std::vector<polynomial<rational, N, Order>> candidate;
  bool have_candidate = false;
};

// True when the candidate maps to exactly the image modulo the current prime.
template <std::size_t N, class Order>
bool agrees(const std::vector<polynomial<rational, N, Order>>& candidate,
            const std::vector<polynomial<residue, N, Order>>& image_basis) {
  if (candidate.size() != image_basis.size()) {
    return false;
  }
  residue::clear_undefined_image();
  for (std::size_t e = 0; e < candidate.size(); ++e) {
    if (image(candidate[e]) != image_basis[e] || residue::undefined_image()) {
      return false;
    }
  }
  return true;
}

// Folds one image into a group, as in the modular gcd.
template <std::size_t N, class Order>
void combine(group<N, Order>& g, const std::vector<polynomial<residue, N, Order>>& basis) {
  const word p = residue::prime();
  const modular_detail::prime_field F{p};
  if (g.coefficients.empty()) {
    g.coefficients.resize(basis.size());
  }
  const word m_inverse = F.inv(modular_detail::reduce(g.modulus, p));
  for (std::size_t e = 0; e < basis.size(); ++e) {
    auto& into = g.coefficients[e];
    std::map<exponents<N>, word> now;
    for (const auto& t : basis[e].terms()) {
      now.emplace(t.mon.exponents(), t.coeff.value());
      into.emplace(t.mon.exponents(), mpz_class(0));
    }
    for (auto& [mon, r] : into) {
      const auto it = now.find(mon);
      const word target = it == now.end() ? 0 : it->second;
      const word delta = F.mul(F.sub(target, modular_detail::reduce(r, p)), m_inverse);
      r += g.modulus * static_cast<unsigned long>(delta);
    }
  }
  g.modulus *= static_cast<unsigned long>(p);
  ++g.primes;
}

template <std::size_t N, class Order>
bool reconstruct(group<N, Order>& g, std::size_t& bits) {
  using poly = polynomial<rational, N, Order>;
  std::vector<poly> out;
  out.reserve(g.coefficients.size());
  std::size_t largest = 0;
  for (const auto& element : g.coefficients) {
    std::vector<typename poly::term> terms;
    terms.reserve(element.size());
    for (const auto& [mon, r] : element) {
      rational c;
      if (!modular_detail::rational_reconstruction(r, g.modulus, c)) {
        return false;
      }
      if (sgn(c) != 0) {
        largest = std::max({largest, mpz_sizeinbase(c.get_num().get_mpz_t(), 2),
                            mpz_sizeinbase(c.get_den().get_mpz_t(), 2)});
        terms.push_back({monomial<N>(mon), c});
      }
    }
    out.emplace_back(std::move(terms));
  }
  g.candidate = std::move(out);
  bits = largest;
  return true;
}

// The exact half of the check: G is a Gröbner basis, and every generator lies
// in the ideal it spans.
template <std::size_t N, class Order>
bool is_basis_containing(const std::vector<polynomial<rational, N, Order>>& G,
                         const std::vector<polynomial<rational, N, Order>>& generators) {
  using mono = monomial<N>;
  for (const auto& f : generators) {
    if (!normal_form(f, G).is_zero()) {
      return false;
    }
  }
  for (std::size_t i = 0; i < G.size(); ++i) {
    for (std::size_t j = i + 1; j < G.size(); ++j) {
      const mono& a = G[i].leading_monomial();
      const mono& b = G[j].leading_monomial();
      if (mono::lcm(a, b) == a * b) {
        continue;  // Buchberger's first criterion
      }
      if (!normal_form(s_polynomial(G[i], G[j]), G).is_zero()) {
        return false;
      }
    }
  }
  return true;
}

// Whether multiplication by h is invertible on the quotient by a
// zero-dimensional <G>, decided modulo a prime: the matrix is exact over Q,
// and a determinant nonzero modulo p is nonzero. A prime that divides a
// denominator, or at which the matrix happens to be singular, is passed over
// for the next; false only if none of a few dozen show it nonsingular.
template <std::size_t N, class Order>
bool is_unit_modulo(const polynomial<rational, N, Order>& h,
                    const std::vector<polynomial<rational, N, Order>>& G,
                    const quotient_basis<N>& quotient) {
  using poly = polynomial<rational, N, Order>;
  const std::size_t d = quotient.dimension();
  std::vector<std::vector<rational>> matrix(d, std::vector<rational>(d, rational(0)));
  for (std::size_t j = 0; j < d; ++j) {
    const poly image = normal_form(h * poly::from_monomial(quotient.monomials[j], rational(1)), G);
    for (const auto& t : image.terms()) {
      const std::size_t row = quotient.index_of(t.mon);
      VARIETAS_ASSERT(row < d);
      matrix[row][j] = t.coeff;
    }
  }
  const word previous = residue::prime();
  bool unit = false;
  for (std::size_t index = 0; index < 32 && !unit; ++index) {
    const word p = modular_detail::nth_prime_below_2_31(index);
    residue::select(p);
    residue::clear_undefined_image();
    const modular_detail::prime_field F{p};
    std::vector<std::vector<word>> a(d, std::vector<word>(d));
    for (std::size_t r = 0; r < d; ++r) {
      for (std::size_t c = 0; c < d; ++c) {
        a[r][c] = residue(matrix[r][c]).value();
      }
    }
    if (residue::undefined_image()) {
      continue;
    }
    bool singular = false;
    for (std::size_t c = 0; c < d && !singular; ++c) {
      std::size_t pivot = c;
      while (pivot < d && a[pivot][c] == 0) {
        ++pivot;
      }
      if (pivot == d) {
        singular = true;
        break;
      }
      std::swap(a[c], a[pivot]);
      const word inverse = F.inv(a[c][c]);
      for (std::size_t r = c + 1; r < d; ++r) {
        if (a[r][c] == 0) {
          continue;
        }
        const word factor = F.mul(a[r][c], inverse);
        for (std::size_t k = c; k < d; ++k) {
          a[r][k] = F.sub(a[r][k], F.mul(factor, a[c][k]));
        }
      }
    }
    unit = !singular;
  }
  residue::select(previous);
  residue::clear_undefined_image();
  return unit;
}

enum class verdict { accepted, rejected, unverifiable };

template <std::size_t N, class Order>
using image_function = std::function<bool(std::vector<polynomial<residue, N, Order>>&)>;

template <std::size_t N, class Order>
using check_function = std::function<verdict(const std::vector<polynomial<rational, N, Order>>&)>;

// The lifting loop, for any ideal whose reduced basis can be computed modulo
// a prime. False if the primes ran out or the check could not be made, and the
// caller then computes the basis exactly.
template <std::size_t N, class Order>
bool lift(const image_function<N, Order>& image_of, const check_function<N, Order>& check,
          std::vector<polynomial<rational, N, Order>>& out, statistics& s, const options& opts) {
  std::vector<group<N, Order>> groups;
  std::vector<polynomial<residue, N, Order>> basis;
  bool result = false;
  const word original = residue::prime();

  std::size_t threads = opts.threads;
  if (threads == 0) {
    threads = std::min<std::size_t>(4, std::thread::hardware_concurrency() / 2);
  }
  threads = std::max<std::size_t>(1, threads);
  struct computed {
    word prime = 0;
    bool usable = false;
    std::vector<polynomial<residue, N, Order>> basis;
  };
  std::vector<computed> batch;
  std::size_t next_in_batch = 0;
  std::size_t next_index = 0;

  for (;;) {
    if (next_in_batch == batch.size()) {
      if (next_index >= opts.max_primes) {
        break;
      }
      // The next batch, one prime per thread. The primes are found here first,
      // since the list of them is shared and filled on demand.
      const std::size_t size = std::min(threads, opts.max_primes - next_index);
      batch.assign(size, computed{});
      for (std::size_t k = 0; k < size; ++k) {
        batch[k].prime = modular_detail::nth_prime_below_2_31(next_index + k);
      }
      next_index += size;
      next_in_batch = 0;
      const auto image_start = clock::now();
      const auto work = [&](std::size_t k) {
        residue::select(batch[k].prime);
        residue::clear_undefined_image();
        batch[k].usable = image_of(batch[k].basis) && !residue::undefined_image();
      };
      if (size == 1) {
        work(0);
      } else {
        std::vector<std::thread> pool;
        pool.reserve(size);
        for (std::size_t k = 0; k < size; ++k) {
          pool.emplace_back(work, k);
        }
        for (std::thread& t : pool) {
          t.join();
        }
      }
      s.image_seconds += since(image_start);
    }
    computed& item = batch[next_in_batch++];
    if (!item.usable) {
      ++s.unusable_primes;
      continue;
    }
    ++s.primes;
    const word p = item.prime;
    residue::select(p);
    residue::clear_undefined_image();
    basis = std::move(item.basis);

    std::vector<exponents<N>> leading;
    leading.reserve(basis.size());
    for (const auto& g : basis) {
      leading.push_back(g.leading_monomial().exponents());
    }
    std::size_t gi = 0;
    while (gi < groups.size() && groups[gi].leading != leading) {
      ++gi;
    }
    if (gi == groups.size()) {
      groups.emplace_back();
      groups.back().leading = leading;
    }
    group<N, Order>& g = groups[gi];

    const auto lifting_start = clock::now();
    // A candidate meets its first prime it was not built from here.
    if (g.have_candidate) {
      g.have_candidate = false;
      residue::select(p);
      if (agrees(g.candidate, basis)) {
        s.lifting_seconds += since(lifting_start);
        const auto check_start = clock::now();
        const verdict v = check(g.candidate);
        s.verification_seconds += since(check_start);
        if (v == verdict::accepted) {
          out = g.candidate;
          ++g.primes;  // the prime that confirmed it, though not folded in
          result = true;
          break;
        }
        if (v == verdict::unverifiable) {
          break;
        }
      }
      ++s.candidates_rejected;
      residue::select(p);
      residue::clear_undefined_image();
    }

    combine(g, basis);
    bool majority = true;
    for (const auto& other : groups) {
      majority = majority && other.primes <= g.primes;
    }
    if (majority && g.primes >= g.next_attempt) {
      ++s.reconstructions;
      g.have_candidate = reconstruct(g, s.coefficient_bits);
      g.next_attempt = g.primes + std::max<std::size_t>(1, g.primes / 4);
    }
    s.lifting_seconds += since(lifting_start);
  }

  s.wasted_images = batch.size() - next_in_batch;
  std::size_t largest = 0;
  for (const auto& g : groups) {
    largest = std::max(largest, g.primes);
  }
  s.outvoted_primes = s.primes - largest;
  residue::select(original);
  residue::clear_undefined_image();
  return result;
}

}  // namespace detail

// The reduced Gröbner basis of the ideal the generators span, lifted from
// Z/p and checked over Q as described above; by Buchberger when it cannot be.
template <std::size_t N, class Order>
std::vector<polynomial<rational, N, Order>> groebner_basis(
    const std::vector<polynomial<rational, N, Order>>& generators,
    statistics* stats = nullptr, const options& opts = options()) {
  using poly = polynomial<rational, N, Order>;
  statistics local;
  statistics& s = stats != nullptr ? *stats : local;

  std::vector<poly> nonzero;
  for (const poly& f : generators) {
    if (!f.is_zero()) {
      nonzero.push_back(f);
    }
  }
  if (nonzero.empty()) {
    return {};
  }

  const detail::image_function<N, Order> image_of =
      [&](std::vector<polynomial<residue, N, Order>>& out) {
        std::vector<polynomial<residue, N, Order>> mapped;
        if (!detail::images(nonzero, mapped)) {
          return false;
        }
        out = varietas::groebner_basis(mapped);
        return true;
      };
  const detail::check_function<N, Order> check = [&](const std::vector<poly>& G) {
    return detail::is_basis_containing(G, nonzero) ? detail::verdict::accepted
                                                   : detail::verdict::rejected;
  };

  std::vector<poly> out;
  if (detail::lift(image_of, check, out, s, opts)) {
    s.lifted = true;
    return out;
  }
  return varietas::groebner_basis(nonzero);
}

// I : h^∞ for I spanned by the generators and h the product of the
// multipliers, as a reduced basis under Order: the saturation the half-angle
// substitution calls for (varietas/core/ideal/saturation.hpp), lifted from
// Z/p. It is checked as described above when the result is zero-dimensional,
// and computed by Buchberger when it is not.
template <std::size_t N, class Order>
std::vector<polynomial<rational, N, Order>> saturate_by_product(
    const std::vector<polynomial<rational, N, Order>>& generators,
    const std::vector<polynomial<rational, N, Order>>& multipliers,
    statistics* stats = nullptr, const options& opts = options()) {
  using poly = polynomial<rational, N, Order>;
  statistics local;
  statistics& s = stats != nullptr ? *stats : local;

  std::vector<poly> nonzero;
  for (const poly& f : generators) {
    if (!f.is_zero()) {
      nonzero.push_back(f);
    }
  }
  if (nonzero.empty()) {
    return {};
  }
  poly h = poly::constant(rational(1));
  for (const poly& m : multipliers) {
    VARIETAS_ASSERT(!m.is_zero());
    h = h * m;
  }

  const detail::image_function<N, Order> image_of =
      [&](std::vector<polynomial<residue, N, Order>>& out) {
        std::vector<polynomial<residue, N, Order>> mapped;
        std::vector<polynomial<residue, N, Order>> mapped_h;
        if (!detail::images(nonzero, mapped) ||
            !detail::images(std::vector<poly>{h}, mapped_h)) {
          return false;
        }
        out = varietas::saturate(mapped, mapped_h.front());
        return true;
      };
  const detail::check_function<N, Order> check = [&](const std::vector<poly>& G) {
    const quotient_basis<N> quotient = standard_monomials(G);
    if (!quotient.is_zero_dimensional) {
      return detail::verdict::unverifiable;
    }
    if (quotient.dimension() == 0) {
      // G = {1}, the unit ideal, which every ideal is contained in: the
      // check is vacuous, and an unreachable pose is exactly as certain as
      // the reverse containment, no more.
      return detail::verdict::accepted;
    }
    if (!detail::is_basis_containing(G, nonzero)) {
      return detail::verdict::rejected;
    }
    return detail::is_unit_modulo(h, G, quotient) ? detail::verdict::accepted
                                                  : detail::verdict::rejected;
  };

  std::vector<poly> out;
  if (detail::lift(image_of, check, out, s, opts)) {
    s.lifted = true;
    return out;
  }
  return varietas::saturate_by_product(nonzero, multipliers);
}

}  // namespace modular
}  // namespace varietas

#endif
