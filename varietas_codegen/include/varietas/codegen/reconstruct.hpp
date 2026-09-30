#ifndef VARIETAS_CODEGEN_RECONSTRUCT_HPP
#define VARIETAS_CODEGEN_RECONSTRUCT_HPP

#include <algorithm>
#include <array>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <map>
#include <random>
#include <utility>
#include <vector>

#include <gmpxx.h>

#include "varietas/codegen/modular_gcd.hpp"
#include "varietas/codegen/prime_field.hpp"
#include "varietas/codegen/rational.hpp"
#include "varietas/codegen/rational_function.hpp"
#include "varietas/core/config.hpp"
#include "varietas/core/monomial.hpp"
#include "varietas/core/order/grevlex.hpp"
#include "varietas/core/polynomial.hpp"

namespace varietas {
namespace reconstruction {

// Rational functions over Q recovered from their values at points of Z/p.
//
// The parametric solve computes rational functions of the pose by carrying
// them symbolically through Buchberger's algorithm, and on arms with offsets
// the functions it carries on the way swell to thousands of terms although the
// ones it finally produces have tens. Reconstruction never forms the
// intermediate ones. What it needs is a black box that, given a pose with
// coordinates in Z/p, returns the value of every wanted function there; for
// the inverse kinematics that is one Grobner basis at a fixed pose, over a
// prime field, with no parameters. From enough such values each function is
// recovered, first over Z/p and then, from several primes, over Q.
//
// The steps, for each prime:
//
//   1. Degrees. The functions are sampled along a random line through the
//      parameter space, where each is a univariate rational function of the
//      line's parameter whose numerator and denominator degrees are, for a
//      line in general position, the total degrees of the multivariate ones.
//      Those are read off by rational reconstruction of the interpolating
//      polynomial, choosing the step of the extended Euclidean algorithm whose
//      quotient has the largest degree (Monagan's maximal quotient rule), and
//      the number of points is raised until every function is predicted
//      correctly at two points it was not built from.
//
//   2. Coefficients. With the degrees known, N(p) - v D(p) = 0 at a sample is
//      one linear equation in the coefficients of N and D, dense over all
//      monomials of those total degrees. For N/D in lowest terms of exactly
//      those degrees the solutions are the multiples of one vector, so once
//      there are enough samples the nullspace is a line, and it is normalised
//      so that the grevlex-leading coefficient of D is one, which is how the
//      result is normalised over Q and so the same vector modulo every prime
//      that is not unlucky.
//
// Across primes the coefficients are combined by the Chinese remainder theorem
// and read back by rational reconstruction, as in the modular gcd, until two
// consecutive primes leave every reconstruction unchanged.
//
// Unlike the modular gcd this has no final exact test, because the functions
// are only available as a black box, and the answer is therefore right with
// high probability rather than certainly: a wrong one requires a sample, a
// line or a prime to fall on the zero set of some fixed nonzero polynomial,
// with every choice made from a set of size about 2^31. The caller is expected
// to check the result against an exact computation at a few rational points,
// which is what reconstructed_ik does.

using word = std::uint64_t;

template <std::size_t P>
using point = std::array<residue, P>;

// Values of every function at a point, for the prime currently selected. False
// if the point is unusable, which the reconstruction treats as a draw to be
// thrown away: a pole, or a pose whose solution set has a different shape
// from the general one.
template <std::size_t P>
using black_box = std::function<bool(const point<P>&, std::vector<residue>&)>;

// Called once per prime, before any sample under it, with the prime already
// selected. False if the prime is unusable for the problem, such as one that
// divides a denominator of the input.
using prime_hook = std::function<bool(word)>;

struct statistics {
  std::size_t primes = 0;
  std::size_t samples = 0;           // black-box calls that returned values
  std::size_t rejected_samples = 0;  // and those that did not
  std::size_t line_samples = 0;      // of the first, the ones spent on degrees
  int largest_numerator_degree = 0;
  int largest_denominator_degree = 0;
  std::size_t largest_system = 0;  // unknowns in one linear solve
  double sampling_seconds = 0.0;   // inside the black box
  double solving_seconds = 0.0;    // in the linear systems
};

struct options {
  std::size_t max_primes = 64;
  std::size_t max_line_points = 400;
  std::uint64_t seed = 0x5eed;
};

namespace detail {

using modular_detail::upoly;

inline modular_detail::prime_field field() { return modular_detail::prime_field{residue::prime()}; }

// Rational reconstruction of the polynomial interpolating (t_i, v_i), by the
// maximal quotient rule. False when no step stands out, which means there are
// not yet enough points.
inline bool univariate(const std::vector<word>& t, const std::vector<word>& v, upoly& numerator,
                       upoly& denominator) {
  const auto F = field();
  // Newton interpolation.
  upoly u;
  upoly basis{1};
  for (std::size_t i = 0; i < t.size(); ++i) {
    const word at = modular_detail::evaluate(F, u, t[i]);
    const word scale = F.mul(F.sub(v[i], at), F.inv(modular_detail::evaluate(F, basis, t[i])));
    u = modular_detail::add(F, u, modular_detail::scale(F, basis, scale));
    basis = modular_detail::multiply(F, basis, upoly{F.neg(t[i]), 1});
  }
  if (u.empty()) {
    numerator.clear();
    denominator = upoly{1};
    return true;
  }
  upoly r0 = basis;  // the product of (x - t_i)
  upoly r1 = u;
  upoly s0;
  upoly s1{1};
  int best = 0;
  while (!r1.empty()) {
    auto [q, r2] = modular_detail::divide(F, r0, r1);
    if (modular_detail::degree(q) > best) {
      best = modular_detail::degree(q);
      numerator = r1;
      denominator = s1;
    }
    upoly s2 = modular_detail::add(
        F, s0, modular_detail::scale(F, modular_detail::multiply(F, q, s1), F.neg(1)));
    r0 = std::move(r1);
    r1 = std::move(r2);
    s0 = std::move(s1);
    s1 = std::move(s2);
  }
  // Two points beyond what the degrees need, at least, before the answer is
  // believed.
  if (best < 3) {
    return false;
  }
  const word lead = F.inv(denominator.back());
  numerator = modular_detail::scale(F, numerator, lead);
  denominator = modular_detail::scale(F, denominator, lead);
  return true;
}

// Every exponent vector of total degree at most d, in decreasing grevlex order.
template <std::size_t P>
std::vector<monomial<P>> monomials_up_to(int d) {
  std::vector<monomial<P>> out;
  std::array<typename monomial<P>::exponent_type, P> e{};
  std::function<void(std::size_t, int)> fill = [&](std::size_t v, int left) {
    if (v == P) {
      out.push_back(monomial<P>(e));
      return;
    }
    for (int k = 0; k <= left; ++k) {
      e[v] = static_cast<typename monomial<P>::exponent_type>(k);
      fill(v + 1, left - k);
    }
    e[v] = 0;
  };
  fill(0, d);
  std::sort(out.begin(), out.end(), [](const monomial<P>& a, const monomial<P>& b) {
    return grevlex::compare(a, b) > 0;
  });
  return out;
}

template <std::size_t P>
word evaluate_monomial(const monomial<P>& m, const point<P>& at) {
  residue value = residue::raw(1);
  for (std::size_t v = 0; v < P; ++v) {
    for (unsigned k = 0; k < m[v]; ++k) {
      value = value * at[v];
    }
  }
  return value.value();
}

// The nullspace of an r x c matrix over Z/p when it is one-dimensional, as a
// vector spanning it; the caller normalises. Empty when the nullity is not one.
inline std::vector<word> kernel_line(std::vector<std::vector<word>> a, std::size_t columns) {
  const auto F = field();
  std::vector<std::size_t> pivot_of_row;
  std::vector<int> pivot_row(columns, -1);
  std::size_t row = 0;
  for (std::size_t c = 0; c < columns && row < a.size(); ++c) {
    std::size_t found = row;
    while (found < a.size() && a[found][c] == 0) {
      ++found;
    }
    if (found == a.size()) {
      continue;
    }
    std::swap(a[row], a[found]);
    const word inverse = F.inv(a[row][c]);
    for (std::size_t k = c; k < columns; ++k) {
      a[row][k] = F.mul(a[row][k], inverse);
    }
    for (std::size_t i = 0; i < a.size(); ++i) {
      if (i == row || a[i][c] == 0) {
        continue;
      }
      const word factor = a[i][c];
      for (std::size_t k = c; k < columns; ++k) {
        if (a[row][k] != 0) {
          a[i][k] = F.sub(a[i][k], F.mul(factor, a[row][k]));
        }
      }
    }
    pivot_row[c] = static_cast<int>(row);
    ++row;
  }
  if (row + 1 != columns) {
    return {};  // nullity other than one
  }
  std::size_t free = columns;
  for (std::size_t c = 0; c < columns; ++c) {
    if (pivot_row[c] < 0) {
      free = c;
      break;
    }
  }
  std::vector<word> x(columns, 0);
  x[free] = 1;
  for (std::size_t c = 0; c < columns; ++c) {
    if (pivot_row[c] >= 0) {
      x[c] = F.neg(a[static_cast<std::size_t>(pivot_row[c])][free]);
    }
  }
  return x;
}

// One function modulo the current prime: coefficients over the monomials of
// its degree bounds, normalised so that the leading coefficient of the
// denominator is one.
template <std::size_t P>
struct modular_image {
  bool zero = false;
  int numerator_degree = 0;
  int denominator_degree = 0;
  std::map<std::array<std::uint16_t, P>, word> numerator;
  std::map<std::array<std::uint16_t, P>, word> denominator;
  std::array<std::uint16_t, P> denominator_leading{};
};

}  // namespace detail

// Reconstructs `count` rational functions of P parameters from the black box.
// False if the budget of primes or points ran out first.
template <std::size_t P>
bool reconstruct(const black_box<P>& black, const prime_hook& prepare, std::size_t count,
                 std::vector<rational_function<P>>& out, statistics* stats = nullptr,
                 const options& opts = options()) {
  using namespace detail;
  using poly = polynomial<rational, P, grevlex>;
  statistics local;
  statistics& s = stats != nullptr ? *stats : local;
  using clock = std::chrono::steady_clock;
  const auto box = [&](const point<P>& x, std::vector<residue>& values) {
    const auto start = clock::now();
    const bool ok = black(x, values);
    s.sampling_seconds += std::chrono::duration<double>(clock::now() - start).count();
    return ok;
  };

  std::vector<int> numerator_degree(count, -1);
  std::vector<int> denominator_degree(count, -1);
  bool degrees_known = false;
  std::vector<std::vector<monomial<P>>> support_numerator(count);
  std::vector<std::vector<monomial<P>>> support_denominator(count);
  bool have_support = false;

  // Residues combined so far, per function.
  struct accumulated {
    std::map<std::array<std::uint16_t, P>, mpz_class> numerator;
    std::map<std::array<std::uint16_t, P>, mpz_class> denominator;
    std::array<std::uint16_t, P> denominator_leading{};
    bool zero = false;
  };
  std::vector<accumulated> acc(count);
  mpz_class modulus = 1;
  std::vector<rational_function<P>> previous;
  bool have_previous = false;

  std::mt19937_64 rng(opts.seed);

  for (std::size_t index = 0; index < opts.max_primes; ++index) {
    const word p = modular_detail::nth_prime_below_2_31(index);
    residue::select(p);
    residue::clear_undefined_image();
    if (!prepare(p) || residue::undefined_image()) {
      continue;
    }
    const auto F = field();
    std::uniform_int_distribution<word> draw(1, p - 1);
    const auto random_point = [&] {
      point<P> x;
      for (auto& c : x) {
        c = residue::raw(draw(rng));
      }
      return x;
    };

    std::vector<point<P>> points;
    std::vector<std::vector<residue>> values;
    std::vector<residue> buffer;

    // Step 1, on the first usable prime: degrees from a line.
    if (!degrees_known) {
      const point<P> base = random_point();
      const point<P> direction = random_point();
      std::vector<word> ts;
      std::vector<std::vector<word>> line_values(count);
      std::size_t next_t = 1;
      bool settled = false;
      while (!settled && ts.size() < opts.max_line_points) {
        const word t = next_t++;
        point<P> x;
        for (std::size_t v = 0; v < P; ++v) {
          x[v] = base[v] + residue::raw(t) * direction[v];
        }
        if (!box(x, buffer)) {
          ++s.rejected_samples;
          continue;
        }
        ++s.samples;
        ++s.line_samples;
        // Kept out of the linear systems below: points on one line do not
        // determine a polynomial in P variables, and a system filled with
        // them has a nullspace larger than a line.
        ts.push_back(t);
        for (std::size_t f = 0; f < count; ++f) {
          line_values[f].push_back(buffer[f].value());
        }
        if (ts.size() < 4 || ts.size() % 2 != 0) {
          continue;
        }
        settled = true;
        for (std::size_t f = 0; f < count && settled; ++f) {
          upoly n;
          upoly d;
          if (!univariate(ts, line_values[f], n, d)) {
            settled = false;
            break;
          }
          numerator_degree[f] = n.empty() ? -1 : modular_detail::degree(n);
          denominator_degree[f] = modular_detail::degree(d);
        }
      }
      if (!settled) {
        return false;
      }
      degrees_known = true;
      for (std::size_t f = 0; f < count; ++f) {
        s.largest_numerator_degree = std::max(s.largest_numerator_degree, numerator_degree[f]);
        s.largest_denominator_degree =
            std::max(s.largest_denominator_degree, denominator_degree[f]);
      }
    }

    // Step 2: enough points for the largest system, and a few over.
    //
    // On the first prime the unknowns are every monomial of the degree bounds.
    // After it, only the monomials that prime found a nonzero coefficient for:
    // a coefficient that is nonzero over Q is nonzero modulo all but finitely
    // many primes, so the support does not change, and the systems, and the
    // number of samples they need, shrink to the size of the answer. A term
    // the first prime missed by vanishing there makes the smaller system
    // inconsistent, the nullspace comes out empty, and the dense systems are
    // used again from the next prime on.
    std::size_t needed = 0;
    std::vector<std::vector<monomial<P>>> numerator_monomials(count);
    std::vector<std::vector<monomial<P>>> denominator_monomials(count);
    for (std::size_t f = 0; f < count; ++f) {
      if (numerator_degree[f] < 0) {
        continue;  // identically zero along a general line, so zero
      }
      if (have_support) {
        numerator_monomials[f] = support_numerator[f];
        denominator_monomials[f] = support_denominator[f];
      } else {
        numerator_monomials[f] = monomials_up_to<P>(numerator_degree[f]);
        denominator_monomials[f] = monomials_up_to<P>(denominator_degree[f]);
      }
      const std::size_t unknowns = numerator_monomials[f].size() + denominator_monomials[f].size();
      needed = std::max(needed, unknowns + 4);
      s.largest_system = std::max(s.largest_system, unknowns);
    }
    while (points.size() < needed) {
      const point<P> x = random_point();
      if (!box(x, buffer)) {
        ++s.rejected_samples;
        continue;
      }
      ++s.samples;
      points.push_back(x);
      values.push_back(buffer);
    }

    const auto solve_start = clock::now();
    std::vector<modular_image<P>> images(count);
    bool prime_ok = true;
    for (std::size_t f = 0; f < count && prime_ok; ++f) {
      modular_image<P>& image = images[f];
      if (numerator_degree[f] < 0) {
        image.zero = true;
        continue;
      }
      const auto& nm = numerator_monomials[f];
      const auto& dm = denominator_monomials[f];
      const std::size_t columns = nm.size() + dm.size();
      const std::size_t rows = columns + 4;
      std::vector<std::vector<word>> a(rows, std::vector<word>(columns, 0));
      for (std::size_t r = 0; r < rows; ++r) {
        const word v = values[r][f].value();
        for (std::size_t c = 0; c < nm.size(); ++c) {
          a[r][c] = evaluate_monomial(nm[c], points[r]);
        }
        for (std::size_t c = 0; c < dm.size(); ++c) {
          a[r][nm.size() + c] = F.neg(F.mul(v, evaluate_monomial(dm[c], points[r])));
        }
      }
      const std::vector<word> x = kernel_line(std::move(a), columns);
      if (x.empty()) {
        prime_ok = false;
        break;
      }
      // dm is in decreasing grevlex order, so the first nonzero is the lead.
      std::size_t lead = dm.size();
      for (std::size_t c = 0; c < dm.size(); ++c) {
        if (x[nm.size() + c] != 0) {
          lead = c;
          break;
        }
      }
      if (lead == dm.size()) {
        prime_ok = false;
        break;
      }
      const word scale = F.inv(x[nm.size() + lead]);
      image.denominator_leading = dm[lead].exponents();
      for (std::size_t c = 0; c < nm.size(); ++c) {
        if (x[c] != 0) {
          image.numerator[nm[c].exponents()] = F.mul(x[c], scale);
        }
      }
      for (std::size_t c = 0; c < dm.size(); ++c) {
        if (x[nm.size() + c] != 0) {
          image.denominator[dm[c].exponents()] = F.mul(x[nm.size() + c], scale);
        }
      }
    }
    s.solving_seconds += std::chrono::duration<double>(clock::now() - solve_start).count();
    if (!prime_ok) {
      have_support = false;  // back to the dense systems
      continue;              // an unlucky prime or unlucky samples; try the next
    }
    if (!have_support) {
      // Kept in the order the monomials were offered in, decreasing grevlex,
      // which is what finds the leading coefficient of a denominator.
      for (std::size_t f = 0; f < count; ++f) {
        support_numerator[f].clear();
        support_denominator[f].clear();
        for (const auto& m : numerator_monomials[f]) {
          if (images[f].numerator.count(m.exponents()) != 0) {
            support_numerator[f].push_back(m);
          }
        }
        for (const auto& m : denominator_monomials[f]) {
          if (images[f].denominator.count(m.exponents()) != 0) {
            support_denominator[f].push_back(m);
          }
        }
      }
      have_support = true;
    }
    ++s.primes;

    // Chinese remaindering, as in the modular gcd.
    const word m_inverse = F.inv(modular_detail::reduce(modulus, p));
    const auto combine = [&](std::map<std::array<std::uint16_t, P>, mpz_class>& into,
                             const std::map<std::array<std::uint16_t, P>, word>& image) {
      for (const auto& [e, c] : image) {
        into.emplace(e, mpz_class(0));
      }
      for (auto& [e, r] : into) {
        const auto it = image.find(e);
        const word target = it == image.end() ? 0 : it->second;
        const word delta = F.mul(F.sub(target, modular_detail::reduce(r, p)), m_inverse);
        r += modulus * static_cast<unsigned long>(delta);
      }
    };
    bool structure_changed = false;
    for (std::size_t f = 0; f < count; ++f) {
      if (images[f].zero) {
        acc[f].zero = true;
        continue;
      }
      if (modulus != 1 && images[f].denominator_leading != acc[f].denominator_leading) {
        structure_changed = true;
        break;
      }
    }
    if (structure_changed) {
      continue;  // this prime disagrees on a leading monomial: unlucky
    }
    for (std::size_t f = 0; f < count; ++f) {
      if (images[f].zero) {
        continue;
      }
      acc[f].denominator_leading = images[f].denominator_leading;
      combine(acc[f].numerator, images[f].numerator);
      combine(acc[f].denominator, images[f].denominator);
    }
    modulus *= static_cast<unsigned long>(p);

    // Rational reconstruction of every coefficient.
    std::vector<rational_function<P>> current;
    current.reserve(count);
    bool reconstructed = true;
    const auto lift = [&](const std::map<std::array<std::uint16_t, P>, mpz_class>& from,
                          poly& to) {
      std::vector<typename poly::term> terms;
      for (const auto& [e, r] : from) {
        rational c;
        if (!modular_detail::rational_reconstruction(r, modulus, c)) {
          return false;
        }
        if (sgn(c) != 0) {
          terms.push_back({monomial<P>(e), c});
        }
      }
      to = poly(std::move(terms));
      return true;
    };
    for (std::size_t f = 0; f < count && reconstructed; ++f) {
      if (acc[f].zero) {
        current.emplace_back(rational(0));
        continue;
      }
      poly n;
      poly d;
      if (!lift(acc[f].numerator, n) || !lift(acc[f].denominator, d) || d.is_zero()) {
        reconstructed = false;
        break;
      }
      current.emplace_back(n, d);
    }
    if (!reconstructed) {
      continue;
    }
    if (have_previous && current.size() == previous.size()) {
      bool same = true;
      for (std::size_t f = 0; f < count && same; ++f) {
        same = current[f].numerator() == previous[f].numerator() &&
               current[f].denominator() == previous[f].denominator();
      }
      if (same) {
        out = std::move(current);
        return true;
      }
    }
    previous = std::move(current);
    have_previous = true;
  }
  return false;
}

}  // namespace reconstruction
}  // namespace varietas

#endif
