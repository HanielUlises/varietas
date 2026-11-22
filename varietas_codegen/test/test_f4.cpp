// F4 against Buchberger.
//
// The reduced Gröbner basis of an ideal is unique for a fixed order, so the
// two algorithms have to return the same polynomials, coefficient for
// coefficient, and that is the whole of what is checked here: on random
// systems that are zero-dimensional, positive-dimensional and inconsistent,
// under grevlex, lex and the block order saturation uses, and on the
// saturated ideal of an arm at a pose, which is what the reconstruction asks
// for hundreds of times per prime.

#include <array>
#include <cstdint>
#include <random>
#include <vector>

#include <gtest/gtest.h>

#include "varietas/codegen/f4.hpp"
#include "varietas/codegen/prime_field.hpp"
#include "varietas/core/ideal/buchberger.hpp"
#include "varietas/core/ideal/saturation.hpp"
#include "varietas/core/order/block.hpp"
#include "varietas/core/order/grevlex.hpp"
#include "varietas/core/order/lex.hpp"

namespace {

using varietas::block_order;
using varietas::grevlex;
using varietas::lex;
using varietas::monomial;
using varietas::polynomial;
using varietas::residue;

constexpr residue::word prime = 2147483629;

template <std::size_t N, class Order>
std::vector<polynomial<residue, N, Order>> by_buchberger(
    const std::vector<polynomial<residue, N, Order>>& generators) {
  auto basis = varietas::buchberger(generators);
  varietas::reduce_groebner_basis(basis);
  return basis;
}

// `count` polynomials of `terms` terms each, every exponent vector of total
// degree at most `degree`, every coefficient nonzero.
template <std::size_t N, class Order>
std::vector<polynomial<residue, N, Order>> random_system(std::mt19937_64& rng, int count,
                                                         int degree, int terms) {
  using poly = polynomial<residue, N, Order>;
  std::uniform_int_distribution<residue::word> nonzero(1, prime - 1);
  std::vector<poly> out;
  for (int k = 0; k < count; ++k) {
    std::vector<typename poly::term> ts;
    for (int t = 0; t < terms; ++t) {
      std::array<std::uint16_t, N> e{};
      int left = std::uniform_int_distribution<int>(0, degree)(rng);
      for (std::size_t v = 0; v < N && left > 0; ++v) {
        const int x = std::uniform_int_distribution<int>(0, left)(rng);
        e[v] = static_cast<std::uint16_t>(x);
        left -= x;
      }
      ts.push_back({monomial<N>(e), residue::raw(nonzero(rng))});
    }
    out.emplace_back(std::move(ts));
  }
  return out;
}

template <std::size_t N, class Order>
void agree_on_random_systems(int trials, int count, int degree, int terms, std::uint64_t seed) {
  residue::select(prime);
  std::mt19937_64 rng(seed);
  for (int k = 0; k < trials; ++k) {
    const auto generators = random_system<N, Order>(rng, count, degree, terms);
    ASSERT_EQ(varietas::groebner_basis(generators), by_buchberger(generators)) << "trial " << k;
  }
}

TEST(F4, AgreesWithBuchbergerOnZeroDimensionalSystems) {
  agree_on_random_systems<3, grevlex>(100, 3, 3, 4, 1);
  agree_on_random_systems<4, grevlex>(40, 4, 2, 4, 2);
}

TEST(F4, AgreesWithBuchbergerUnderLex) { agree_on_random_systems<3, lex>(60, 3, 2, 3, 3); }

TEST(F4, AgreesWithBuchbergerOnPositiveDimensionalSystems) {
  agree_on_random_systems<3, grevlex>(60, 2, 3, 3, 4);
}

TEST(F4, AgreesWithBuchbergerUnderTheSaturationOrder) {
  agree_on_random_systems<4, block_order<1, lex, grevlex>>(60, 3, 3, 3, 5);
}

TEST(F4, FindsTheUnitIdeal) {
  residue::select(prime);
  using poly = polynomial<residue, 2, grevlex>;
  const poly x = poly::variable(0);
  const poly y = poly::variable(1);
  const poly one = poly::constant(residue::raw(1));
  // x y = 1 and x = 0 have no common solution.
  const auto basis = varietas::groebner_basis(std::vector<poly>{x * y - one, x});
  ASSERT_EQ(basis.size(), 1u);
  EXPECT_TRUE(basis.front().leading_monomial().is_one());
  EXPECT_TRUE(varietas::is_unit_ideal(basis));
}

TEST(F4, LeavesTheZeroIdealEmpty) {
  residue::select(prime);
  using poly = polynomial<residue, 2, grevlex>;
  EXPECT_TRUE(varietas::groebner_basis(std::vector<poly>{}).empty());
  EXPECT_TRUE(varietas::groebner_basis(std::vector<poly>{poly()}).empty());
}

// The saturation of a 3R arm's position residuals by the half-angle
// denominators, at random poses over Z/p. The arm has three pairwise skew axes,
// the one the reconstruction spends the most samples on.
TEST(F4, AgreesWithBuchbergerOnTheSaturatedIdealOfAnArm) {
  residue::select(prime);
  using poly = polynomial<residue, 3, grevlex>;
  const poly t1 = poly::variable(0);
  const poly t2 = poly::variable(1);
  const poly t3 = poly::variable(2);
  const poly one = poly::constant(residue::raw(1));
  const auto c = [](long long n, long long d = 1) { return poly::constant(residue(n) / residue(d)); };

  // cos and sin of each joint, times 1 + t^2.
  const poly c1 = one - t1 * t1, s1 = c(2) * t1, d1 = one + t1 * t1;
  const poly c2 = one - t2 * t2, s2 = c(2) * t2, d2 = one + t2 * t2;
  const poly c3 = one - t3 * t3, s3 = c(2) * t3, d3 = one + t3 * t3;
  // A position map of the shape the rationalisation produces: polynomial in
  // the cosines and sines, cleared by the product of the denominators. The
  // coefficients are those of no arm in particular; what matters is that no
  // structure makes the ideal special.
  const poly x = c1 * c2 * c3 * c(3, 5) + s1 * s2 * d3 * c(1, 2) - c1 * s3 * d2 * c(2, 7) +
                 d1 * d2 * d3 * c(1, 3) + s1 * c2 * s3 * c(4, 9);
  const poly y = s1 * c2 * c3 * c(2, 5) - c1 * s2 * d3 * c(1, 3) + s1 * s3 * d2 * c(5, 7) +
                 c1 * c2 * s3 * c(1, 6) - d1 * s2 * c3 * c(3, 4);
  const poly z = s2 * c3 * d1 * c(4, 5) + c2 * s3 * d1 * c(3, 7) - d1 * d2 * d3 * c(1, 5) +
                 s1 * s2 * s3 * c(2, 3);
  const poly denominator = d1 * d2 * d3;

  std::mt19937_64 rng(6);
  std::uniform_int_distribution<residue::word> draw(1, prime - 1);
  using sat = block_order<1, lex, grevlex>;
  using lifted = polynomial<residue, 4, sat>;
  for (int k = 0; k < 10; ++k) {
    std::vector<poly> residuals{x - denominator * poly::constant(residue::raw(draw(rng))),
                                y - denominator * poly::constant(residue::raw(draw(rng))),
                                z - denominator * poly::constant(residue::raw(draw(rng)))};
    std::vector<lifted> generators;
    for (const poly& r : residuals) {
      generators.push_back(varietas::embed_polynomial<sat, 4>(r, 1));
    }
    generators.push_back(lifted::constant(residue::raw(1)) -
                         lifted::variable(0) * varietas::embed_polynomial<sat, 4>(denominator, 1));
    ASSERT_EQ(varietas::groebner_basis(generators), by_buchberger(generators)) << "pose " << k;
  }
}

// The statistics are filled in the same terms as Buchberger's, so that a
// caller reading them need not know which engine ran.
TEST(F4, ReportsItsPairs) {
  residue::select(prime);
  std::mt19937_64 rng(8);
  const auto generators = random_system<3, grevlex>(rng, 3, 3, 4);
  varietas::buchberger_statistics s;
  varietas::groebner_basis(generators, &s);
  EXPECT_GT(s.pairs_generated, 0u);
  EXPECT_LE(s.pairs_discarded_coprime + s.pairs_discarded_chain, s.pairs_generated);
  EXPECT_LE(s.reductions_to_zero, s.pairs_reduced);
}

}  // namespace
