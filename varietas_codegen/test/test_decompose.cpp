// The decomposition along factors, on varieties whose pieces are known.
//
// Two things are checked separately, as for the factorisation. The pieces have
// to be the expected ones, which the small examples here fix by hand. And
// their union has to be V(I) exactly, which is checked without trusting the
// algorithm: every piece's ideal contains I, so each piece lies inside V(I),
// and a product taking one basis element from each piece lies in the radical
// of I, so every point of V(I) lies in some piece. Membership in the radical
// is decided by Rabinowitsch's trick, f in sqrt(I) exactly when I : f^inf is
// the unit ideal, which is saturation and has nothing to do with factoring.

#include <algorithm>
#include <cstddef>
#include <vector>

#include <gtest/gtest.h>

#include "varietas/codegen/decompose.hpp"
#include "varietas/codegen/rational.hpp"
#include "varietas/core/ideal/buchberger.hpp"
#include "varietas/core/ideal/saturation.hpp"
#include "varietas/core/order/grevlex.hpp"
#include "varietas/core/order/lex.hpp"
#include "varietas/core/polynomial.hpp"
#include "varietas/core/quotient/quotient_basis.hpp"

namespace {

using varietas::grevlex;
using varietas::make_rational;
using varietas::rational;

using poly = varietas::polynomial<rational, 3, grevlex>;

poly x() { return poly::variable(0); }
poly y() { return poly::variable(1); }
poly z() { return poly::variable(2); }
poly c(std::int64_t n, std::int64_t d = 1) { return poly::constant(make_rational(n, d)); }

bool in_radical(const poly& f, const std::vector<poly>& generators) {
  return varietas::is_unit_ideal(varietas::saturate(generators, f));
}

template <class Pieces>
void expect_exact_union(const std::vector<poly>& generators, const Pieces& pieces) {
  const auto basis = varietas::groebner_basis(generators);
  // Each piece inside V(I).
  for (const auto& piece : pieces) {
    for (const poly& g : basis) {
      EXPECT_TRUE(varietas::is_member(g, piece.basis));
    }
  }
  // V(I) inside the union: every product of one element from each piece
  // vanishes on V(I).
  std::vector<std::size_t> index(pieces.size(), 0);
  while (true) {
    poly product = c(1);
    for (std::size_t k = 0; k < pieces.size(); ++k) {
      product = product * pieces[k].basis[index[k]];
    }
    EXPECT_TRUE(in_radical(product, generators));
    std::size_t k = 0;
    while (k < pieces.size() && ++index[k] == pieces[k].basis.size()) {
      index[k] = 0;
      ++k;
    }
    if (k == pieces.size()) {
      break;
    }
  }
}

// Whether some piece has exactly this ideal.
template <class Pieces>
bool has_piece(const Pieces& pieces, const std::vector<poly>& generators) {
  const auto basis = varietas::groebner_basis(generators);
  return std::any_of(pieces.begin(), pieces.end(),
                     [&](const auto& piece) { return piece.basis == basis; });
}

TEST(Decompose, ALineAndACircle) {
  const std::vector<poly> generators{(x() - y()) * (x() * x() + y() * y() - c(1)), z()};
  const auto pieces = varietas::decompose(generators);
  ASSERT_EQ(pieces.size(), 2u);
  EXPECT_TRUE(has_piece(pieces, {x() - y(), z()}));
  EXPECT_TRUE(has_piece(pieces, {x() * x() + y() * y() - c(1), z()}));
  for (const auto& piece : pieces) {
    EXPECT_EQ(piece.dimension.dimension, 1u);
  }
  expect_exact_union(generators, pieces);
}

// x y = x z = 0 is the plane x = 0 together with the line y = z = 0, and no
// generator says so by itself: x y factors into x and y, and the branch on y
// still has x z, which factors again.
TEST(Decompose, APlaneAndALineOfDifferentDimension) {
  const std::vector<poly> generators{x() * y(), x() * z()};
  const auto pieces = varietas::decompose(generators);
  ASSERT_EQ(pieces.size(), 2u);
  EXPECT_EQ(pieces[0].basis, varietas::groebner_basis(std::vector<poly>{x()}));
  EXPECT_EQ(pieces[0].dimension.dimension, 2u);
  EXPECT_EQ(pieces[1].basis, varietas::groebner_basis(std::vector<poly>{y(), z()}));
  EXPECT_EQ(pieces[1].dimension.dimension, 1u);
  expect_exact_union(generators, pieces);
}

// x (x - 1) = 0 splits into the line x = 0, where x y vanishes already, and
// the plane x = 1, where x y = 0 leaves only the point y = 0.
TEST(Decompose, LaterBranchesStayOffEarlierOnes) {
  const std::vector<poly> generators{x() * (x() - c(1)), x() * y(), z()};
  const auto pieces = varietas::decompose(generators);
  ASSERT_EQ(pieces.size(), 2u);
  EXPECT_TRUE(has_piece(pieces, {x(), z()}));
  EXPECT_TRUE(has_piece(pieces, {x() - c(1), y(), z()}));
  expect_exact_union(generators, pieces);
}

TEST(Decompose, AnEmptyVarietyHasNoPieces) {
  EXPECT_TRUE(varietas::decompose(std::vector<poly>{x() * y(), x() - c(1), y() - c(1)}).empty());
}

// The decomposition is of the variety: z^2 = 0 is the plane z = 0, and the
// square the ideal recorded is not carried into the piece.
TEST(Decompose, RepeatedFactorsDescribeTheSamePoints) {
  const std::vector<poly> generators{z() * z() * (x() - c(2))};
  const auto pieces = varietas::decompose(generators);
  ASSERT_EQ(pieces.size(), 2u);
  EXPECT_TRUE(has_piece(pieces, {z()}));
  EXPECT_TRUE(has_piece(pieces, {x() - c(2)}));
}

// y^2 - 2 x^2 is two planes over Q(sqrt 2) and one over Q, where the
// factorisation happens. The decomposition says one, which is the limitation
// stated, not a bug.
TEST(Decompose, AFactorisationOverAnExtensionIsNotFound) {
  const std::vector<poly> generators{y() * y() - c(2) * x() * x()};
  const auto pieces = varietas::decompose(generators);
  ASSERT_EQ(pieces.size(), 1u);
  EXPECT_EQ(pieces[0].basis, varietas::groebner_basis(generators));
}

// Nothing factors: the decomposition is the ideal itself.
TEST(Decompose, AnIrreducibleVarietyIsOnePiece) {
  const std::vector<poly> generators{x() * x() + y() * y() + z() * z() - c(1), x() + y() + z()};
  const auto pieces = varietas::decompose(generators);
  ASSERT_EQ(pieces.size(), 1u);
  EXPECT_EQ(pieces[0].basis, varietas::groebner_basis(generators));
  EXPECT_EQ(pieces[0].dimension.dimension, 1u);
}

// Points of the twisted cubic (t, t^2, t^3) where z = 1, z = 8 or z^2 = 27,
// none of which is one point. z = 1 is t^3 = 1, which is t = 1 and the two
// roots of t^2 + t + 1; z = 8 is t = 2 and the two roots of t^2 + 2t + 4. On
// z^2 = 27, y^3 = z^2 = 27 by the cubic's own equations, and y^3 - 27 is
// (y - 3)(y^2 + 3y + 9): two points where y = 3 and four where y is a complex
// cube root of 27. Those factors are in no generator; they appear only in the
// bases of the branches.
TEST(Decompose, PointsOnATwistedCubic) {
  using lpoly = varietas::polynomial<rational, 3, varietas::lex>;
  const lpoly X = lpoly::variable(0);
  const lpoly Y = lpoly::variable(1);
  const lpoly Z = lpoly::variable(2);
  const auto k = [](std::int64_t n) { return lpoly::constant(rational(n)); };
  const std::vector<lpoly> generators{
      X * X - Y, X * Y - Z, (Z - k(1)) * (Z - k(8)) * (Z * Z - k(27))};
  varietas::decomposition_statistics statistics;
  const auto pieces = varietas::decompose(generators, &statistics);
  ASSERT_EQ(pieces.size(), 6u);
  std::vector<std::size_t> points;
  for (const auto& piece : pieces) {
    EXPECT_EQ(piece.dimension.dimension, 0u);
    // The number of points, with multiplicity, is the dimension of the quotient.
    points.push_back(varietas::standard_monomials(piece.basis).dimension());
  }
  std::sort(points.begin(), points.end());
  EXPECT_EQ(points, (std::vector<std::size_t>{1, 1, 2, 2, 2, 4}));
  EXPECT_GE(statistics.factorisations, 4u);
}

}  // namespace
