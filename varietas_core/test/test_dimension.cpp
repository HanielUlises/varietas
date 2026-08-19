// Dimension of a variety, and determinants over a ring that is not a field.
//
// The two are tested together because they arrive together: the singular locus
// is cut out by minors of a polynomial matrix, and the only useful thing to say
// about the set that results is how big it is. The cases below are chosen so
// that the answer is known independently of any Gröbner computation — a point,
// a line, a plane curve, a hypersurface — and so that the two verdicts the
// library now offers about the size of a variety, the finiteness criterion and
// the dimension, are checked against each other where they overlap.

#include <cstdint>
#include <vector>

#include <gtest/gtest.h>

#include "varietas/core/ideal/dimension.hpp"
#include "varietas/core/ideal/ideal.hpp"
#include "varietas/core/minors.hpp"
#include "varietas/core/order/grevlex.hpp"
#include "varietas/core/polynomial.hpp"

namespace {

using varietas::grevlex;

template <std::size_t N>
using poly = varietas::polynomial<double, N, grevlex>;

template <std::size_t N>
poly<N> var(std::size_t i, std::uint8_t power = 1) {
  return poly<N>::variable(i, power);
}

template <std::size_t N>
poly<N> constant(double c) {
  return poly<N>::constant(c);
}

template <std::size_t N>
varietas::affine_dimension<N> dimension_of(const std::vector<poly<N>>& generators) {
  return varietas::ideal<double, N, grevlex>(generators).dimension();
}

// --- dimension -------------------------------------------------------------

TEST(dimension, the_zero_ideal_is_the_whole_affine_space) {
  const auto d = dimension_of<3>({});
  EXPECT_FALSE(d.is_empty);
  EXPECT_EQ(d.dimension, 3u);
  for (std::size_t i = 0; i < 3; ++i) {
    EXPECT_TRUE(d.independent[i]);
  }
}

TEST(dimension, the_unit_ideal_is_empty_rather_than_zero_dimensional) {
  // x and x - 1 have no common zero, so the ideal is the whole ring.
  const auto d = dimension_of<2>({var<2>(0), var<2>(0) - constant<2>(1.0)});
  EXPECT_TRUE(d.is_empty);
  EXPECT_FALSE(d.is_zero_dimensional())
      << "an empty variety is not a finite nonempty one, and the distinction is "
         "the whole reason is_empty is reported";
}

TEST(dimension, a_point_is_zero_dimensional_and_agrees_with_the_finiteness_criterion) {
  const varietas::ideal<double, 3, grevlex> point(
      {var<3>(0) - constant<3>(1.0), var<3>(1) - constant<3>(2.0),
       var<3>(2) - constant<3>(3.0)});

  const auto d = point.dimension();
  EXPECT_FALSE(d.is_empty);
  EXPECT_EQ(d.dimension, 0u);
  EXPECT_TRUE(d.is_zero_dimensional());
  EXPECT_TRUE(point.is_zero_dimensional()) << "the two verdicts must agree at dimension zero";
}

TEST(dimension, a_line_in_three_space_is_one_dimensional_and_names_its_free_variable) {
  // x = 0, y = 0: the z axis. The ideal leaves z unconstrained and nothing else.
  const auto d = dimension_of<3>({var<3>(0), var<3>(1)});
  EXPECT_FALSE(d.is_empty);
  EXPECT_EQ(d.dimension, 1u);
  EXPECT_FALSE(d.independent[0]);
  EXPECT_FALSE(d.independent[1]);
  EXPECT_TRUE(d.independent[2]) << "the independent set must be the coordinate the line moves in";
}

TEST(dimension, a_hypersurface_has_dimension_one_less_than_its_space) {
  // The unit sphere in three variables.
  const auto d = dimension_of<3>(
      {var<3>(0, 2) + var<3>(1, 2) + var<3>(2, 2) - constant<3>(1.0)});
  EXPECT_EQ(d.dimension, 2u);
  EXPECT_FALSE(d.is_zero_dimensional());
}

TEST(dimension, a_curve_cut_out_by_two_surfaces_in_three_space) {
  // The circle x^2 + y^2 = 1 in the plane z = 0.
  const auto d = dimension_of<3>({var<3>(2), var<3>(0, 2) + var<3>(1, 2) - constant<3>(1.0)});
  EXPECT_EQ(d.dimension, 1u);
}

TEST(dimension, is_read_from_leading_terms_and_so_ignores_embedded_structure) {
  // The union of a plane and a line meeting it: x·z = 0, y·z = 0 is the plane
  // z = 0 together with the z axis. Dimension is the maximum over components,
  // which is two, and the report says nothing about the line — a dimension is
  // one number, and separating the pieces is what splitting is for.
  const auto d = dimension_of<3>({var<3>(0) * var<3>(2), var<3>(1) * var<3>(2)});
  EXPECT_EQ(d.dimension, 2u);
}

// --- determinants over a polynomial ring -----------------------------------

TEST(minors, expands_a_determinant_over_a_ring_without_dividing) {
  // [[x, y], [z, x]] has determinant x^2 - y z.
  varietas::dense_matrix<poly<3>> m(2, 2);
  m(0, 0) = var<3>(0);
  m(0, 1) = var<3>(1);
  m(1, 0) = var<3>(2);
  m(1, 1) = var<3>(0);

  EXPECT_EQ(varietas::determinant(m), var<3>(0, 2) - var<3>(1) * var<3>(2));
}

TEST(minors, a_three_by_three_determinant_matches_the_cofactor_expansion) {
  // A matrix with a known determinant: the circulant on x, y, z with rows
  // (x, y, z), (y, z, x), (z, x, y). That is one row swap away from the
  // standard circulant, so its determinant is -(x^3 + y^3 + z^3 - 3xyz).
  varietas::dense_matrix<poly<3>> m(3, 3);
  const std::array<std::size_t, 3> row0 = {0, 1, 2};
  for (std::size_t i = 0; i < 3; ++i) {
    for (std::size_t j = 0; j < 3; ++j) {
      m(i, j) = var<3>(row0[(j + i) % 3]);
    }
  }

  const poly<3> expected = constant<3>(3.0) * var<3>(0) * var<3>(1) * var<3>(2) -
                           var<3>(0, 3) - var<3>(1, 3) - var<3>(2, 3);
  EXPECT_EQ(varietas::determinant(m), expected);
}

TEST(minors, of_a_wide_matrix_are_the_maximal_ones_and_cut_out_rank_deficiency) {
  // [[x, y, z], [1, 1, 1]] drops rank exactly where its three 2x2 minors vanish,
  // that is where x = y = z.
  varietas::dense_matrix<poly<3>> m(2, 3);
  for (std::size_t j = 0; j < 3; ++j) {
    m(0, j) = var<3>(j);
    m(1, j) = constant<3>(1.0);
  }

  const auto minors = varietas::maximal_minors(m);
  ASSERT_EQ(minors.size(), 3u);

  const varietas::ideal<double, 3, grevlex> rank_deficient(minors);
  EXPECT_TRUE(rank_deficient.contains(var<3>(0) - var<3>(1)));
  EXPECT_TRUE(rank_deficient.contains(var<3>(1) - var<3>(2)));
  EXPECT_FALSE(rank_deficient.contains(var<3>(0)));
  EXPECT_EQ(rank_deficient.dimension().dimension, 1u)
      << "the diagonal line is one dimensional";
}

TEST(minors, a_zero_minor_is_dropped_rather_than_returned) {
  // A matrix with a repeated row: every maximal minor vanishes identically.
  varietas::dense_matrix<poly<2>> m(2, 2);
  m(0, 0) = var<2>(0);
  m(0, 1) = var<2>(1);
  m(1, 0) = var<2>(0);
  m(1, 1) = var<2>(1);

  EXPECT_TRUE(varietas::maximal_minors(m).empty());
}

}  // namespace
