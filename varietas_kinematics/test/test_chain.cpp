#include <cmath>
#include <cstdint>
#include <string>

#include <gtest/gtest.h>

#include "varietas/codegen/rational.hpp"
#include "varietas/kinematics/chain.hpp"
#include "varietas/kinematics/rigid_transform.hpp"

namespace {

using varietas::chain;
using varietas::chain_status;
using varietas::fixed_joint;
using varietas::joint;
using varietas::joint_type;
using varietas::matrix3;
using varietas::prismatic_joint;
using varietas::revolute_joint;
using varietas::rigid_transform;
using varietas::rotation_about_axis;
using varietas::rotation_from_quaternion;
using varietas::vector3;

// Builds a coefficient from an integer fraction in whichever field is under
// test, mirroring varietas_codegen's fixture so that one body runs over both.
template <class Coeff>
struct factory;

template <>
struct factory<double> {
  static double of(std::int64_t n, std::int64_t d = 1) {
    return static_cast<double>(n) / static_cast<double>(d);
  }
  static double of_double(double x) { return x; }
};

template <>
struct factory<varietas::rational> {
  static varietas::rational of(std::int64_t n, std::int64_t d = 1) {
    return varietas::make_rational(n, d);
  }
  // The binary value actually stored, promoted to the rational that equals it.
  static varietas::rational of_double(double x) { return varietas::rational(x); }
};

template <class Coeff>
Coeff c(std::int64_t n, std::int64_t d = 1) {
  return factory<Coeff>::of(n, d);
}

template <class Coeff>
vector3<Coeff> translation(std::int64_t x, std::int64_t y, std::int64_t z) {
  return vector3<Coeff>(c<Coeff>(x), c<Coeff>(y), c<Coeff>(z));
}

// Agreement in whatever sense the field supports: identical over the exact
// field, and to within rounding over double, where composing two transforms and
// applying them in turn are two different summation orders and need not give
// the same bits.
template <class Coeff>
void expect_agrees(const vector3<Coeff>& a, const vector3<Coeff>& b) {
  using traits = varietas::coefficient_traits<Coeff>;
  if (traits::is_exact) {
    EXPECT_EQ(a, b);
    return;
  }
  for (std::size_t i = 0; i < 3; ++i) {
    EXPECT_NEAR(traits::to_double(a[i]), traits::to_double(b[i]), 1e-12);
  }
}

// A planar chain of unit-length links rotating about z, the running example of
// the library and the case whose inverse kinematics is known in closed form.
template <class Coeff>
chain<Coeff> planar_two_link() {
  chain<Coeff> robot("planar_2r");
  robot.add_joint(revolute_joint<Coeff>("q1", vector3<Coeff>::unit(2),
                                        rigid_transform<Coeff>::identity()));
  robot.add_joint(revolute_joint<Coeff>(
      "q2", vector3<Coeff>::unit(2),
      rigid_transform<Coeff>::translation_only(translation<Coeff>(1, 0, 0))));
  robot.set_tool(
      rigid_transform<Coeff>::translation_only(translation<Coeff>(1, 0, 0)));
  return robot;
}

template <class Coeff>
class chain_test : public ::testing::Test {};

using coefficient_types = ::testing::Types<double, varietas::rational>;
TYPED_TEST_SUITE(chain_test, coefficient_types);

// --- rigid transforms ------------------------------------------------------

TYPED_TEST(chain_test, quaternion_rotation_is_orthogonal_without_normalisation) {
  using Coeff = TypeParam;

  // (1, 1, 0, 0) is the quarter turn about x, of norm sqrt(2): the Rodrigues
  // matrix is rational even though the quaternion is not of unit length.
  const matrix3<Coeff> r = rotation_from_quaternion<Coeff>(
      c<Coeff>(1), c<Coeff>(1), c<Coeff>(0), c<Coeff>(0));

  EXPECT_DOUBLE_EQ(orthogonality_defect(r), 0.0);
  EXPECT_EQ(r, (matrix3<Coeff>::from_rows({c<Coeff>(1), c<Coeff>(0), c<Coeff>(0),
                                           c<Coeff>(0), c<Coeff>(0), c<Coeff>(-1),
                                           c<Coeff>(0), c<Coeff>(1), c<Coeff>(0)})));
}

TYPED_TEST(chain_test, quaternion_rotation_has_unit_determinant) {
  using Coeff = TypeParam;

  // A quaternion with no special structure, to check the general formula.
  const matrix3<Coeff> r = rotation_from_quaternion<Coeff>(
      c<Coeff>(2), c<Coeff>(-1), c<Coeff>(3), c<Coeff>(5));

  EXPECT_DOUBLE_EQ(orthogonality_defect(r), 0.0);
  EXPECT_NEAR(varietas::coefficient_traits<Coeff>::to_double(r.determinant()), 1.0,
              1e-12);
}

TYPED_TEST(chain_test, rodrigues_quarter_turn_about_z) {
  using Coeff = TypeParam;

  const matrix3<Coeff> r = rotation_about_axis<Coeff>(vector3<Coeff>::unit(2),
                                                      c<Coeff>(0), c<Coeff>(1));

  const vector3<Coeff> x = vector3<Coeff>::unit(0);
  EXPECT_EQ(r * x, vector3<Coeff>::unit(1));
  EXPECT_DOUBLE_EQ(orthogonality_defect(r), 0.0);
}

TYPED_TEST(chain_test, composition_and_inverse_close_in_se3) {
  using Coeff = TypeParam;

  const rigid_transform<Coeff> a(
      rotation_from_quaternion<Coeff>(c<Coeff>(1), c<Coeff>(0), c<Coeff>(0),
                                      c<Coeff>(1)),
      translation<Coeff>(2, -3, 5));
  const rigid_transform<Coeff> b(
      rotation_from_quaternion<Coeff>(c<Coeff>(3), c<Coeff>(1), c<Coeff>(-2),
                                      c<Coeff>(0)),
      translation<Coeff>(-1, 4, 0));

  // The inverse uses the transpose rather than a solve, so it costs no
  // division and is exact in both fields.
  EXPECT_EQ(a * a.inverse(), rigid_transform<Coeff>::identity());

  // Composition agrees with applying the two transforms in turn.
  const vector3<Coeff> p = translation<Coeff>(1, 2, 3);
  expect_agrees<Coeff>((a * b).apply(p), a.apply(b.apply(p)));
}

// --- chain structure -------------------------------------------------------

TYPED_TEST(chain_test, planar_two_link_validates) {
  using Coeff = TypeParam;

  const chain<Coeff> robot = planar_two_link<Coeff>();
  const auto diagnostic = robot.validate();

  EXPECT_TRUE(diagnostic.ok()) << to_string(diagnostic.status);
  EXPECT_EQ(robot.degrees_of_freedom(), 2u);
  EXPECT_EQ(robot.size(), 2u);
}

TYPED_TEST(chain_test, variable_indices_skip_fixed_joints) {
  using Coeff = TypeParam;

  chain<Coeff> robot("mixed");
  robot.add_joint(fixed_joint<Coeff>("base", rigid_transform<Coeff>::identity()));
  robot.add_joint(revolute_joint<Coeff>("q1", vector3<Coeff>::unit(2),
                                        rigid_transform<Coeff>::identity()));
  robot.add_joint(fixed_joint<Coeff>("spacer", rigid_transform<Coeff>::identity()));
  robot.add_joint(prismatic_joint<Coeff>("d2", vector3<Coeff>::unit(0),
                                         rigid_transform<Coeff>::identity()));

  EXPECT_EQ(robot.degrees_of_freedom(), 2u);
  EXPECT_EQ(robot.variable_of_joint(0), chain<Coeff>::npos);
  EXPECT_EQ(robot.variable_of_joint(1), 0u);
  EXPECT_EQ(robot.variable_of_joint(2), chain<Coeff>::npos);
  EXPECT_EQ(robot.variable_of_joint(3), 1u);

  EXPECT_EQ(robot.joint_of_variable(0), 1u);
  EXPECT_EQ(robot.joint_of_variable(1), 3u);
  EXPECT_EQ(robot.joint_of_variable(2), chain<Coeff>::npos);
}

TYPED_TEST(chain_test, folding_fixed_joints_preserves_the_geometry) {
  using Coeff = TypeParam;

  chain<Coeff> robot("with_fixed");
  robot.add_joint(fixed_joint<Coeff>(
      "base", rigid_transform<Coeff>::translation_only(translation<Coeff>(0, 0, 1))));
  robot.add_joint(revolute_joint<Coeff>("q1", vector3<Coeff>::unit(2),
                                        rigid_transform<Coeff>::identity()));
  robot.add_joint(fixed_joint<Coeff>(
      "wrist",
      rigid_transform<Coeff>::translation_only(translation<Coeff>(1, 0, 0))));
  robot.set_tool(
      rigid_transform<Coeff>::translation_only(translation<Coeff>(0, 2, 0)));

  const chain<Coeff> folded = robot.fold_fixed_joints();

  EXPECT_EQ(folded.size(), 1u);
  EXPECT_EQ(folded.degrees_of_freedom(), robot.degrees_of_freedom());

  // The leading fixed joint has moved into the origin of the joint it preceded,
  // and the trailing one into the tool frame.
  EXPECT_EQ(folded.joints()[0].origin.translation(), translation<Coeff>(0, 0, 1));
  EXPECT_EQ(folded.tool(), robot.trailing_transform());
  EXPECT_EQ(folded.tool().translation(), translation<Coeff>(1, 2, 0));
}

TYPED_TEST(chain_test, trailing_transform_stops_at_the_last_moving_joint) {
  using Coeff = TypeParam;

  const chain<Coeff> robot = planar_two_link<Coeff>();
  EXPECT_EQ(robot.trailing_transform(), robot.tool());
}

// --- named failure modes ---------------------------------------------------

TYPED_TEST(chain_test, rejects_a_non_unit_axis) {
  using Coeff = TypeParam;

  chain<Coeff> robot("scaled_axis");
  robot.add_joint(revolute_joint<Coeff>(
      "q1", vector3<Coeff>(c<Coeff>(0), c<Coeff>(0), c<Coeff>(2)),
      rigid_transform<Coeff>::identity()));

  const auto diagnostic = robot.validate();
  EXPECT_FALSE(diagnostic.ok());
  EXPECT_EQ(diagnostic.status, chain_status::axis_not_unit);
  EXPECT_EQ(diagnostic.joint_index, 0u);
}

TYPED_TEST(chain_test, rejects_a_zero_axis) {
  using Coeff = TypeParam;

  chain<Coeff> robot("zero_axis");
  robot.add_joint(revolute_joint<Coeff>("q1", vector3<Coeff>::zero(),
                                        rigid_transform<Coeff>::identity()));

  EXPECT_EQ(robot.validate().status, chain_status::axis_degenerate);
}

TYPED_TEST(chain_test, rejects_a_reflection_as_an_origin) {
  using Coeff = TypeParam;

  // Orthogonal, so the Gram test passes, but of determinant -1.
  const matrix3<Coeff> reflection = matrix3<Coeff>::from_rows(
      {c<Coeff>(1), c<Coeff>(0), c<Coeff>(0), c<Coeff>(0), c<Coeff>(1), c<Coeff>(0),
       c<Coeff>(0), c<Coeff>(0), c<Coeff>(-1)});

  chain<Coeff> robot("reflected");
  robot.add_joint(revolute_joint<Coeff>(
      "q1", vector3<Coeff>::unit(2),
      rigid_transform<Coeff>::rotation_only(reflection)));

  EXPECT_DOUBLE_EQ(orthogonality_defect(reflection), 0.0);
  EXPECT_EQ(robot.validate().status, chain_status::origin_improper);
}

TYPED_TEST(chain_test, rejects_inverted_limits) {
  using Coeff = TypeParam;

  joint<Coeff> j = revolute_joint<Coeff>("q1", vector3<Coeff>::unit(2),
                                         rigid_transform<Coeff>::identity());
  j.has_limits = true;
  j.lower = 1.0;
  j.upper = -1.0;

  chain<Coeff> robot("bad_limits");
  robot.add_joint(std::move(j));

  EXPECT_EQ(robot.validate().status, chain_status::limits_inverted);
}

TYPED_TEST(chain_test, rejects_a_chain_with_nothing_to_solve_for) {
  using Coeff = TypeParam;

  chain<Coeff> robot("all_fixed");
  robot.add_joint(fixed_joint<Coeff>("base", rigid_transform<Coeff>::identity()));

  EXPECT_EQ(robot.validate().status, chain_status::no_actuated_joints);
}

// --- the asymmetry between the two fields ----------------------------------

// A rotation built from the floating point cosine and sine of thirty degrees is
// orthogonal to within rounding and no further. Over double that is accepted,
// because there is no stronger statement available; over the exact field the
// same entries describe a matrix that is provably not a rotation, and it is
// rejected. This is the same distinction test_exact_groebner draws for ideal
// membership, at the point where geometry enters the pipeline: an origin
// accepted here becomes exact input to Buchberger, so a rounded one would make
// the resulting basis an exact statement about a robot that does not exist.
template <class Coeff>
matrix3<Coeff> rounded_thirty_degrees_about_z() {
  const double cosine = std::cos(M_PI / 6.0);
  const double sine = std::sin(M_PI / 6.0);
  return rotation_about_axis<Coeff>(vector3<Coeff>::unit(2),
                                    factory<Coeff>::of_double(cosine),
                                    factory<Coeff>::of_double(sine));
}

TEST(chain_exactness, rounded_rotation_is_accepted_over_double) {
  using Coeff = double;

  const matrix3<Coeff> r = rounded_thirty_degrees_about_z<Coeff>();
  chain<Coeff> robot("rounded");
  robot.add_joint(revolute_joint<Coeff>(
      "q1", vector3<Coeff>::unit(2), rigid_transform<Coeff>::rotation_only(r)));

  EXPECT_TRUE(robot.validate().ok());
  EXPECT_LT(orthogonality_defect(r), 1e-12);
}

TEST(chain_exactness, rounded_rotation_is_rejected_over_the_rational_field) {
  using Coeff = varietas::rational;

  const matrix3<Coeff> r = rounded_thirty_degrees_about_z<Coeff>();
  chain<Coeff> robot("rounded");
  robot.add_joint(revolute_joint<Coeff>(
      "q1", vector3<Coeff>::unit(2), rigid_transform<Coeff>::rotation_only(r)));

  const auto diagnostic = robot.validate();
  EXPECT_FALSE(diagnostic.ok());
  EXPECT_EQ(diagnostic.status, chain_status::origin_not_orthogonal);

  // The defect is of the order of the rounding error and is nonetheless fatal:
  // c^2 + s^2 differs from one by a rational of magnitude about 1e-17.
  EXPECT_GT(diagnostic.defect, 0.0);
  EXPECT_LT(diagnostic.defect, 1e-15);

  // The tolerance argument is not a way round it.
  EXPECT_FALSE(robot.validate(1e-3).ok());
}

// The exact quarter turns and half turns that a URDF is overwhelmingly made of
// pass the same gate, so the strictness costs nothing on real models.
TEST(chain_exactness, right_angle_rotations_are_exact_over_the_rational_field) {
  using Coeff = varietas::rational;

  const Coeff zero = c<Coeff>(0);
  const Coeff one = c<Coeff>(1);
  const Coeff minus_one = c<Coeff>(-1);

  for (const auto& cs : {std::make_pair(one, zero), std::make_pair(zero, one),
                         std::make_pair(minus_one, zero),
                         std::make_pair(zero, minus_one)}) {
    for (std::size_t axis = 0; axis < 3; ++axis) {
      const matrix3<Coeff> r =
          rotation_about_axis<Coeff>(vector3<Coeff>::unit(axis), cs.first, cs.second);
      EXPECT_DOUBLE_EQ(orthogonality_defect(r), 0.0);
      EXPECT_EQ(r.determinant(), one);
    }
  }
}

}  // namespace
