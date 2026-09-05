#ifndef VARIETAS_KINEMATICS_TRIGONOMETRIC_HPP
#define VARIETAS_KINEMATICS_TRIGONOMETRIC_HPP

#include <array>
#include <cmath>
#include <cstddef>
#include <vector>

#include "varietas/core/config.hpp"
#include "varietas/core/polynomial.hpp"
#include "varietas/kinematics/chain.hpp"
#include "varietas/kinematics/rigid_transform.hpp"

namespace varietas {

// The second rationalisation: cosine and sine as independent variables, tied by
// the circle relation.
//
// The tangent half-angle substitution in rationalize.hpp buys one variable per
// joint and pays for it with a denominator. Writing c_i = cos q_i and
// s_i = sin q_i instead makes the trade the other way: two variables per joint,
// no denominators anywhere, and the trigonometric identity carried explicitly
// as a generator,
//
//     c_i^2 + s_i^2 - 1.
//
// Rodrigues' formula is then already polynomial,
//
//     R = I + s [u] + (1 - c) [u]^2,
//
// so a chain composes as an ordinary product of polynomial matrices. Nothing
// has to be cleared, so nothing spurious is attached, so there is nothing to
// saturate away: the loci t_i = ±i that the half-angle substitution invents
// have no counterpart here, because no denominator was ever multiplied through.
// The circle relations are genuine equations of the problem rather than damage
// repair, and V(c^2 + s^2 - 1) is exactly the parameter space of a revolute
// joint including the configuration q = π that the half-angle map sends to
// infinity.
//
// Which formulation is better depends on the question, and the two are kept
// side by side rather than one replacing the other.
//
//   - Inverse kinematics is zero dimensional, and there the variable count
//     dominates: n variables against 2n, and the saturation is a single extra
//     elimination paid once. The half-angle form wins, and it is what the
//     solver and the emitter are built on.
//
//   - Implicitization eliminates the joint variables entirely, and there the
//     denominators dominate: the Rabinowitsch variable joins the block being
//     eliminated, and the cleared residuals carry the denominator's degree into
//     every S-polynomial. The trigonometric form wins, and by a margin large
//     enough to change what is computable, measured on the torus of
//     test_workspace, seventy seconds against ten milliseconds for the
//     identical quartic.
//
// So this header is not a replacement for rationalize.hpp. It is the
// formulation the elimination-shaped constructions use, workspace
// implicitization now and the singular locus next, and rationalize.hpp remains
// the one the fibre-shaped constructions use.

// The number of ring variables a chain needs in this formulation: two for each
// revolute joint, one for each prismatic joint, none for a fixed joint. It is a
// runtime quantity because the joint types are, and the ring size is a template
// parameter, so callers assert the two agree exactly as they do for the
// half-angle ring.
template <class Coeff>
std::size_t trigonometric_variable_count(const chain<Coeff>& robot) {
  std::size_t count = 0;
  for (const joint<Coeff>& j : robot.joints()) {
    switch (j.type) {
      case joint_type::revolute:
        count += 2;
        break;
      case joint_type::prismatic:
        count += 1;
        break;
      case joint_type::fixed:
        break;
    }
  }
  return count;
}

// The index of the cosine variable of each revolute joint, in chain order. The
// sine is the next index. Prismatic joints occupy one index and appear in no
// entry, since they satisfy no circle relation.
template <class Coeff>
std::vector<std::size_t> cosine_variable_indices(const chain<Coeff>& robot) {
  std::vector<std::size_t> indices;
  std::size_t variable = 0;
  for (const joint<Coeff>& j : robot.joints()) {
    switch (j.type) {
      case joint_type::revolute:
        indices.push_back(variable);
        variable += 2;
        break;
      case joint_type::prismatic:
        variable += 1;
        break;
      case joint_type::fixed:
        break;
    }
  }
  return indices;
}

// A rigid transform whose entries are polynomials in the joint variables. The
// contrast with rational_transform is the whole point of the header: there is
// no denominator field, because there is no denominator.
template <class Coeff, std::size_t N, class Order>
class trigonometric_transform {
 public:
  using polynomial_type = polynomial<Coeff, N, Order>;
  using traits = coefficient_traits<Coeff>;

  static constexpr std::size_t num_vars = N;

  trigonometric_transform() {
    for (std::size_t i = 0; i < 3; ++i) {
      rotation_[3 * i + i] = polynomial_type::constant(traits::one());
    }
  }

  static trigonometric_transform identity() { return trigonometric_transform(); }

  static trigonometric_transform constant(const rigid_transform<Coeff>& t) {
    trigonometric_transform r;
    for (std::size_t i = 0; i < 3; ++i) {
      for (std::size_t j = 0; j < 3; ++j) {
        r.rotation_[3 * i + j] = polynomial_type::constant(t.rotation()(i, j));
      }
      r.translation_[i] = polynomial_type::constant(t.translation()[i]);
    }
    return r;
  }

  // Rodrigues' formula about the unit axis u, with the cosine in variable k and
  // the sine in variable k + 1. As in the half-angle form, [u]^2 = u u^T - I
  // holds because u is a unit vector, which chain::validate has checked.
  static trigonometric_transform revolute(const vector3<Coeff>& u, std::size_t k) {
    VARIETAS_ASSERT(k + 1 < N);

    const polynomial_type c = polynomial_type::variable(k);
    const polynomial_type s = polynomial_type::variable(k + 1);
    const polynomial_type one = polynomial_type::constant(traits::one());
    const polynomial_type one_minus_c = one - c;

    trigonometric_transform r;
    for (std::size_t i = 0; i < 3; ++i) {
      for (std::size_t j = 0; j < 3; ++j) {
        Coeff outer = u[i] * u[j];
        if (i == j) {
          outer = outer - traits::one();
        }
        const Coeff skew = skew_entry(u, i, j);

        polynomial_type entry = polynomial_type::constant(skew) * s +
                                polynomial_type::constant(outer) * one_minus_c;
        if (i == j) {
          entry = entry + one;
        }
        r.rotation_[3 * i + j] = entry;
      }
    }
    return r;
  }

  // A prismatic joint is already polynomial in its displacement and is
  // unaffected by the choice of rationalisation.
  static trigonometric_transform prismatic(const vector3<Coeff>& u, std::size_t k) {
    VARIETAS_ASSERT(k < N);

    trigonometric_transform r;
    const polynomial_type d = polynomial_type::variable(k);
    for (std::size_t i = 0; i < 3; ++i) {
      r.translation_[i] = polynomial_type::constant(u[i]) * d;
    }
    return r;
  }

  const polynomial_type& rotation(std::size_t i, std::size_t j) const {
    VARIETAS_ASSERT(i < 3 && j < 3);
    return rotation_[3 * i + j];
  }

  const polynomial_type& translation(std::size_t i) const {
    VARIETAS_ASSERT(i < 3);
    return translation_[i];
  }

  // Ordinary composition of rigid transforms, with no bookkeeping beside it.
  friend trigonometric_transform operator*(const trigonometric_transform& a,
                                           const trigonometric_transform& b) {
    trigonometric_transform r;
    for (std::size_t i = 0; i < 3; ++i) {
      for (std::size_t j = 0; j < 3; ++j) {
        polynomial_type sum;
        for (std::size_t k = 0; k < 3; ++k) {
          sum = sum + a.rotation(i, k) * b.rotation(k, j);
        }
        r.rotation_[3 * i + j] = sum;
      }

      polynomial_type sum = a.translation(i);
      for (std::size_t k = 0; k < 3; ++k) {
        sum = sum + a.rotation(i, k) * b.translation(k);
      }
      r.translation_[i] = sum;
    }
    return r;
  }

  // The transform at a point of the parameter space. The caller supplies the
  // cosine and sine themselves; a point off the circle is not a configuration,
  // and the result there is not a rigid transform, which is precisely what the
  // circle relations exist to exclude.
  rigid_transform<Coeff> evaluate(const std::array<Coeff, N>& point) const {
    matrix3<Coeff> rotation;
    vector3<Coeff> translation;
    for (std::size_t i = 0; i < 3; ++i) {
      for (std::size_t j = 0; j < 3; ++j) {
        rotation(i, j) = rotation_[3 * i + j].evaluate(point);
      }
      translation[i] = translation_[i].evaluate(point);
    }
    return rigid_transform<Coeff>(rotation, translation);
  }

 private:
  static Coeff skew_entry(const vector3<Coeff>& u, std::size_t i, std::size_t j) {
    if (i == j) {
      return traits::zero();
    }
    const std::size_t k = 3 - i - j;
    const bool cyclic = (j == (i + 1) % 3);
    return cyclic ? traits::negate(u[k]) : u[k];
  }

  std::array<polynomial_type, 9> rotation_{};
  std::array<polynomial_type, 3> translation_{};
};

// Forward kinematics as a polynomial map. As with the half-angle form the chain
// must already be validated and folded.
template <std::size_t N, class Order, class Coeff>
trigonometric_transform<Coeff, N, Order> trigonometric_forward_kinematics(
    const chain<Coeff>& robot) {
  using transform = trigonometric_transform<Coeff, N, Order>;
  VARIETAS_ASSERT(trigonometric_variable_count(robot) == N);

  transform running;
  std::size_t variable = 0;
  for (const joint<Coeff>& j : robot.joints()) {
    running = running * transform::constant(j.origin);
    switch (j.type) {
      case joint_type::revolute:
        running = running * transform::revolute(j.axis, variable);
        variable += 2;
        break;
      case joint_type::prismatic:
        running = running * transform::prismatic(j.axis, variable);
        variable += 1;
        break;
      case joint_type::fixed:
        break;
    }
  }
  return running * transform::constant(robot.tool());
}

// c_k^2 + s_k^2 - 1, one per revolute joint. These are generators of the ideal
// on the same footing as the pose residuals, not a correction applied to it.
template <std::size_t N, class Order, class Coeff>
std::vector<polynomial<Coeff, N, Order>> circle_relations(const chain<Coeff>& robot) {
  using poly = polynomial<Coeff, N, Order>;

  std::vector<poly> relations;
  for (const std::size_t k : cosine_variable_indices(robot)) {
    VARIETAS_ASSERT(k + 1 < N);
    relations.push_back(poly::variable(k, 2) + poly::variable(k + 1, 2) -
                        poly::constant(coefficient_traits<Coeff>::one()));
  }
  return relations;
}

// The point of the parameter space a configuration corresponds to, and the
// angles a point corresponds to. Unlike the half-angle map this is a bijection
// onto the product of circles with no exceptional configuration, which is why
// the recovery needs the two-argument arctangent and not a rational inverse.
template <class Coeff>
std::vector<Coeff> trigonometric_point(const chain<Coeff>& robot,
                                       const std::vector<Coeff>& cosines,
                                       const std::vector<Coeff>& sines) {
  const auto indices = cosine_variable_indices(robot);
  VARIETAS_ASSERT(cosines.size() == indices.size());
  VARIETAS_ASSERT(sines.size() == indices.size());

  std::vector<Coeff> point(trigonometric_variable_count(robot),
                           coefficient_traits<Coeff>::zero());
  for (std::size_t i = 0; i < indices.size(); ++i) {
    point[indices[i]] = cosines[i];
    point[indices[i] + 1] = sines[i];
  }
  return point;
}

inline double angle_from_cosine_sine(double c, double s) { return std::atan2(s, c); }

}  // namespace varietas

#endif
