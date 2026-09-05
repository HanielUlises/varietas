#ifndef VARIETAS_KINEMATICS_RATIONALIZE_HPP
#define VARIETAS_KINEMATICS_RATIONALIZE_HPP

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>

#include "varietas/core/config.hpp"
#include "varietas/core/ideal/saturation.hpp"
#include "varietas/core/polynomial.hpp"
#include "varietas/kinematics/chain.hpp"
#include "varietas/kinematics/rigid_transform.hpp"

namespace varietas {

// The tangent half-angle substitution, applied to a whole chain.
//
// Forward kinematics is transcendental as stated: it needs the cosine and sine
// of a real angle. Writing t_i = tan(q_i / 2) makes it rational,
//
//     cos q = (1 - t^2) / (1 + t^2),    sin q = 2t / (1 + t^2),
//
// and Rodrigues' formula for a rotation about a unit axis u then reads
//
//     R = I + sin q [u] + (1 - cos q) [u]^2
//       = (1 / (1 + t^2)) [ (1 + t^2) I + 2t [u] + 2t^2 [u]^2 ],
//
// a matrix of quadratics over a single denominator 1 + t^2. That shape is what
// makes this layer cheap, and it is worth stating precisely because it is the
// property the whole representation rests on: every factor of the chain is a
// polynomial matrix over a denominator of the form (1 + t_i^2)^{e_i}, and that
// form is closed under composition. So a transform can be carried as a
// numerator together with the exponent vector e, and no general rational
// function arithmetic, with no common denominators to compute and no gcd, is
// required.
//
// What this deliberately is not is an element of SE(3) over a field. The
// numerator entries are polynomials, and polynomials are not invertible, so
// matrix3<Coeff> cannot be reused: it asks coefficient_traits for inverses.
// Having to write a separate type is the type system saying something true.

// A rigid transform whose entries are rational functions of the joint
// variables, held as numerator over ∏ (1 + t_i^2)^{e_i}.
template <class Coeff, std::size_t N, class Order>
class rational_transform {
 public:
  using polynomial_type = polynomial<Coeff, N, Order>;
  using traits = coefficient_traits<Coeff>;
  using exponent_type = std::uint16_t;

  static constexpr std::size_t num_vars = N;

  rational_transform() {
    for (std::size_t i = 0; i < 3; ++i) {
      rotation_[3 * i + i] = polynomial_type::constant(traits::one());
    }
  }

  static rational_transform identity() { return rational_transform(); }

  // A fixed placement, which is a rational transform with denominator one.
  static rational_transform constant(const rigid_transform<Coeff>& t) {
    rational_transform r;
    for (std::size_t i = 0; i < 3; ++i) {
      for (std::size_t j = 0; j < 3; ++j) {
        r.rotation_[3 * i + j] = polynomial_type::constant(t.rotation()(i, j));
      }
      r.translation_[i] = polynomial_type::constant(t.translation()[i]);
    }
    return r;
  }

  // The displacement of a revolute joint about the unit axis u, parameterised
  // by t = tan(q / 2) in the variable at index k. The numerator is Rodrigues'
  // formula cleared of its denominator, using [u]^2 = u u^T - I, which holds
  // because u is a unit vector, a fact chain::validate has already checked.
  static rational_transform revolute(const vector3<Coeff>& u, std::size_t k) {
    VARIETAS_ASSERT(k < N);

    const polynomial_type t = polynomial_type::variable(k);
    const polynomial_type t2 = polynomial_type::variable(k, 2);
    const polynomial_type one = polynomial_type::constant(traits::one());
    const Coeff two = traits::one() + traits::one();

    rational_transform r;
    for (std::size_t i = 0; i < 3; ++i) {
      for (std::size_t j = 0; j < 3; ++j) {
        // [u]^2 = u u^T - I.
        Coeff outer = u[i] * u[j];
        if (i == j) {
          outer = outer - traits::one();
        }
        // 2t [u], the skew part.
        const Coeff skew = skew_entry(u, i, j);

        polynomial_type entry = (two * t) * polynomial_type::constant(skew) +
                                (two * t2) * polynomial_type::constant(outer);
        if (i == j) {
          entry = entry + one + t2;
        }
        r.rotation_[3 * i + j] = entry;
      }
    }
    r.denominator_[k] = 1;
    return r;
  }

  // A prismatic joint contributes its displacement directly: the variable is
  // the displacement itself, no substitution is involved, and the denominator
  // stays one.
  static rational_transform prismatic(const vector3<Coeff>& u, std::size_t k) {
    VARIETAS_ASSERT(k < N);

    rational_transform r;
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

  const std::array<exponent_type, N>& denominator_exponents() const noexcept {
    return denominator_;
  }

  // ∏ (1 + t_i^2)^{e_i}, expanded. Needed to clear the denominators when the
  // pose equations are formed, and as the multiplier the ideal is saturated by.
  polynomial_type denominator() const {
    polynomial_type d = polynomial_type::constant(traits::one());
    for (std::size_t i = 0; i < N; ++i) {
      for (exponent_type e = 0; e < denominator_[i]; ++e) {
        d = d * half_angle_denominator<Coeff, N, Order>(i);
      }
    }
    return d;
  }

  // Composition in SE(3), carried through the representation. Writing
  // T = (1/D)[R | p], the product (1/Da)[Ra | pa] · (1/Db)[Rb | pb] has
  // rotation Ra Rb over Da Db and translation (Ra pb + pa Db) over Da Db, so
  // the numerators multiply and the exponent vectors add. Nothing else can
  // arise, which is the closure property the whole design rests on.
  friend rational_transform operator*(const rational_transform& a,
                                      const rational_transform& b) {
    rational_transform r;
    const polynomial_type db = b.denominator();

    for (std::size_t i = 0; i < 3; ++i) {
      for (std::size_t j = 0; j < 3; ++j) {
        polynomial_type sum;
        for (std::size_t k = 0; k < 3; ++k) {
          sum = sum + a.rotation(i, k) * b.rotation(k, j);
        }
        r.rotation_[3 * i + j] = sum;
      }

      polynomial_type sum = a.translation(i) * db;
      for (std::size_t k = 0; k < 3; ++k) {
        sum = sum + a.rotation(i, k) * b.translation(k);
      }
      r.translation_[i] = sum;
    }

    for (std::size_t i = 0; i < N; ++i) {
      r.denominator_[i] =
          static_cast<exponent_type>(a.denominator_[i] + b.denominator_[i]);
    }
    return r;
  }

  // The transform at a configuration, by evaluating the rational functions.
  // This is what ties the algebra back to the map it represents: substituting
  // t_i = tan(q_i / 2) must reproduce the numerical forward kinematics.
  rigid_transform<Coeff> evaluate(const std::array<Coeff, N>& point) const {
    Coeff d = traits::one();
    for (std::size_t i = 0; i < N; ++i) {
      Coeff factor = traits::one() + point[i] * point[i];
      for (exponent_type e = 0; e < denominator_[i]; ++e) {
        d = d * factor;
      }
    }
    VARIETAS_ASSERT(!traits::is_zero(d));
    const Coeff inverse = traits::inverse(d);

    matrix3<Coeff> rotation;
    vector3<Coeff> translation;
    for (std::size_t i = 0; i < 3; ++i) {
      for (std::size_t j = 0; j < 3; ++j) {
        rotation(i, j) = rotation_[3 * i + j].evaluate(point) * inverse;
      }
      translation[i] = translation_[i].evaluate(point) * inverse;
    }
    return rigid_transform<Coeff>(rotation, translation);
  }

 private:
  static Coeff skew_entry(const vector3<Coeff>& u, std::size_t i, std::size_t j) {
    // [u] is the matrix with rows (0, -uz, uy), (uz, 0, -ux), (-uy, ux, 0).
    if (i == j) {
      return traits::zero();
    }
    const std::size_t k = 3 - i - j;  // the index neither i nor j
    // Sign is +1 on the cyclic permutations (0,1,2), (1,2,0), (2,0,1).
    const bool cyclic = (j == (i + 1) % 3);
    return cyclic ? traits::negate(u[k]) : u[k];
  }

  std::array<polynomial_type, 9> rotation_{};
  std::array<polynomial_type, 3> translation_{};
  std::array<exponent_type, N> denominator_{};
};

// Forward kinematics as a rational map. The chain must already have been
// validated and folded; fixed joints left in place would multiply a constant
// matrix into every generator for no reason, which is what fold_fixed_joints
// exists to prevent.
template <std::size_t N, class Order, class Coeff>
rational_transform<Coeff, N, Order> rational_forward_kinematics(
    const chain<Coeff>& robot) {
  VARIETAS_ASSERT(robot.degrees_of_freedom() == N);

  rational_transform<Coeff, N, Order> running;
  std::size_t variable = 0;
  for (const joint<Coeff>& j : robot.joints()) {
    running = running * rational_transform<Coeff, N, Order>::constant(j.origin);
    switch (j.type) {
      case joint_type::revolute:
        running = running * rational_transform<Coeff, N, Order>::revolute(j.axis, variable++);
        break;
      case joint_type::prismatic:
        running = running * rational_transform<Coeff, N, Order>::prismatic(j.axis, variable++);
        break;
      case joint_type::fixed:
        break;
    }
  }
  return running * rational_transform<Coeff, N, Order>::constant(robot.tool());
}

// The residuals of the pose equations, denominators cleared: twelve
// polynomials, nine from the rotation and three from the translation. The
// system is overdetermined, which is the usual formulation and which
// Buchberger handles without special treatment: the redundancy is information
// about the rotation group, not an inconsistency.
template <std::size_t N, class Order, class Coeff>
std::vector<polynomial<Coeff, N, Order>> pose_residuals(
    const rational_transform<Coeff, N, Order>& map,
    const rigid_transform<Coeff>& target) {
  using poly = polynomial<Coeff, N, Order>;

  const poly denominator = map.denominator();

  std::vector<poly> residuals;
  residuals.reserve(12);
  for (std::size_t i = 0; i < 3; ++i) {
    for (std::size_t j = 0; j < 3; ++j) {
      residuals.push_back(map.rotation(i, j) -
                          denominator * poly::constant(target.rotation()(i, j)));
    }
    residuals.push_back(map.translation(i) -
                        denominator * poly::constant(target.translation()[i]));
  }
  return residuals;
}

// The residuals of the position equations alone. Asking only where the tool is,
// and not how it is turned, is the classical inverse kinematics problem for an
// arm with fewer joints than SE(3) has dimensions: a planar two-link arm has
// two elbow solutions for a reachable point, and none of them for the pose,
// since the orientation is then whatever the position leaves it.
template <std::size_t N, class Order, class Coeff>
std::vector<polynomial<Coeff, N, Order>> position_residuals(
    const rational_transform<Coeff, N, Order>& map, const vector3<Coeff>& target) {
  using poly = polynomial<Coeff, N, Order>;

  const poly denominator = map.denominator();

  std::vector<poly> residuals;
  residuals.reserve(3);
  for (std::size_t i = 0; i < 3; ++i) {
    residuals.push_back(map.translation(i) - denominator * poly::constant(target[i]));
  }
  return residuals;
}

// One half-angle denominator per revolute joint. A prismatic joint contributes
// none, having needed no substitution.
template <std::size_t N, class Order, class Coeff>
std::vector<polynomial<Coeff, N, Order>> half_angle_multipliers(
    const chain<Coeff>& robot) {
  std::vector<polynomial<Coeff, N, Order>> multipliers;
  std::size_t variable = 0;
  for (const joint<Coeff>& j : robot.joints()) {
    if (!j.is_actuated()) {
      continue;
    }
    if (j.type == joint_type::revolute) {
      multipliers.push_back(half_angle_denominator<Coeff, N, Order>(variable));
    }
    ++variable;
  }
  return multipliers;
}

// The saturated ideal of a set of kinematic residuals.
//
// Clearing the denominators attaches to the variety the loci 1 + t_i^2 = 0,
// that is t_i = ±i, which are the images of the configurations q_i = π that the
// substitution sent to infinity. They are not configurations of the arm. Over
// the reals they are invisible, and that is exactly why they have to be removed
// here rather than noticed later: dim_k A counts them regardless, and that
// count is the certificate that no branch of the inverse kinematics was missed.
template <std::size_t N, class Order, class Coeff>
std::vector<polynomial<Coeff, N, Order>> kinematic_ideal_generators(
    const chain<Coeff>& robot, const std::vector<polynomial<Coeff, N, Order>>& residuals,
    buchberger_statistics* statistics = nullptr) {
  const auto multipliers = half_angle_multipliers<N, Order>(robot);
  if (multipliers.empty()) {
    return groebner_basis(residuals, statistics);
  }
  return saturate_by_product(residuals, multipliers, statistics);
}

// Generators of the inverse kinematics ideal at a full pose.
template <std::size_t N, class Order, class Coeff>
std::vector<polynomial<Coeff, N, Order>> pose_ideal_generators(
    const chain<Coeff>& robot, const rigid_transform<Coeff>& target,
    buchberger_statistics* statistics = nullptr) {
  const auto map = rational_forward_kinematics<N, Order>(robot);
  return kinematic_ideal_generators<N, Order>(robot, pose_residuals(map, target),
                                              statistics);
}

// Generators of the inverse kinematics ideal at a tool position.
template <std::size_t N, class Order, class Coeff>
std::vector<polynomial<Coeff, N, Order>> position_ideal_generators(
    const chain<Coeff>& robot, const vector3<Coeff>& target,
    buchberger_statistics* statistics = nullptr) {
  const auto map = rational_forward_kinematics<N, Order>(robot);
  return kinematic_ideal_generators<N, Order>(robot, position_residuals(map, target),
                                              statistics);
}

// The angle a joint variable stands for, and its inverse. The substitution is
// a bijection away from q = π, which is the point it sends to infinity and
// which saturation removes; recovering an angle from a solution of the variety
// is therefore always well defined.
inline double angle_from_variable(double t) { return 2.0 * std::atan(t); }
inline double variable_from_angle(double q) { return std::tan(0.5 * q); }

}  // namespace varietas

#endif
