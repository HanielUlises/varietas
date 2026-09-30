#ifndef VARIETAS_IK_RECONSTRUCTED_IK_HPP
#define VARIETAS_IK_RECONSTRUCTED_IK_HPP

#include <array>
#include <cstddef>
#include <optional>
#include <random>
#include <string>
#include <vector>

#include "varietas/codegen/parametric_solution.hpp"
#include "varietas/codegen/prime_field.hpp"
#include "varietas/codegen/rational.hpp"
#include "varietas/codegen/rational_function.hpp"
#include "varietas/codegen/reconstruct.hpp"
#include "varietas/core/ideal/division.hpp"
#include "varietas/core/monomial.hpp"
#include "varietas/core/order/grevlex.hpp"
#include "varietas/core/order/order_id.hpp"
#include "varietas/core/quotient/quotient_basis.hpp"
#include "varietas/ik/cast_chain.hpp"
#include "varietas/ik/parametric_ik.hpp"
#include "varietas/kinematics/chain.hpp"
#include "varietas/kinematics/rationalize.hpp"

// The parametric position solve, recovered from fixed poses instead of carried
// out over Q(p).
//
// parametric_position_ik computes one Grobner basis over the field of rational
// functions of the pose. With the modular gcd that is feasible for three
// parameters, but on arms with offsets the rational functions it carries on
// the way swell to thousands of terms, while the action matrices it finally
// produces have entries of a dozen. This computes the same matrices without
// ever forming the intermediate functions: it solves the arm at many poses,
// each one an ordinary Grobner basis over a prime field with no parameters,
// and recovers every matrix entry as a rational function of the pose from its
// values (varietas/codegen/reconstruct.hpp).
//
// The result is the same codegen::parametric_solution the symbolic solve
// produces, so emit() writes the same kind of header from it, and on the arms
// where both finish the two agree entry for entry.
//
// What is given up is certainty. The symbolic solve is exact; this one is right
// with high probability, since reconstruction from samples has no final exact
// test. So the result is checked here, exactly, at a few rational poses: the
// fixed-pose basis is computed over Q there and its action matrices compared
// with the reconstructed ones evaluated at the same pose. A disagreement
// anywhere is reported as a failure rather than returned.
namespace varietas {
namespace ik {

struct reconstruction_report {
  reconstruction::statistics sampling;
  std::size_t exact_checks = 0;  // rational poses the result was checked at
  bool checks_passed = false;
};

namespace detail {

// The standard monomials and the action data at one pose, over whatever field
// the chain is written in.
template <std::size_t N, class Coeff>
struct pose_solution {
  quotient_basis<N> quotient;
  // N matrices of d x d, column-major as in parametric_matrix, then N rows of
  // d variable coordinates.
  std::vector<Coeff> values;
};

template <std::size_t N, std::size_t P, class Coeff>
bool solve_at_pose(const chain<Coeff>& robot,
                   const rational_transform<Coeff, N, grevlex>& map,
                   const std::array<std::size_t, P>& coordinates,
                   const std::array<Coeff, P>& pose, pose_solution<N, Coeff>& out) {
  using poly = polynomial<Coeff, N, grevlex>;
  using traits = coefficient_traits<Coeff>;
  const poly denominator = map.denominator();
  std::vector<poly> residuals;
  residuals.reserve(P);
  for (std::size_t k = 0; k < P; ++k) {
    residuals.push_back(map.translation(coordinates[k]) - denominator * poly::constant(pose[k]));
  }
  const auto basis = kinematic_ideal_generators<N, grevlex>(robot, residuals);
  out.quotient = standard_monomials(basis);
  if (!out.quotient.is_zero_dimensional || out.quotient.dimension() == 0) {
    return false;
  }
  const std::size_t d = out.quotient.dimension();
  out.values.assign(N * d * d + N * d, traits::zero());
  for (std::size_t i = 0; i < N; ++i) {
    for (std::size_t j = 0; j < d; ++j) {
      const poly image =
          normal_form(poly::variable(i) * poly::from_monomial(out.quotient.monomials[j], traits::one()),
                      basis);
      for (const auto& t : image.terms()) {
        const std::size_t row = out.quotient.index_of(t.mon);
        VARIETAS_ASSERT(row < d);
        out.values[i * d * d + j * d + row] = t.coeff;
      }
    }
    const poly coordinate = normal_form(poly::variable(i), basis);
    for (const auto& t : coordinate.terms()) {
      const std::size_t k = out.quotient.index_of(t.mon);
      VARIETAS_ASSERT(k < d);
      out.values[N * d * d + i * d + k] = t.coeff;
    }
  }
  return true;
}

}  // namespace detail

template <std::size_t N, std::size_t P>
parametric_ik_result<N, P> reconstructed_position_ik(const chain<rational>& robot,
                                                     const std::array<std::size_t, P>& coordinates,
                                                     reconstruction_report* report = nullptr,
                                                     std::size_t exact_checks = 2) {
  static_assert(P >= 1 && P <= 3, "a tool position has three coordinates");
  using field = rational_function<P>;

  parametric_ik_result<N, P> result;
  reconstruction_report local;
  reconstruction_report& rep = report != nullptr ? *report : local;

  // The same refusals, in the same order, as the symbolic solve.
  if (P < N) {
    result.status = parametric_ik_status::underdetermined;
    return result;
  }
  if (P > N) {
    result.status = parametric_ik_status::overdetermined;
    return result;
  }
  std::array<bool, 3> constrained{};
  for (std::size_t k = 0; k < P; ++k) {
    if (coordinates[k] >= 3 || constrained[coordinates[k]]) {
      result.status = parametric_ik_status::bad_coordinates;
      result.offending_coordinate = coordinates[k];
      return result;
    }
    constrained[coordinates[k]] = true;
  }
  const chain<rational> exact = robot.fold_fixed_joints();
  if (exact.degrees_of_freedom() != N) {
    result.status = parametric_ik_status::wrong_degrees_of_freedom;
    return result;
  }
  const auto exact_map = rational_forward_kinematics<N, grevlex>(exact);
  for (std::size_t i = 0; i < 3; ++i) {
    if (!constrained[i] && !exact_map.translation(i).is_zero()) {
      result.status = parametric_ik_status::dropped_coordinate_is_not_identically_zero;
      result.offending_coordinate = i;
      return result;
    }
  }

  // The shape of the solution set at a general pose, from three rational poses
  // drawn at random: the standard monomials two of them agree on. A pose on a
  // special fibre gives a different shape, and two such poses agreeing with
  // each other is not a case worth planning for.
  std::mt19937_64 rng(0x1dea1);
  std::uniform_int_distribution<long> small(-40, 40);
  const auto random_rational_pose = [&] {
    std::array<rational, P> pose;
    for (auto& c : pose) {
      c = make_rational(small(rng), 7 + (small(rng) + 40) % 11);
    }
    return pose;
  };
  std::vector<detail::pose_solution<N, rational>> exact_samples;
  std::vector<std::array<rational, P>> exact_poses;
  for (int i = 0; i < 3; ++i) {
    exact_poses.push_back(random_rational_pose());
    detail::pose_solution<N, rational> sample;
    if (!detail::solve_at_pose<N, P>(exact, exact_map, coordinates, exact_poses.back(), sample)) {
      sample.quotient.monomials.clear();
    }
    exact_samples.push_back(std::move(sample));
  }
  std::optional<quotient_basis<N>> generic;
  for (int i = 0; i < 3 && !generic; ++i) {
    for (int j = i + 1; j < 3; ++j) {
      if (!exact_samples[i].quotient.monomials.empty() &&
          exact_samples[i].quotient.monomials == exact_samples[j].quotient.monomials) {
        generic = exact_samples[i].quotient;
        break;
      }
    }
  }
  if (!generic) {
    result.status = parametric_ik_status::not_zero_dimensional;
    return result;
  }
  const std::size_t d = generic->dimension();
  const std::size_t count = N * d * d + N * d;

  // The black box, over whichever prime is selected.
  chain<residue> modular_chain;
  std::optional<rational_transform<residue, N, grevlex>> modular_map;
  const reconstruction::prime_hook prepare = [&](reconstruction::word) {
    modular_chain = cast_chain<residue>(exact);
    if (residue::undefined_image()) {
      return false;
    }
    modular_map = rational_forward_kinematics<N, grevlex>(modular_chain);
    return true;
  };
  const reconstruction::black_box<P> box = [&](const reconstruction::point<P>& pose,
                                               std::vector<residue>& values) {
    detail::pose_solution<N, residue> sample;
    if (!detail::solve_at_pose<N, P>(modular_chain, *modular_map, coordinates, pose, sample)) {
      return false;
    }
    if (sample.quotient.monomials != generic->monomials) {
      return false;  // a special fibre
    }
    values = std::move(sample.values);
    return true;
  };

  std::vector<field> functions;
  if (!reconstruction::reconstruct<P>(box, prepare, count, functions, &rep.sampling)) {
    result.status = parametric_ik_status::not_zero_dimensional;
    return result;
  }
  residue::select(modular_detail::nth_prime_below_2_31(0));

  // The exact check, at the rational poses already solved plus fresh ones if
  // more were asked for.
  rep.checks_passed = true;
  rep.exact_checks = 0;
  for (std::size_t c = 0; c < exact_checks; ++c) {
    if (c >= exact_samples.size()) {
      exact_poses.push_back(random_rational_pose());
      detail::pose_solution<N, rational> sample;
      detail::solve_at_pose<N, P>(exact, exact_map, coordinates, exact_poses.back(), sample);
      exact_samples.push_back(std::move(sample));
    }
    const auto& sample = exact_samples[c];
    if (sample.quotient.monomials != generic->monomials) {
      continue;  // a special pose; it says nothing about the general one
    }
    for (std::size_t f = 0; f < count && rep.checks_passed; ++f) {
      const auto& g = functions[f];
      if (!g.is_zero() && sgn(g.denominator().evaluate(exact_poses[c])) == 0) {
        rep.checks_passed = false;  // a pole of the result where Q has a value
        break;
      }
      const rational value = g.is_zero() ? rational(0) : g.evaluate(exact_poses[c]);
      if (value != sample.values[f]) {
        rep.checks_passed = false;
      }
    }
    ++rep.exact_checks;
  }
  if (!rep.checks_passed || rep.exact_checks == 0) {
    result.status = parametric_ik_status::not_zero_dimensional;
    return result;
  }

  // Assembled exactly as the symbolic solve assembles it.
  result.branches = d;
  auto& solution = result.solution;
  solution.order = order_id::grevlex;
  solution.quotient = *generic;
  solution.one_index = generic->index_of(monomial<N>::one());
  for (const joint<rational>& j : exact.joints()) {
    if (j.is_actuated()) {
      solution.unknown_names.push_back(detail::variable_name(j));
    }
  }
  for (std::size_t k = 0; k < P; ++k) {
    solution.parameter_names.push_back(detail::coordinate_name(coordinates[k]));
  }
  for (std::size_t i = 0; i < N; ++i) {
    codegen::parametric_matrix<P> m;
    m.dimension = d;
    m.entries.assign(functions.begin() + static_cast<std::ptrdiff_t>(i * d * d),
                     functions.begin() + static_cast<std::ptrdiff_t>((i + 1) * d * d));
    solution.action.push_back(std::move(m));
  }
  solution.variable_coordinates.assign(N, std::vector<field>(d));
  for (std::size_t i = 0; i < N; ++i) {
    for (std::size_t k = 0; k < d; ++k) {
      solution.variable_coordinates[i][k] = functions[N * d * d + i * d + k];
    }
  }
  for (std::size_t k = 0; k < P; ++k) {
    solution.residual_numerators.push_back(exact_map.translation(coordinates[k]));
  }
  solution.residual_denominator = exact_map.denominator();
  VARIETAS_ASSERT(solution.is_well_formed());
  return result;
}

}  // namespace ik
}  // namespace varietas

#endif
