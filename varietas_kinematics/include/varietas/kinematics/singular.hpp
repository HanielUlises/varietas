#ifndef VARIETAS_KINEMATICS_SINGULAR_HPP
#define VARIETAS_KINEMATICS_SINGULAR_HPP

#include <array>
#include <cstddef>
#include <vector>

#include "varietas/core/config.hpp"
#include "varietas/core/embed.hpp"
#include "varietas/core/ideal/buchberger.hpp"
#include "varietas/core/ideal/dimension.hpp"
#include "varietas/core/ideal/ideal.hpp"
#include "varietas/core/minors.hpp"
#include "varietas/core/order/grevlex.hpp"
#include "varietas/core/polynomial.hpp"
#include "varietas/kinematics/chain.hpp"
#include "varietas/kinematics/trigonometric.hpp"
#include "varietas/kinematics/workspace.hpp"

namespace varietas {

// The singular locus, as a variety.
//
// A configuration is singular when the differential of the forward kinematics
// map drops rank there: the arm loses, instantaneously, the ability to move its
// tool in some direction, and every numerical inverse kinematics method in
// existence degrades in a neighbourhood of it. The usual treatment is a
// condition number: evaluate the Jacobian at a configuration, take its
// smallest singular value, and call anything below a threshold "near
// singular". That reports a number about one configuration. It does not say
// what the singular set *is*, how many pieces it has, or where they go in the
// workspace, and it cannot, because those are questions about a variety and a
// condition number is a measurement at a point.
//
// Algebraically the question is elementary. Rank is not a polynomial condition
// (the rank of a matrix is not a continuous function of its entries at all),
// but rank *deficiency* is: a matrix has rank less than k exactly when every
// k by k minor vanishes. The minors of the Jacobian are polynomials in the
// joint variables, and adjoined to the equations of the parameter space they
// generate the ideal whose variety is the singular set. Everything the library
// already has then applies to it: dimension says whether the singularities form
// a surface or a curve or a finite set, splitting separates the branches, and
// elimination pushes the set forward into the workspace, where it is the
// obstruction the arm's motion planner actually meets.
//
// This is the second consumer of the trigonometric formulation and the one that
// settles the choice. The construction is elimination-shaped, so it inherits
// the cost argument from workspace.hpp; but it has its own reason as well. The
// differential of a map restricted to a variety is not the differential of the
// polynomials that define it. Over the half-angle ring one would need the
// Jacobian of the cleared numerators, corrected for the denominators that were
// multiplied through, and the correction is a quotient rule in every entry. In
// the trigonometric ring the differential is the geometric Jacobian, unchanged
// from the textbook, and it is already polynomial in c and s. The reason is
// worth stating exactly, since it is what makes the whole header legitimate:
// along the circle c^2 + s^2 = 1 the tangent direction at (c, s) is (-s, c),
// which is precisely (dc/dq, ds/dq), so differentiating the polynomial map
// along the constraint reproduces d/dq, and the geometric Jacobian is the
// differential of the map restricted to the parameter variety, with no
// Lagrange multiplier and no correction term.

// Which rows of the geometric Jacobian the task uses. The rows are the six
// components of the tool twist: 0-2 the linear velocity, 3-5 the angular.
//
// The choice is not cosmetic. An arm is singular *for a task*, and asking for
// the wrong one is how a singularity analysis comes to disagree with the robot:
// a planar three-link arm is nowhere singular for the position task alone,
// since two joints suffice to move a point in the plane, and is singular on a
// whole surface for the planar pose task that also asks for the tool's
// heading. The cost is chosen here too, since the number of minors is
// C(rows, k) C(joints, k).
inline std::vector<std::size_t> position_rows() { return {0, 1, 2}; }

inline std::vector<std::size_t> orientation_rows() { return {3, 4, 5}; }

inline std::vector<std::size_t> pose_rows() { return {0, 1, 2, 3, 4, 5}; }

// Position in the plane and heading within it, for an arm that moves in z = 0
// about axes parallel to z.
inline std::vector<std::size_t> planar_pose_rows() { return {0, 1, 5}; }

// The geometric Jacobian as a matrix of polynomials: six rows, one column per
// actuated joint.
//
// Column i of a revolute joint is the twist a_i × (p - p_i), a_i, with a_i the
// joint axis in the base frame, p_i a point on it and p the tool position; a
// prismatic joint contributes a_i, 0. Every one of those is polynomial in the
// trigonometric variables, since a_i and p_i are the axis and origin of the
// frame the joint acts in and that frame is a product of the polynomial
// matrices trigonometric_transform composes.
//
// Note which frame each column is read in: a_i and p_i are taken *before* the
// joint's own motion is applied. That is not an approximation but an identity:
// a rotation fixes its own axis, and the origin of the joint frame does not
// move under it. It is why the walk below multiplies the origin in, reads
// the column, and only then applies the joint.
template <std::size_t V, class Order, class Coeff>
dense_matrix<polynomial<Coeff, V, Order>> trigonometric_jacobian(const chain<Coeff>& robot) {
  using poly = polynomial<Coeff, V, Order>;
  using transform = trigonometric_transform<Coeff, V, Order>;
  VARIETAS_ASSERT(trigonometric_variable_count(robot) == V);

  const transform map = trigonometric_forward_kinematics<V, Order>(robot);
  const std::array<poly, 3> tool = {map.translation(0), map.translation(1),
                                    map.translation(2)};

  dense_matrix<poly> jacobian(6, robot.degrees_of_freedom());

  transform running;
  std::size_t variable = 0;
  std::size_t column = 0;
  for (const joint<Coeff>& j : robot.joints()) {
    running = running * transform::constant(j.origin);
    if (j.type == joint_type::fixed) {
      continue;
    }

    // The joint axis in the base frame, and the origin of its frame.
    std::array<poly, 3> axis;
    for (std::size_t i = 0; i < 3; ++i) {
      poly sum;
      for (std::size_t k = 0; k < 3; ++k) {
        sum = sum + running.rotation(i, k) * poly::constant(j.axis[k]);
      }
      axis[i] = sum;
    }

    if (j.type == joint_type::revolute) {
      const std::array<poly, 3> arm = {tool[0] - running.translation(0),
                                       tool[1] - running.translation(1),
                                       tool[2] - running.translation(2)};
      for (std::size_t i = 0; i < 3; ++i) {
        const std::size_t a = (i + 1) % 3;
        const std::size_t b = (i + 2) % 3;
        jacobian(i, column) = axis[a] * arm[b] - axis[b] * arm[a];
        jacobian(3 + i, column) = axis[i];
      }
      running = running * transform::revolute(j.axis, variable);
      variable += 2;
    } else {
      for (std::size_t i = 0; i < 3; ++i) {
        jacobian(i, column) = axis[i];
        jacobian(3 + i, column) = poly();
      }
      running = running * transform::prismatic(j.axis, variable);
      variable += 1;
    }
    ++column;
  }

  return jacobian;
}

// Generators of the ideal of the singular locus, in the joint variables: the
// circle relations, which say that the point is a configuration at all, and the
// maximal minors of the task rows of the Jacobian, which say that the arm is
// singular there.
//
// The minor size is min(rows, joints), so a redundant arm, one with more joints
// than the task has dimensions, is singular where it cannot span the task space,
// which is the right statement and the one the maximal minors make without
// being told.
template <std::size_t V, class Order, class Coeff>
std::vector<polynomial<Coeff, V, Order>> singular_generators(
    const chain<Coeff>& robot, const std::vector<std::size_t>& rows) {
  using poly = polynomial<Coeff, V, Order>;

  const dense_matrix<poly> jacobian = trigonometric_jacobian<V, Order>(robot);

  std::vector<std::size_t> columns(jacobian.cols());
  for (std::size_t i = 0; i < columns.size(); ++i) {
    columns[i] = i;
  }
  for (const std::size_t r : rows) {
    VARIETAS_ASSERT(r < 6);
  }

  std::vector<poly> generators = circle_relations<V, Order>(robot);
  for (poly& minor : maximal_minors(jacobian.submatrix(rows, columns))) {
    generators.push_back(std::move(minor));
  }
  return generators;
}

// The reduced Gröbner basis of that ideal.
//
// An empty variety, reported by is_unit_ideal and equivalently by
// ideal_dimension's is_empty, means the arm has no singular configuration at
// all, over the algebraic closure and therefore in particular over the reals.
// The converse does not hold and the tests make a point of it: a nonempty
// singular variety may have no real points, in which case the arm is in fact
// everywhere nonsingular and the ideal is merely recording where it would be
// singular if its joint angles were allowed to be complex. Reality is a
// property of points, not of ideals, exactly as it was for unreachability.
template <std::size_t V, class Coeff>
std::vector<polynomial<Coeff, V, grevlex>> singular_ideal(
    const chain<Coeff>& robot, const std::vector<std::size_t>& rows,
    buchberger_statistics* statistics = nullptr) {
  return groebner_basis(singular_generators<V, grevlex>(robot, rows), statistics);
}

// The dimension of the singular locus, as a set of configurations. Zero means
// isolated singular configurations, one a curve of them, and the independent
// set names which combination of joint variables is free to move along it.
template <std::size_t V, class Coeff>
affine_dimension<V> singular_dimension(const chain<Coeff>& robot,
                                       const std::vector<std::size_t>& rows) {
  return ideal_dimension(singular_ideal<V>(robot, rows));
}

// The image of the singular locus in the workspace: the equations satisfied by
// every tool position the arm can reach singularly.
//
// This is the construction the singular locus is for. A set of bad
// configurations is of limited use to a planner, which works in the space the
// task is posed in; the image of that set is the surface the arm cannot cross
// without passing through a singularity, and for an ordinary arm it is the
// boundary of the reachable workspace. The two-link planar arm below returns
// the circle of radius l1 + l2 and the point at radius |l1 - l2|, which is
// exactly the outer boundary and the inner hole.
//
// It is the same elimination workspace_relations performs, with the minors
// adjoined to the graph before the joint variables are eliminated. The Closure
// Theorem applies verbatim and carries the same warning: what comes back cuts
// out the Zariski closure of the singular image, so a semialgebraic piece of it
// (an arc rather than a whole circle) is not something these equations can
// distinguish.
template <std::size_t V, class Coeff>
std::vector<polynomial<Coeff, 3, grevlex>> singular_workspace_relations(
    const chain<Coeff>& robot, const std::vector<std::size_t>& rows,
    buchberger_statistics* statistics = nullptr) {
  using order = trigonometric_workspace_order<V>;
  using layout = trigonometric_workspace_layout<V>;
  using enlarged = trigonometric_workspace_polynomial<Coeff, V>;

  std::vector<enlarged> generators = trigonometric_workspace_generators<V>(robot);

  // The circle relations are already among the workspace generators; only the
  // minors have to be lifted into the enlarged ring.
  const dense_matrix<polynomial<Coeff, V, grevlex>> jacobian =
      trigonometric_jacobian<V, grevlex>(robot);
  std::vector<std::size_t> columns(jacobian.cols());
  for (std::size_t i = 0; i < columns.size(); ++i) {
    columns[i] = i;
  }
  for (const auto& minor : maximal_minors(jacobian.submatrix(rows, columns))) {
    generators.push_back(embed_polynomial<order, V + 3>(minor, layout::first_joint));
  }

  const auto graph = groebner_basis(generators, statistics);
  const auto eliminated = eliminated_generators(graph, layout::eliminated_block);

  std::vector<polynomial<Coeff, 3, grevlex>> relations;
  relations.reserve(eliminated.size());
  for (const auto& g : eliminated) {
    relations.push_back(project_polynomial<grevlex, 3>(g, layout::first_pose));
  }
  reduce_groebner_basis(relations);
  return relations;
}

}  // namespace varietas

#endif
