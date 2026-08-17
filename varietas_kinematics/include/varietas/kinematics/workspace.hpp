#ifndef VARIETAS_KINEMATICS_WORKSPACE_HPP
#define VARIETAS_KINEMATICS_WORKSPACE_HPP

#include <cstddef>
#include <vector>

#include "varietas/core/config.hpp"
#include "varietas/core/embed.hpp"
#include "varietas/core/ideal/ideal.hpp"
#include "varietas/core/ideal/saturation.hpp"
#include "varietas/core/order/block.hpp"
#include "varietas/core/order/grevlex.hpp"
#include "varietas/core/polynomial.hpp"
#include "varietas/kinematics/chain.hpp"
#include "varietas/kinematics/rationalize.hpp"
#include "varietas/kinematics/trigonometric.hpp"

namespace varietas {

// Implicit equations of the reachable workspace, by elimination.
//
// Everything so far has treated the pose as a constant: a target goes in, an
// ideal in the joint variables comes out. Putting the pose coordinates into the
// ring instead turns the same residuals into a description of the graph of the
// forward kinematics map, and eliminating the joint variables projects that
// graph onto the pose coordinates. By the Closure Theorem the result is the
// ideal of the Zariski closure of the image — the equations every reachable
// point satisfies.
//
// The word closure carries real weight and the library should not pretend
// otherwise. The reachable set of an arm is semialgebraic: an annulus, a shell,
// something with a boundary, and boundaries are cut out by inequalities that no
// ideal can express. Zariski closure is the smallest algebraic set containing
// it, so when the workspace is full dimensional in the pose coordinates the
// closure is everything and the elimination ideal says only what was true
// identically — that a planar arm stays in its plane, for instance. The
// construction is informative exactly when the map is not dominant, that is,
// when the arm has fewer degrees of freedom than the space it moves in, and
// then it is very informative: a two-joint arm with perpendicular axes traces a
// torus, and elimination returns the quartic that cuts it out.

// The enlarged ring: the Rabinowitsch variable first, then the joint variables,
// then the pose coordinates, with a block order that eliminates everything but
// the pose in a single pass.
//
// Saturating and then eliminating are both eliminations, and performing them
// separately means computing a Gröbner basis of the saturation in full — an
// object nobody wants — and then computing another one to project it.
// Eliminating y and the joint variables together asks Buchberger for the only
// thing actually needed.
//
// The saving is smaller than that argument suggests: on the torus below it is
// about a tenth, seventy-eight seconds against seventy. The cost lies in the
// parameterisation rather than in the arrangement of the passes, which is why
// this path is no longer the one workspace_relations takes — see the
// trigonometric construction further down, where the same torus eliminates in
// sixty milliseconds and returns the identical quartic. It is kept because the
// agreement of the two is a test rather than an assumption, and because the
// contrast is the clearest statement of what each formulation costs.
template <std::size_t N>
using workspace_order = block_order<N + 1, grevlex, grevlex>;

template <class Coeff, std::size_t N>
using workspace_polynomial = polynomial<Coeff, N + 4, workspace_order<N>>;

// Variable layout of the enlarged ring.
template <std::size_t N>
struct workspace_layout {
  static constexpr std::size_t inverse_denominator = 0;
  static constexpr std::size_t first_joint = 1;
  static constexpr std::size_t first_pose = N + 1;
  static constexpr std::size_t eliminated_block = N + 1;
};

// The residuals of the position equations with the pose left symbolic: the tool
// coordinate minus the pose variable, denominators cleared. Their zero set is
// the graph of the forward kinematics map, away from the loci the substitution
// invented.
template <std::size_t N, class Coeff>
std::vector<workspace_polynomial<Coeff, N>> workspace_residuals(
    const chain<Coeff>& robot) {
  using enlarged = workspace_polynomial<Coeff, N>;
  using order = workspace_order<N>;
  using layout = workspace_layout<N>;

  const auto map = rational_forward_kinematics<N, grevlex>(robot);
  const auto denominator =
      embed_polynomial<order, N + 4>(map.denominator(), layout::first_joint);

  std::vector<enlarged> residuals;
  residuals.reserve(3);
  for (std::size_t i = 0; i < 3; ++i) {
    const auto translation =
        embed_polynomial<order, N + 4>(map.translation(i), layout::first_joint);
    residuals.push_back(translation -
                        denominator * enlarged::variable(layout::first_pose + i));
  }
  return residuals;
}

// Generators of the ideal of the Zariski closure of the reachable set of tool
// positions, in the three pose variables alone.
//
// The joint variables are saturated away first for the reason they always are:
// clearing the half-angle denominators attached the loci t_i = ±i to the graph,
// and projecting those would contribute components to the answer that no
// configuration of the arm reaches.
template <std::size_t N, class Coeff>
std::vector<polynomial<Coeff, 3, grevlex>> workspace_relations_half_angle(
    const chain<Coeff>& robot, buchberger_statistics* statistics = nullptr) {
  using enlarged = workspace_polynomial<Coeff, N>;
  using order = workspace_order<N>;
  using layout = workspace_layout<N>;

  std::vector<enlarged> generators = workspace_residuals<N>(robot);

  // 1 - y·h, with h the product of the half-angle denominators: the same
  // Rabinowitsch generator saturation would introduce, written here directly so
  // that its variable is eliminated in the same pass as the joint variables.
  const auto multipliers = half_angle_multipliers<N, grevlex>(robot);
  if (!multipliers.empty()) {
    polynomial<Coeff, N, grevlex> product =
        polynomial<Coeff, N, grevlex>::constant(coefficient_traits<Coeff>::one());
    for (const auto& h : multipliers) {
      product = product * h;
    }
    const enlarged lifted =
        embed_polynomial<order, N + 4>(product, layout::first_joint);
    generators.push_back(enlarged::constant(coefficient_traits<Coeff>::one()) -
                         enlarged::variable(layout::inverse_denominator) * lifted);
  }

  const auto graph = groebner_basis(generators, statistics);

  // The elements free of the eliminated block generate the elimination ideal,
  // and are already a Gröbner basis of it under the trailing block's order.
  const auto eliminated = eliminated_generators(graph, layout::eliminated_block);

  std::vector<polynomial<Coeff, 3, grevlex>> relations;
  relations.reserve(eliminated.size());
  for (const auto& g : eliminated) {
    relations.push_back(project_polynomial<grevlex, 3>(g, layout::first_pose));
  }
  reduce_groebner_basis(relations);
  return relations;
}

// The same construction in the trigonometric formulation, which is the one the
// library uses.
//
// The enlarged ring holds the joint variables first and the three pose
// coordinates after, with a block order eliminating the former. There is no
// Rabinowitsch variable and no saturation, because there is no denominator to
// have cleared: the circle relations enter as ordinary generators and the
// spurious components the half-angle substitution has to be cleaned of are
// never created. What is eliminated is exactly the parameter block, and the
// Closure Theorem applies verbatim.
//
// V is the size of the joint block — two variables per revolute joint and one
// per prismatic, which trigonometric_variable_count computes from the chain.
template <std::size_t V>
using trigonometric_workspace_order = block_order<V, grevlex, grevlex>;

template <class Coeff, std::size_t V>
using trigonometric_workspace_polynomial =
    polynomial<Coeff, V + 3, trigonometric_workspace_order<V>>;

template <std::size_t V>
struct trigonometric_workspace_layout {
  static constexpr std::size_t first_joint = 0;
  static constexpr std::size_t first_pose = V;
  static constexpr std::size_t eliminated_block = V;
};

// The residuals with the pose left symbolic. Nothing is cleared, so each is
// simply the tool coordinate minus the pose variable, and their zero set
// together with the circle relations is the graph of the forward kinematics map
// over the actual parameter space of the joints.
template <std::size_t V, class Coeff>
std::vector<trigonometric_workspace_polynomial<Coeff, V>>
trigonometric_workspace_generators(const chain<Coeff>& robot) {
  using enlarged = trigonometric_workspace_polynomial<Coeff, V>;
  using order = trigonometric_workspace_order<V>;
  using layout = trigonometric_workspace_layout<V>;

  const auto map = trigonometric_forward_kinematics<V, grevlex>(robot);

  std::vector<enlarged> generators;
  for (const auto& relation : circle_relations<V, grevlex>(robot)) {
    generators.push_back(embed_polynomial<order, V + 3>(relation, layout::first_joint));
  }
  for (std::size_t i = 0; i < 3; ++i) {
    const auto translation =
        embed_polynomial<order, V + 3>(map.translation(i), layout::first_joint);
    generators.push_back(translation - enlarged::variable(layout::first_pose + i));
  }
  return generators;
}

// Generators of the ideal of the Zariski closure of the reachable set of tool
// positions, in the three pose variables alone.
template <std::size_t V, class Coeff>
std::vector<polynomial<Coeff, 3, grevlex>> workspace_relations(
    const chain<Coeff>& robot, buchberger_statistics* statistics = nullptr) {
  using layout = trigonometric_workspace_layout<V>;

  const auto graph = groebner_basis(trigonometric_workspace_generators<V>(robot), statistics);
  const auto eliminated = eliminated_generators(graph, layout::eliminated_block);

  std::vector<polynomial<Coeff, 3, grevlex>> relations;
  relations.reserve(eliminated.size());
  for (const auto& g : eliminated) {
    relations.push_back(project_polynomial<grevlex, 3>(g, layout::first_pose));
  }
  reduce_groebner_basis(relations);
  return relations;
}

// True when the arm's tool positions fill their space, so that the Zariski
// closure is everything and the elimination ideal is zero. Reported rather than
// left for the caller to infer from an empty vector, because it is the case in
// which the construction has nothing to say and saying so is the honest answer.
template <class Coeff>
bool workspace_is_dense(const std::vector<polynomial<Coeff, 3, grevlex>>& relations) {
  return relations.empty();
}

}  // namespace varietas

#endif
