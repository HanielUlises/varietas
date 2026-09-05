// The cost of a parametric Grobner basis, over Q(p) and over F_p(p).
//
// Measures the same computation twice with only the coefficient field changed,
// to separate the cost of arithmetic on rational coefficients from the number
// and size of the polynomials the algorithm produces.
//
// Build (from the repository root):
//   g++ -std=c++17 -O2 tools/experiments/field_cost.cpp -o field_cost \
//     -Ivarietas_core/include -Ivarietas_codegen/include \
//     -Ivarietas_kinematics/include -Itools/experiments \
//     -I/usr/include/eigen3 -lgmpxx -lgmp

#include <array>
#include <chrono>
#include <cstddef>
#include <cstdio>
#include <string>
#include <vector>

#include "modular_field.hpp"

#include "varietas/codegen/rational_function.hpp"
#include "varietas/core/order/grevlex.hpp"
#include "varietas/core/quotient/quotient_basis.hpp"
#include "varietas/kinematics/chain.hpp"
#include "varietas/kinematics/rationalize.hpp"

using varietas::grevlex;

namespace {

// n as an element of whatever field is in play, without naming its type.
template <class Coeff>
Coeff scalar(int n) {
  using traits = varietas::coefficient_traits<Coeff>;
  Coeff result = traits::zero();
  const Coeff one = traits::one();
  for (int i = 0; i < (n < 0 ? -n : n); ++i) {
    result = result + one;
  }
  return n < 0 ? traits::negate(result) : result;
}

template <class Coeff>
varietas::rigid_transform<Coeff> along_x(int length) {
  using traits = varietas::coefficient_traits<Coeff>;
  return varietas::rigid_transform<Coeff>::translation_only(
      varietas::vector3<Coeff>(scalar<Coeff>(length), traits::zero(), traits::zero()));
}

// The planar two-link arm: two joints turning about z, unit links.
template <class Coeff>
varietas::chain<Coeff> planar_2r() {
  varietas::chain<Coeff> robot("planar_2r");
  robot.add_joint(varietas::revolute_joint<Coeff>("q1", varietas::vector3<Coeff>::unit(2),
                                                  varietas::rigid_transform<Coeff>::identity()));
  robot.add_joint(varietas::revolute_joint<Coeff>("q2", varietas::vector3<Coeff>::unit(2),
                                                  along_x<Coeff>(1)));
  robot.set_tool(along_x<Coeff>(1));
  return robot;
}

// The anthropomorphic arm: a base yawing about z, a shoulder and an elbow
// pitching about y, unit links.
template <class Coeff>
varietas::chain<Coeff> anthropomorphic_3r() {
  varietas::chain<Coeff> robot("anthropomorphic_3r");
  robot.add_joint(varietas::revolute_joint<Coeff>("q1", varietas::vector3<Coeff>::unit(2),
                                                  varietas::rigid_transform<Coeff>::identity()));
  robot.add_joint(varietas::revolute_joint<Coeff>("q2", varietas::vector3<Coeff>::unit(1),
                                                  varietas::rigid_transform<Coeff>::identity()));
  robot.add_joint(varietas::revolute_joint<Coeff>("q3", varietas::vector3<Coeff>::unit(1),
                                                  along_x<Coeff>(1)));
  robot.set_tool(along_x<Coeff>(1));
  return robot;
}

// The anthropomorphic arm with its base joint removed: the reduced problem the
// decoupling produces, in two joints and two parameters.
template <class Coeff>
varietas::chain<Coeff> shoulder_and_elbow() {
  varietas::chain<Coeff> robot("shoulder_and_elbow");
  robot.add_joint(varietas::revolute_joint<Coeff>("q2", varietas::vector3<Coeff>::unit(1),
                                                  varietas::rigid_transform<Coeff>::identity()));
  robot.add_joint(varietas::revolute_joint<Coeff>("q3", varietas::vector3<Coeff>::unit(1),
                                                  along_x<Coeff>(1)));
  robot.set_tool(along_x<Coeff>(1));
  return robot;
}

struct measurement {
  bool completed = false;
  double forward_map_seconds = 0.0;
  double residual_seconds = 0.0;
  double seconds = 0.0;
  std::size_t basis_size = 0;
  std::size_t quotient_dimension = 0;
  bool zero_dimensional = false;
  std::size_t largest_basis_polynomial = 0;  // terms
};

// The parametric position problem, posed and solved over the given field.
template <std::size_t N, std::size_t P, class Field>
measurement solve(const varietas::chain<Field>& robot,
                  const std::array<std::size_t, P>& coordinates) {
  using poly = varietas::polynomial<Field, N, grevlex>;

  // Timed in three phases, because the parametric pipeline does substantial
  // coefficient arithmetic before Buchberger is entered at all: the forward map
  // is a product of transforms whose entries are elements of the coefficient
  // field, and every one of those products is normalised.
  const auto phase0 = std::chrono::steady_clock::now();
  const auto map = varietas::rational_forward_kinematics<N, grevlex>(robot);
  const poly denominator = map.denominator();
  const auto phase1 = std::chrono::steady_clock::now();
  std::fprintf(stderr, "  [forward map done in %.1f s]\n",
               std::chrono::duration<double>(phase1 - phase0).count());
  std::fflush(stderr);

  std::vector<poly> residuals;
  residuals.reserve(P);
  for (std::size_t k = 0; k < P; ++k) {
    residuals.push_back(map.translation(coordinates[k]) -
                        denominator * poly::constant(Field::parameter(k)));
  }
  const auto phase2 = std::chrono::steady_clock::now();
  std::fprintf(stderr, "  [residuals done in %.1f s]\n",
               std::chrono::duration<double>(phase2 - phase1).count());
  std::fflush(stderr);

  const auto start = phase2;
  const auto basis = varietas::kinematic_ideal_generators<N, grevlex>(robot, residuals);
  const auto quotient = varietas::standard_monomials(basis);
  const auto stop = std::chrono::steady_clock::now();

  measurement m;
  m.completed = true;
  m.forward_map_seconds = std::chrono::duration<double>(phase1 - phase0).count();
  m.residual_seconds = std::chrono::duration<double>(phase2 - phase1).count();
  m.seconds = std::chrono::duration<double>(stop - start).count();
  m.basis_size = basis.size();
  m.zero_dimensional = quotient.is_zero_dimensional;
  m.quotient_dimension = quotient.dimension();
  for (const auto& g : basis) {
    if (g.size() > m.largest_basis_polynomial) {
      m.largest_basis_polynomial = g.size();
    }
  }
  return m;
}

void report(const char* system, const char* field, const measurement& m) {
  std::printf("%-28s %-12s %8.3f %8.3f %10.3f %7zu %8zu %8zu\n", system, field,
              m.forward_map_seconds, m.residual_seconds, m.seconds, m.basis_size,
              m.quotient_dimension, m.largest_basis_polynomial);
  std::fflush(stdout);
}

}  // namespace

int main(int argc, char** argv) {
  const std::string which = argc > 1 ? argv[1] : "all";

  std::printf("%-28s %-12s %8s %8s %10s %7s %8s %8s\n", "system", "field", "fwd map",
              "residual", "groebner", "basis", "dim_k A", "maxterm");
  std::fflush(stdout);

  if (which == "all" || which == "planar") {
    report("planar 2R (N=2, P=2)", "F_p(x,y)",
           solve<2, 2, varietas::modular_function<2>>(planar_2r<varietas::modular_function<2>>(),
                                                      {0, 1}));
    report("planar 2R (N=2, P=2)", "Q(x,y)",
           solve<2, 2, varietas::rational_function<2>>(planar_2r<varietas::rational_function<2>>(),
                                                       {0, 1}));
  }

  if (which == "all" || which == "reduced") {
    report("reduced 2R (N=2, P=2)", "F_p(r,z)",
           solve<2, 2, varietas::modular_function<2>>(
               shoulder_and_elbow<varietas::modular_function<2>>(), {0, 2}));
    report("reduced 2R (N=2, P=2)", "Q(r,z)",
           solve<2, 2, varietas::rational_function<2>>(
               shoulder_and_elbow<varietas::rational_function<2>>(), {0, 2}));
  }

  if (which == "all" || which == "spatial-fp") {
    report("anthropomorphic 3R (N=3, P=3)", "F_p(x,y,z)",
           solve<3, 3, varietas::modular_function<3>>(
               anthropomorphic_3r<varietas::modular_function<3>>(), {0, 1, 2}));
  }

  if (which == "spatial-q") {
    report("anthropomorphic 3R (N=3, P=3)", "Q(x,y,z)",
           solve<3, 3, varietas::rational_function<3>>(
               anthropomorphic_3r<varietas::rational_function<3>>(), {0, 1, 2}));
  }

  return 0;
}
