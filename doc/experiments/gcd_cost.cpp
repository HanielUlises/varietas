// Cost of one subresultant gcd on trivariate polynomials, over F_p and over Q.
//
// Isolates the two contributions to the cost of the cancellation that
// rational_function performs after every coefficient operation: the size of the
// coefficients, which the choice of field controls, and the number of terms in
// the polynomials, which it does not.
//
// Build (from the repository root):
//   g++ -std=c++17 -O2 tools/experiments/gcd_cost.cpp -o gcd_cost \
//     -Ivarietas_core/include -Ivarietas_codegen/include -Itools/experiments \
//     -I/usr/include/eigen3 -lgmpxx -lgmp
#include <chrono>
#include <cstdio>
#include <random>
#include <array>
#include <vector>

#include "modular_field.hpp"

#include "varietas/codegen/rational.hpp"
#include "varietas/core/gcd.hpp"
#include "varietas/core/order/grevlex.hpp"

using varietas::grevlex;

template <class Coeff>
varietas::polynomial<Coeff, 3, grevlex> dense(int degree, std::mt19937& rng) {
  using poly = varietas::polynomial<Coeff, 3, grevlex>;
  using mon = varietas::monomial<3>;
  std::vector<typename poly::term> terms;
  std::uniform_int_distribution<int> coeff(1, 50);
  for (int i = 0; i <= degree; ++i)
    for (int j = 0; i + j <= degree; ++j)
      for (int k = 0; i + j + k <= degree; ++k)
        terms.push_back({mon(std::array<mon::exponent_type, 3>{
                             (mon::exponent_type)i, (mon::exponent_type)j, (mon::exponent_type)k}),
                         Coeff(coeff(rng))});
  return poly(std::move(terms));
}

template <class Coeff>
void bench(const char* field, int degree) {
  std::mt19937 rng(7);
  const auto g = dense<Coeff>(degree / 2, rng);
  const auto a = dense<Coeff>(degree / 2, rng) * g;
  const auto b = dense<Coeff>(degree / 2, rng) * g;

  const auto t0 = std::chrono::steady_clock::now();
  const auto result = varietas::polynomial_gcd(a, b);
  const double sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  std::printf("%-10s deg %2d   inputs %5zu / %5zu terms   gcd %5zu terms   %9.3f s\n", field,
              degree, a.size(), b.size(), result.size(), sec);
  std::fflush(stdout);
}

int main() {
  for (int d : {4, 6, 8, 10}) {
    bench<varietas::modular>("F_p", d);
    bench<varietas::rational>("Q", d);
  }
  return 0;
}
