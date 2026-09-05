// Cost of one subresultant gcd on dense trivariate polynomials.
//
// Isolates the two contributions to the cost of the cancellation that
// rational_function performs after every coefficient operation: the size of the
// coefficients, which the choice of coefficient field controls, and the number
// of terms in the operands, which it does not.
//
// Two dense polynomials of equal total degree are formed with a planted common
// factor, so the gcd is nontrivial and the remainder sequence runs to its full
// length rather than terminating early on coprime inputs. The same inputs, up
// to reduction of the coefficients, are used over both fields.
//
// Build (from the repository root):
//   g++ -std=c++17 -O2 doc/experiments/gcd_cost.cpp -o gcd_cost \
//     -Ivarietas_core/include -Ivarietas_codegen/include -Idoc/experiments \
//     -I/usr/include/eigen3 -lgmpxx -lgmp
//
// Output is one row per (field, degree) in a form a plotting script can read,
// followed by a least-squares fit of log t against log m.

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <random>
#include <vector>

#include "modular_field.hpp"

#include "varietas/codegen/rational.hpp"
#include "varietas/core/gcd.hpp"
#include "varietas/core/order/grevlex.hpp"

using varietas::grevlex;

namespace {

// Every monomial of total degree at most `degree` in three variables, with
// coefficients drawn from a fixed sequence so that the two fields see the same
// integers.
template <class Coeff>
varietas::polynomial<Coeff, 3, grevlex> dense(int degree, std::mt19937& rng) {
  using poly = varietas::polynomial<Coeff, 3, grevlex>;
  using mon = varietas::monomial<3>;

  std::vector<typename poly::term> terms;
  std::uniform_int_distribution<int> coefficient(1, 50);
  for (int i = 0; i <= degree; ++i) {
    for (int j = 0; i + j <= degree; ++j) {
      for (int k = 0; i + j + k <= degree; ++k) {
        terms.push_back({mon(std::array<mon::exponent_type, 3>{
                             static_cast<mon::exponent_type>(i),
                             static_cast<mon::exponent_type>(j),
                             static_cast<mon::exponent_type>(k)}),
                         Coeff(coefficient(rng))});
      }
    }
  }
  return poly(std::move(terms));
}

struct sample {
  double seconds = 0.0;
  std::size_t operand_terms = 0;
  std::size_t gcd_terms = 0;
};

template <class Coeff>
sample measure(int degree, int repeats) {
  std::vector<double> times;
  sample s;
  for (int r = 0; r < repeats; ++r) {
    std::mt19937 rng(7 + r);
    const auto factor = dense<Coeff>(degree / 2, rng);
    const auto a = dense<Coeff>(degree / 2, rng) * factor;
    const auto b = dense<Coeff>(degree / 2, rng) * factor;

    const auto start = std::chrono::steady_clock::now();
    const auto result = varietas::polynomial_gcd(a, b);
    const auto stop = std::chrono::steady_clock::now();

    times.push_back(std::chrono::duration<double>(stop - start).count());
    s.operand_terms = a.size();
    s.gcd_terms = result.size();
  }
  std::sort(times.begin(), times.end());
  s.seconds = times[times.size() / 2];
  return s;
}

// Least squares on log t = alpha log m + beta, returning alpha.
double fitted_exponent(const std::vector<sample>& samples) {
  double sx = 0.0, sy = 0.0, sxx = 0.0, sxy = 0.0;
  const double n = static_cast<double>(samples.size());
  for (const auto& s : samples) {
    const double x = std::log(static_cast<double>(s.operand_terms));
    const double y = std::log(s.seconds);
    sx += x;
    sy += y;
    sxx += x * x;
    sxy += x * y;
  }
  return (n * sxy - sx * sy) / (n * sxx - sx * sx);
}

template <class Coeff>
std::vector<sample> sweep(const char* field, int lowest, int highest) {
  std::vector<sample> samples;
  for (int degree = lowest; degree <= highest; degree += 2) {
    // Cheap points are repeated and the median taken; expensive points are
    // measured once, since the quantity of interest spans four orders of
    // magnitude and is not sensitive to run-to-run variation.
    const int repeats = degree <= 8 ? 5 : 1;
    const sample s = measure<Coeff>(degree, repeats);
    samples.push_back(s);
    std::printf("%-6s %3d %7zu %7zu %12.6f\n", field, degree, s.operand_terms, s.gcd_terms,
                s.seconds);
    std::fflush(stdout);
  }
  return samples;
}

}  // namespace

int main() {
  std::printf("%-6s %3s %7s %7s %12s\n", "field", "deg", "terms", "gcd", "seconds");

  // The prime field is carried further, being cheaper; the rational field stops
  // where a single measurement would begin to dominate the run.
  const auto modular = sweep<varietas::modular>("Fp", 4, 16);
  const auto rational = sweep<varietas::rational>("Q", 4, 12);

  std::printf("\nfitted exponent of terms, Fp: %.2f\n", fitted_exponent(modular));
  std::printf("fitted exponent of terms, Q:  %.2f\n", fitted_exponent(rational));
  return 0;
}
