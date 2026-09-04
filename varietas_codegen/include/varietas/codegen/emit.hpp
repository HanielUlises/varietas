#ifndef VARIETAS_CODEGEN_EMIT_HPP
#define VARIETAS_CODEGEN_EMIT_HPP

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <sstream>
#include <string>
#include <vector>

#include "varietas/codegen/parametric_solution.hpp"
#include "varietas/codegen/rational_function.hpp"
#include "varietas/core/config.hpp"
#include "varietas/core/order/order_id.hpp"

namespace varietas {
namespace codegen {

// How much of the runtime the generated header carries.
enum class runtime_kind {
  // Action matrices and the coordinates of the variables, and nothing else.
  // The header includes <cstddef> and <cstdint> and can be dropped into a
  // project that has never heard of Eigen or of this library.
  matrices_only,
  // The above plus solve(), which eigendecomposes the matrices and returns the
  // points. This pulls Eigen into the generated header, which is the price of
  // the header being a solver rather than a description of one.
  eigen,
};

// What the generated header is called and where it lives.
struct emit_options {
  std::string name = "solver";           // the struct's name
  std::string name_space = "varietas_generated";
  std::string guard;                     // defaulted from the namespace and name
  std::string source_note;               // free text recorded in the banner
  runtime_kind runtime = runtime_kind::matrices_only;
};

namespace detail {

// A rational written so that the compiler does the rounding, once.
//
// The alternative is a decimal literal, which rounds here and then again when
// the compiler parses it, and which cannot represent 1/3 at all. Writing the
// quotient of two exact integer literals hands the compiler the same rational
// the offline computation held and lets it produce the correctly rounded double
// in one step. Values too large for the mantissa to distinguish are the reason
// the fallback exists, not a reason to prefer decimals generally.
inline std::string emit_rational(const rational& q) {
  const mpz_class& numerator = q.get_num();
  const mpz_class& denominator = q.get_den();

  const bool fits = numerator.fits_slong_p() && denominator.fits_slong_p();
  if (fits) {
    std::ostringstream out;
    out << numerator.get_str();
    if (denominator != 1) {
      out << ".0 / " << denominator.get_str() << ".0";
    } else {
      out << ".0";
    }
    return out.str();
  }

  // Beyond long the literals would not survive the compiler's own parsing, so
  // the value is rounded here instead, at the full precision of a double.
  std::ostringstream out;
  out.precision(17);
  out << std::scientific << q.get_d();
  return out.str();
}

// One monomial in the pose parameters, as a product. Powers are written out
// rather than sent through pow: the exponents here are small, and a product of
// two or three multiplications is both faster and exactly what the compiler can
// fold when the pose is known.
template <std::size_t P>
std::string emit_monomial(const monomial<P>& m, const std::string& pose) {
  std::string text;
  for (std::size_t i = 0; i < P; ++i) {
    for (unsigned e = 0; e < static_cast<unsigned>(m[i]); ++e) {
      if (!text.empty()) {
        text += " * ";
      }
      text += pose + "[" + std::to_string(i) + "]";
    }
  }
  return text;
}

// A polynomial in the pose parameters. The constant polynomial 1 is the common
// case for a denominator and is written as a bare literal rather than as an
// empty product.
template <std::size_t P>
std::string emit_polynomial(const polynomial<rational, P, typename rational_function<P>::parameter_order>& p,
                            const std::string& pose) {
  if (p.is_zero()) {
    return "0.0";
  }

  std::ostringstream out;
  bool first = true;
  for (const auto& t : p.terms()) {
    if (!first) {
      out << " + ";
    }
    first = false;

    const std::string factors = emit_monomial<P>(t.mon, pose);
    const std::string coefficient = emit_rational(t.coeff);
    if (factors.empty()) {
      out << "(" << coefficient << ")";
    } else if (t.coeff == 1) {
      out << factors;
    } else {
      out << "(" << coefficient << ") * " << factors;
    }
  }
  return out.str();
}

template <std::size_t P>
std::string emit_function(const rational_function<P>& f, const std::string& pose) {
  const std::string numerator = emit_polynomial<P>(f.numerator(), pose);
  if (f.denominator().degree() == 0 && f.denominator().leading_coefficient() == 1) {
    return numerator;
  }
  return "(" + numerator + ") / (" + emit_polynomial<P>(f.denominator(), pose) + ")";
}

inline std::string default_guard(const emit_options& options) {
  std::string guard = options.name_space + "_" + options.name + "_GENERATED_HPP";
  for (char& c : guard) {
    c = (c >= 'a' && c <= 'z') ? static_cast<char>(c - 'a' + 'A')
                               : ((c >= 'A' && c <= 'Z') || (c >= '0' && c <= '9') ? c : '_');
  }
  return guard;
}

inline std::string name_or(const std::vector<std::string>& names, std::size_t i,
                           const std::string& prefix) {
  return i < names.size() && !names[i].empty() ? names[i] : prefix + std::to_string(i);
}

}  // namespace detail

// Writes a self-contained C++ header for a solved parametric system.
//
// Self-contained is meant strictly: the result includes <cstddef> and <cstdint>
// and nothing else, names no varietas type, and can be dropped into a project
// that has never heard of this library. That is the point of emitting at all.
// Everything expensive — the Gröbner basis, the quotient, the cancellation of
// the rational functions — has already happened by the time this is called, and
// what is left in the header is arithmetic on doubles.
//
// What the header does not do is solve the eigenproblem. The action matrices
// are what the offline half can compute once and for all; recovering the
// configurations from them is an eigendecomposition, which needs a linear
// algebra library and a choice about how to handle complex spectra, and baking
// either of those into generated source would make it much less portable than
// the arithmetic it is otherwise made of.
template <std::size_t N, std::size_t P>
std::string emit(const parametric_solution<N, P>& solution, const emit_options& options = {}) {
  VARIETAS_ASSERT(solution.is_well_formed());

  const std::string guard = options.guard.empty() ? detail::default_guard(options) : options.guard;
  const std::size_t d = solution.dimension();
  const std::string pose = "pose";

  std::ostringstream out;

  out << "// Generated by varietas. Do not edit.\n"
      << "//\n"
      << "// The action matrices of a zero-dimensional ideal whose coefficients are\n"
      << "// rational functions of the pose, evaluated at a pose supplied at runtime.\n"
      << "// The eigenvalues of these matrices are the coordinates of the solutions:\n"
      << "// for each unknown, the spectrum of its matrix is the set of values that\n"
      << "// unknown takes over the solution set, and a simultaneous eigenvector ties\n"
      << "// the coordinates of one solution together.\n";
  if (!options.source_note.empty()) {
    out << "//\n// " << options.source_note << "\n";
  }
  out << "\n#ifndef " << guard << "\n#define " << guard << "\n\n"
      << "#include <cstddef>\n#include <cstdint>\n";
  if (options.runtime == runtime_kind::eigen) {
    out << "#include <cmath>\n#include <complex>\n\n"
        << "#include <Eigen/Core>\n#include <Eigen/Eigenvalues>\n";
  }
  out << "\n"
      << "namespace " << options.name_space << " {\n\n"
      << "struct " << options.name << " {\n";

  out << "  // The monomial order the basis was computed under. A runtime that\n"
      << "  // assumes a different one is asking a different question, so this is\n"
      << "  // recorded rather than left to convention.\n"
      << "  static constexpr std::uint8_t order_id = "
      << static_cast<unsigned>(static_cast<std::uint8_t>(solution.order)) << ";\n"
      << "  static constexpr const char* order_name = \"" << to_string(solution.order) << "\";\n\n"
      << "  static constexpr std::size_t num_unknowns = " << N << ";\n"
      << "  static constexpr std::size_t num_parameters = " << P << ";\n\n"
      << "  // dim_k A: the number of solutions counted with multiplicity, and the\n"
      << "  // size of every action matrix below.\n"
      << "  static constexpr std::size_t dimension = " << d << ";\n\n"
      << "  // Where the monomial 1 sits in the standard basis. The eigenvalue\n"
      << "  // method divides by the eigenvector's component here.\n"
      << "  static constexpr std::size_t one_index = " << solution.one_index << ";\n\n";

  out << "  // The unknowns, in the order their matrices are indexed:\n";
  for (std::size_t i = 0; i < N; ++i) {
    out << "  //   " << i << ": " << detail::name_or(solution.unknown_names, i, "x") << "\n";
  }
  out << "  // The pose parameters, in the order they are read from `pose`:\n";
  for (std::size_t i = 0; i < P; ++i) {
    out << "  //   " << i << ": " << detail::name_or(solution.parameter_names, i, "p") << "\n";
  }
  out << "\n";

  out << "  // Writes the action matrix of unknown `variable` at `" << pose << "` into\n"
      << "  // `out`, column-major, dimension * dimension entries. Returns false if\n"
      << "  // `variable` is out of range or a denominator vanished at this pose,\n"
      << "  // which is the pose lying on the locus the parametric basis does not\n"
      << "  // describe; `out` is not meaningful in that case.\n"
      << "  static bool action_matrix(std::size_t variable, const double* " << pose
      << ", double* out) {\n"
      << "    switch (variable) {\n";

  for (std::size_t v = 0; v < N; ++v) {
    const auto& matrix = solution.action[v];
    out << "      case " << v << ": {  // "
        << detail::name_or(solution.unknown_names, v, "x") << "\n";

    // Denominators are checked before anything is written, so a pose on the
    // locus the parametric basis does not describe leaves `out` untouched
    // rather than half filled with infinities.
    //
    // The test is written inline rather than through a named temporary: the
    // same denominator commonly appears in many entries, and a temporary would
    // need a distinct name for each. Distinct denominators are deduplicated for
    // the same reason, so the guard costs one comparison per pole and not one
    // per entry.
    std::vector<std::string> guards;
    for (std::size_t j = 0; j < d; ++j) {
      for (std::size_t i = 0; i < d; ++i) {
        const auto& f = matrix(i, j);
        if (f.is_zero() || f.denominator().degree() == 0) {
          continue;  // a nonzero constant denominator cannot vanish
        }
        std::string g = detail::emit_polynomial<P>(f.denominator(), pose);
        if (std::find(guards.begin(), guards.end(), g) == guards.end()) {
          guards.push_back(std::move(g));
        }
      }
    }
    for (const std::string& g : guards) {
      out << "        if ((" << g << ") == 0.0) { return false; }\n";
    }

    for (std::size_t j = 0; j < d; ++j) {
      for (std::size_t i = 0; i < d; ++i) {
        const auto& f = matrix(i, j);
        out << "        out[" << (j * d + i) << "] = ";
        if (f.is_zero()) {
          out << "0.0";
        } else {
          out << detail::emit_function<P>(f, pose);
        }
        out << ";  // (" << i << ", " << j << ")\n";
      }
    }
    out << "        return true;\n      }\n";
  }

  out << "      default:\n        return false;\n"
      << "    }\n  }\n";

  // The coordinates of each variable's normal form, which is what turns an
  // eigenvector into a point. Emitted whatever the runtime, because a caller
  // writing its own eigensolver needs them exactly as much as ours does.
  out << "\n  // Writes the coordinates of the normal form of unknown `variable`\n"
      << "  // in the standard basis into `out`, `dimension` entries. A variable is\n"
      << "  // usually not a standard monomial itself, so its value at a point is\n"
      << "  // this combination evaluated there, not a coordinate read off.\n"
      << "  static bool variable_coordinates(std::size_t variable, const double* " << pose
      << ", double* out) {\n"
      << "    switch (variable) {\n";
  for (std::size_t v = 0; v < N; ++v) {
    const auto& row = solution.variable_coordinates[v];
    out << "      case " << v << ": {  // "
        << detail::name_or(solution.unknown_names, v, "x") << "\n";

    std::vector<std::string> guards;
    for (const auto& f : row) {
      if (f.is_zero() || f.denominator().degree() == 0) {
        continue;
      }
      std::string g = detail::emit_polynomial<P>(f.denominator(), pose);
      if (std::find(guards.begin(), guards.end(), g) == guards.end()) {
        guards.push_back(std::move(g));
      }
    }
    for (const std::string& g : guards) {
      out << "        if ((" << g << ") == 0.0) { return false; }\n";
    }
    for (std::size_t m = 0; m < d; ++m) {
      out << "        out[" << m << "] = "
          << (row[m].is_zero() ? std::string("0.0") : detail::emit_function<P>(row[m], pose))
          << ";\n";
    }
    out << "        return true;\n      }\n";
  }
  out << "      default:\n        return false;\n    }\n  }\n";

  if (options.runtime == runtime_kind::eigen) {
    out << R"CODE(
  // Solve, by the eigenvalue method of Stetter and Moller.
  //
  // The action matrices commute, so they share eigenvectors. A left eigenvector
  // of the multiplication operator is the evaluation functional at a point of
  // the variety: it sends each standard monomial m to m(p). Left eigenvectors of
  // M are right eigenvectors of its transpose, which is what is decomposed here.
  //
  // A single variable will not do as the operator. If two distinct solutions
  // give it the same value the eigenvalue is repeated, the eigenvector is an
  // arbitrary element of a plane, and it carries no point. A linear combination
  // with the fixed coefficients below separates the solutions except on a proper
  // algebraic subset of the coefficient space; the coefficients are fixed rather
  // than random so that two runs on the same pose return the same points in the
  // same order.
  //
  // Having the functional, x_i(p) is the normal form of x_i evaluated at it,
  // divided by what it gives to 1 — the division that makes the functional an
  // evaluation rather than a multiple of one.
  enum class status {
    ok = 0,
    bad_pose,             // a denominator vanished: this pose is off the chart
    eigensolver_failed,   // Eigen did not converge
    deficient,            // an eigenvector carried no point; the set is not certified
  };

  // Writes the real solutions, `num_unknowns` doubles each, into `out`, which
  // must hold at least `capacity * num_unknowns`. Returns the number written,
  // or -1 on failure with `*state` set to why.
  //
  // Real, because a joint angle that is complex is not a configuration a
  // manipulator can be commanded to. Complex points are found and then
  // discarded; a caller wanting all of them can lift this loop out.
  static int solve(const double* pose, double* out, int capacity, status* state = nullptr,
                   double tolerance = 1e-8) {
    const auto fail = [&](status s) {
      if (state != nullptr) { *state = s; }
      return -1;
    };
    if (state != nullptr) { *state = status::ok; }

    using matrix_type = Eigen::Matrix<double, dimension, dimension>;

    // The separating form, as a combination of the action matrices.
    matrix_type separating = matrix_type::Zero();
    double single[dimension * dimension];
    for (std::size_t i = 0; i < num_unknowns; ++i) {
      if (!action_matrix(i, pose, single)) { return fail(status::bad_pose); }
      const double c = 1.0 + 0.37 * static_cast<double>(i)
                     + 0.11 * static_cast<double>(i * i);
      separating += c * Eigen::Map<const matrix_type>(single);
    }

    double coordinates[num_unknowns][dimension];
    for (std::size_t i = 0; i < num_unknowns; ++i) {
      if (!variable_coordinates(i, pose, coordinates[i])) { return fail(status::bad_pose); }
    }

    Eigen::EigenSolver<matrix_type> solver(separating.transpose(), true);
    if (solver.info() != Eigen::Success) { return fail(status::eigensolver_failed); }

    const auto vectors = solver.eigenvectors();
    int written = 0;
    bool discarded = false;
    for (Eigen::Index k = 0; k < vectors.cols(); ++k) {
      const auto v = vectors.col(k);
      const std::complex<double> scale = v(static_cast<Eigen::Index>(one_index));
      if (std::abs(scale) <= tolerance) { discarded = true; continue; }

      std::complex<double> point[num_unknowns];
      bool real = true;
      for (std::size_t i = 0; i < num_unknowns; ++i) {
        std::complex<double> value{0.0, 0.0};
        for (std::size_t m = 0; m < dimension; ++m) {
          if (coordinates[i][m] != 0.0) {
            value += coordinates[i][m] * v(static_cast<Eigen::Index>(m));
          }
        }
        point[i] = value / scale;
        if (std::abs(point[i].imag()) > tolerance) { real = false; }
      }
      if (!real) { continue; }
      if (written >= capacity) { break; }
      for (std::size_t i = 0; i < num_unknowns; ++i) {
        out[static_cast<std::size_t>(written) * num_unknowns + i] = point[i].real();
      }
      ++written;
    }

    if (discarded && state != nullptr) { *state = status::deficient; }
    return written;
  }
)CODE";
  }

  out << "};\n\n"
      << "}  // namespace " << options.name_space << "\n\n"
      << "#endif\n";

  return out.str();
}

}  // namespace codegen
}  // namespace varietas

#endif
