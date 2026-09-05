#ifndef VARIETAS_IK_EMIT_DECOUPLED_HPP
#define VARIETAS_IK_EMIT_DECOUPLED_HPP

#include <cstddef>
#include <sstream>
#include <string>

#include "varietas/codegen/emit.hpp"
#include "varietas/core/config.hpp"
#include "varietas/ik/decoupled_ik.hpp"

// A header for an arm solved by sweeping its first joint out.
//
// The reduced problem is emitted the way any solved system is. What this adds
// is the other half of the answer: the base angle, which is an arctangent of
// the target rather than an eigenvalue, and the pairing of each reduced
// solution with the two turns of the plane that reach the same point. Without
// it the generated header is only most of a solver, and the caller is left to
// rediscover a decomposition that was already made here.
namespace varietas {
namespace ik {

// The generated wrapper returns angles, not the half-angle variables the
// reduced solver returns.
//
// This is a deliberate break with the convention elsewhere. The variables of
// the ring are t = tan(q/2), and for the swept joint that substitution is
// exactly the wrong representation: the second family of solutions turns the
// plane half a turn further, so it produces base angles near pi routinely, and
// tan(q/2) is unbounded there. Returning the angle costs one arctangent per
// joint and is defined everywhere.
namespace detail {

inline std::string coordinate_field(std::size_t coordinate) {
  static const char* const names[3] = {"0", "1", "2"};
  return coordinate < 3 ? names[coordinate] : "0";
}

}  // namespace detail

// The wrapper, as text to be handed to emit_options::epilogue.
template <std::size_t N>
std::string decoupled_epilogue(const decoupled_solution<N>& solution,
                               const std::string& reduced_name,
                               const std::string& wrapper_name) {
  VARIETAS_ASSERT(solution.ok());

  const std::size_t reduced_unknowns = N - 1;
  const std::size_t reduced_dimension = solution.reduced.dimension();
  const std::string a = detail::coordinate_field(solution.frame.radial);
  const std::string b = detail::coordinate_field(solution.frame.swept);
  const std::string k = detail::coordinate_field(solution.frame.axis);
  const std::string sign = solution.frame.reversed ? "-" : "";

  std::ostringstream out;
  out << "// The whole arm: the reduced problem above, with its swept joint put back.\n"
      << "//\n"
      << "// Turning the base sweeps a plane around its axis, so a target is reached by\n"
      << "// turning that plane to face it and then solving in the plane, and also by\n"
      << "// turning half a revolution further and reaching backwards, which is why every\n"
      << "// solution of the reduced problem yields two configurations of the arm.\n"
      << "//\n"
      << "// Unlike the reduced solver above, this returns joint ANGLES in radians. The\n"
      << "// base angle routinely lands near pi, where tan(q/2) is unbounded, so the\n"
      << "// half-angle variable is not a usable output for it.\n"
      << "struct " << wrapper_name << " {\n"
      << "  static constexpr std::size_t num_joints = " << N << ";\n"
      << "  static constexpr std::size_t max_configurations = " << 2 * reduced_dimension << ";\n"
      << "\n"
      << "  using status = " << reduced_name << "::status;\n"
      << "\n"
      << "  // Writes up to `capacity` configurations into `out`, num_joints doubles each,\n"
      << "  // the swept joint first. Returns how many were written, or -1 if the target\n"
      << "  // lies on the locus the reduced basis does not describe.\n"
      << "  static int solve(const double* target, double* out, int capacity,\n"
      << "                   status* state = nullptr) {\n"
      << "    constexpr double half_turn = 3.14159265358979323846;\n"
      << "\n"
      << "    const double radius = std::hypot(target[" << a << "], target[" << b << "]);\n"
      << "    const double heading = std::atan2(target[" << b << "], target[" << a << "]);\n"
      << "\n"
      << "    // Facing the target, and reversed half a turn away with the radius\n"
      << "    // measured the other way. Both are genuine postures of the arm.\n"
      << "    const double families[2][2] = {\n"
      << "      { radius, " << sign << "heading },\n"
      << "      { -radius, " << sign << "(heading + half_turn) },\n"
      << "    };\n"
      << "\n"
      << "    int written = 0;\n"
      << "    bool described = false;\n"
      << "    for (std::size_t f = 0; f < 2; ++f) {\n"
      << "      const double pose[2] = { families[f][0], target[" << k << "] };\n"
      << "      double reduced_out[" << reduced_dimension * reduced_unknowns << "];\n"
      << "      status reduced_state{};\n"
      << "      const int found = " << reduced_name << "::solve(pose, reduced_out, "
      << reduced_dimension << ", &reduced_state);\n"
      << "      if (found < 0) { continue; }  // off the chart for this family only\n"
      << "      described = true;\n"
      << "      for (int s = 0; s < found; ++s) {\n"
      << "        if (written >= capacity) { break; }\n"
      << "        double* row = out + static_cast<std::size_t>(written) * num_joints;\n"
      << "        row[0] = families[f][1];\n"
      << "        for (std::size_t i = 0; i < " << reduced_unknowns << "; ++i) {\n"
      << "          // q = 2 arctan(t), the inverse of the half-angle substitution.\n"
      << "          row[1 + i] = 2.0 * std::atan(\n"
      << "              reduced_out[static_cast<std::size_t>(s) * " << reduced_unknowns
      << " + i]);\n"
      << "        }\n"
      << "        ++written;\n"
      << "      }\n"
      << "    }\n"
      << "\n"
      << "    // Neither turn of the plane was describable, so the target is on the\n"
      << "    // locus rather than merely out of reach. Out of reach returns zero.\n"
      << "    if (!described) {\n"
      << "      if (state != nullptr) { *state = status::bad_pose; }\n"
      << "      return -1;\n"
      << "    }\n"
      << "    if (state != nullptr) { *state = status::ok; }\n"
      << "    return written;\n"
      << "  }\n"
      << "};\n";
  return out.str();
}

// The complete header: the reduced solver and the wrapper that completes it.
template <std::size_t N>
std::string emit_decoupled(const decoupled_solution<N>& solution,
                           codegen::emit_options options = {}) {
  VARIETAS_ASSERT(solution.ok());
  // The wrapper calls the reduced solver, so the reduced solver has to exist.
  VARIETAS_ASSERT(options.runtime == codegen::runtime_kind::eigen);

  const std::string wrapper_name = options.name;
  const std::string reduced_name = options.name + "_reduced";

  options.epilogue = decoupled_epilogue(solution, reduced_name, wrapper_name);
  options.name = reduced_name;
  // The guard is defaulted from the struct's name, which has just changed; let
  // it be derived from the name the caller actually asked for instead.
  if (options.guard.empty()) {
    codegen::emit_options for_guard = options;
    for_guard.name = wrapper_name;
    options.guard = codegen::detail::default_guard(for_guard);
  }
  return codegen::emit(solution.reduced, options);
}

}  // namespace ik
}  // namespace varietas

#endif
