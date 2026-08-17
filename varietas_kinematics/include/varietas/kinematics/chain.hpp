#ifndef VARIETAS_KINEMATICS_CHAIN_HPP
#define VARIETAS_KINEMATICS_CHAIN_HPP

#include <cstddef>
#include <limits>
#include <string>
#include <utility>
#include <vector>

#include "varietas/core/config.hpp"
#include "varietas/kinematics/rigid_transform.hpp"

namespace varietas {

// The description a robot is reduced to before any algebra happens.
//
// Nothing here knows about URDF, and that is the point: the rationalisation,
// the code emitter, the workspace implicitization and the singular locus all
// consume this structure, so a URDF front end is later a translation into it
// and not a dependency of the pipeline. Hand-written chains and DH tables enter
// by the same door, which is also what makes the algebra testable without a
// parser.
//
// A chain is a list of joints in order from base to tip. Each joint carries the
// fixed placement of its frame relative to the previous joint's frame, together
// with the axis it moves along or about; the tool transform closes the chain at
// the end-effector frame. Fixed joints are kept rather than folded into their
// neighbours so that the description remains in correspondence with the robot
// model it came from; folding is a transformation on the chain, offered below.

enum class joint_type {
  fixed,
  revolute,
  prismatic,
};

// A revolute joint contributes one variable t = tan(q/2) to the polynomial
// ring, a prismatic joint contributes its displacement directly, and a fixed
// joint contributes none. Limits are carried as doubles: they are inequalities,
// they play no part in any ideal, and they are consulted only when solutions of
// the variety are filtered down to the reachable ones.
template <class Coeff>
struct joint {
  std::string name;
  joint_type type = joint_type::fixed;
  vector3<Coeff> axis = vector3<Coeff>::unit(2);
  rigid_transform<Coeff> origin = rigid_transform<Coeff>::identity();

  bool has_limits = false;
  double lower = -std::numeric_limits<double>::infinity();
  double upper = std::numeric_limits<double>::infinity();

  bool is_actuated() const noexcept { return type != joint_type::fixed; }
};

enum class chain_status {
  ok,
  // The rotation of a joint origin is not orthogonal, so the placement is not
  // an element of SE(3) and every pose computed through it is wrong.
  origin_not_orthogonal,
  // Orthogonal but of determinant -1: a reflection, not a rotation.
  origin_improper,
  // A moving joint whose axis is the zero vector: it names no motion.
  axis_degenerate,
  // A moving joint whose axis is not of unit length. The half-angle
  // substitution presumes a unit axis, since Rodrigues' formula is a rotation
  // only then; scaling the axis silently rescales the joint variable.
  axis_not_unit,
  // Limits given the wrong way round.
  limits_inverted,
  // No moving joint at all, hence no variables and nothing to solve for.
  no_actuated_joints,
};

// A named failure with the joint it belongs to, in preference to a bare bool:
// the caller needs to know which of a dozen joints is malformed, and how.
struct chain_diagnostic {
  chain_status status = chain_status::ok;
  std::size_t joint_index = 0;
  double defect = 0.0;

  bool ok() const noexcept { return status == chain_status::ok; }
  explicit operator bool() const noexcept { return ok(); }
};

inline const char* to_string(chain_status status) noexcept {
  switch (status) {
    case chain_status::ok:
      return "ok";
    case chain_status::origin_not_orthogonal:
      return "joint origin rotation is not orthogonal";
    case chain_status::origin_improper:
      return "joint origin rotation is a reflection, not a rotation";
    case chain_status::axis_degenerate:
      return "actuated joint has a zero axis";
    case chain_status::axis_not_unit:
      return "actuated joint axis is not of unit length";
    case chain_status::limits_inverted:
      return "joint limits are inverted";
    case chain_status::no_actuated_joints:
      return "chain has no actuated joints";
  }
  return "unknown";
}

template <class Coeff>
class chain {
 public:
  using traits = coefficient_traits<Coeff>;
  using joint_type_alias = joint<Coeff>;

  static constexpr std::size_t npos = static_cast<std::size_t>(-1);

  chain() = default;

  explicit chain(std::string name) : name_(std::move(name)) {}

  const std::string& name() const noexcept { return name_; }
  void set_name(std::string name) { name_ = std::move(name); }

  chain& add_joint(joint<Coeff> j) {
    joints_.push_back(std::move(j));
    return *this;
  }

  const std::vector<joint<Coeff>>& joints() const noexcept { return joints_; }
  std::size_t size() const noexcept { return joints_.size(); }

  const rigid_transform<Coeff>& tool() const noexcept { return tool_; }
  void set_tool(rigid_transform<Coeff> t) { tool_ = std::move(t); }

  // The number of variables the polynomial ring will have, which for a chain of
  // revolute joints is the N that the rest of the library is templated on.
  std::size_t degrees_of_freedom() const noexcept {
    std::size_t count = 0;
    for (const joint<Coeff>& j : joints_) {
      count += j.is_actuated() ? 1 : 0;
    }
    return count;
  }

  // Position of a joint among the actuated ones, that is, the index of its
  // variable in the ring; npos for a fixed joint. The two directions are kept
  // explicit because generated code indexes by variable while a robot model
  // indexes by joint, and conflating them is how a solution ends up applied to
  // the wrong axis.
  std::size_t variable_of_joint(std::size_t joint_index) const {
    VARIETAS_ASSERT(joint_index < joints_.size());
    if (!joints_[joint_index].is_actuated()) {
      return npos;
    }
    std::size_t variable = 0;
    for (std::size_t i = 0; i < joint_index; ++i) {
      variable += joints_[i].is_actuated() ? 1 : 0;
    }
    return variable;
  }

  std::size_t joint_of_variable(std::size_t variable_index) const {
    std::size_t seen = 0;
    for (std::size_t i = 0; i < joints_.size(); ++i) {
      if (!joints_[i].is_actuated()) {
        continue;
      }
      if (seen == variable_index) {
        return i;
      }
      ++seen;
    }
    return npos;
  }

  // Checks that every joint describes an element of SE(3) and a motion the
  // rationalisation can represent. The tolerance is ignored over an exact
  // field, where the only admissible defect is zero: an origin that is merely
  // nearly orthogonal has already lost the exactness the offline pipeline
  // exists to provide, and accepting it would make the Gröbner basis an exact
  // statement about the wrong robot.
  chain_diagnostic validate(double tolerance = 1e-12) const {
    const double allowed = traits::is_exact ? 0.0 : tolerance;

    for (std::size_t i = 0; i < joints_.size(); ++i) {
      const joint<Coeff>& j = joints_[i];

      const double defect = orthogonality_defect(j.origin.rotation());
      if (defect > allowed) {
        return {chain_status::origin_not_orthogonal, i, defect};
      }

      const double det = traits::to_double(j.origin.rotation().determinant());
      if (det < 0.0) {
        return {chain_status::origin_improper, i, det};
      }

      if (j.is_actuated()) {
        if (j.axis.is_zero()) {
          return {chain_status::axis_degenerate, i, 0.0};
        }
        double norm_defect =
            traits::to_double(j.axis.squared_norm() - traits::one());
        if (norm_defect < 0.0) {
          norm_defect = -norm_defect;
        }
        if (norm_defect > allowed) {
          return {chain_status::axis_not_unit, i, norm_defect};
        }
      }

      if (j.has_limits && j.lower > j.upper) {
        return {chain_status::limits_inverted, i, j.lower - j.upper};
      }
    }

    if (degrees_of_freedom() == 0) {
      return {chain_status::no_actuated_joints, npos, 0.0};
    }
    return {};
  }

  // The fixed transform accumulated from the tip of the chain, that is, the
  // product of the origins of the trailing fixed joints with the tool frame.
  // The rationalisation needs it as the constant factor that multiplies the
  // last moving joint's contribution.
  rigid_transform<Coeff> trailing_transform() const {
    rigid_transform<Coeff> t = tool_;
    for (std::size_t i = joints_.size(); i-- > 0;) {
      if (joints_[i].is_actuated()) {
        break;
      }
      t = joints_[i].origin * t;
    }
    return t;
  }

  // The same chain with runs of fixed joints multiplied into the origin of the
  // joint that follows them, and any trailing run folded into the tool frame.
  // The two chains have the same forward kinematics by associativity in SE(3),
  // and the folded one is what the rationalisation should be handed: every
  // fixed joint left in place is a matrix product repeated in every generator.
  chain fold_fixed_joints() const {
    chain folded(name_);
    rigid_transform<Coeff> pending = rigid_transform<Coeff>::identity();

    for (const joint<Coeff>& j : joints_) {
      if (!j.is_actuated()) {
        pending = pending * j.origin;
        continue;
      }
      joint<Coeff> moved = j;
      moved.origin = pending * j.origin;
      folded.add_joint(std::move(moved));
      pending = rigid_transform<Coeff>::identity();
    }

    folded.set_tool(pending * tool_);
    return folded;
  }

 private:
  std::string name_;
  std::vector<joint<Coeff>> joints_;
  rigid_transform<Coeff> tool_ = rigid_transform<Coeff>::identity();
};

// Convenience constructors for the joints a URDF actually contains, where the
// axis is a coordinate direction and the origin a pure translation.

template <class Coeff>
joint<Coeff> revolute_joint(std::string name, const vector3<Coeff>& axis,
                            const rigid_transform<Coeff>& origin) {
  joint<Coeff> j;
  j.name = std::move(name);
  j.type = joint_type::revolute;
  j.axis = axis;
  j.origin = origin;
  return j;
}

template <class Coeff>
joint<Coeff> prismatic_joint(std::string name, const vector3<Coeff>& axis,
                             const rigid_transform<Coeff>& origin) {
  joint<Coeff> j;
  j.name = std::move(name);
  j.type = joint_type::prismatic;
  j.axis = axis;
  j.origin = origin;
  return j;
}

template <class Coeff>
joint<Coeff> fixed_joint(std::string name, const rigid_transform<Coeff>& origin) {
  joint<Coeff> j;
  j.name = std::move(name);
  j.type = joint_type::fixed;
  j.origin = origin;
  return j;
}

}  // namespace varietas

#endif
