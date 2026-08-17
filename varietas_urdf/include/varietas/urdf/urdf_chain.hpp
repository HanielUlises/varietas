#ifndef VARIETAS_URDF_URDF_CHAIN_HPP
#define VARIETAS_URDF_URDF_CHAIN_HPP

#include <algorithm>
#include <string>
#include <vector>

#include <urdf_model/model.h>

#include "varietas/codegen/rational.hpp"
#include "varietas/kinematics/chain.hpp"
#include "varietas/urdf/rational_approximation.hpp"

namespace varietas {
namespace urdf_import {

// The front end: a URDF, which describes geometry in decimal, becomes a chain
// over the exact field, or is refused with a reason.
//
// Refusal is the interesting half. Every rotation in the file is recovered as
// an exact rational quaternion and the recovery is only accepted if it lies
// within a stated tolerance of what the file said; a joint whose placement is
// genuinely oblique — not a multiple of a right angle, not a fraction of small
// denominator — cannot be represented exactly and is reported rather than
// rounded. The tolerance is thus a statement about the model, not about the
// arithmetic: it says how far from the file the exact robot is allowed to be,
// and the offline pipeline is exact about that robot from there on.

enum class import_status {
  ok,
  tip_link_not_found,
  root_link_not_found,
  // The walk from the tip reached a link with no parent joint before reaching
  // the requested root: the two are not on a common chain.
  root_not_on_chain,
  // Floating and planar joints move in more than one degree of freedom and have
  // no single axis to substitute for.
  unsupported_joint_type,
  // The joint axis is not a rational unit vector of bounded denominator.
  axis_not_rational_unit,
  // The exact rotation nearest the stated one is further away than allowed.
  rotation_deviation_exceeded,
  translation_deviation_exceeded,
  no_actuated_joints,
  // The imported chain failed varietas::chain::validate, which should not
  // happen and means the importer, not the file, is at fault.
  chain_invalid,
};

inline const char* to_string(import_status status) noexcept {
  switch (status) {
    case import_status::ok:
      return "ok";
    case import_status::tip_link_not_found:
      return "tip link not found in the model";
    case import_status::root_link_not_found:
      return "root link not found in the model";
    case import_status::root_not_on_chain:
      return "root link is not an ancestor of the tip link";
    case import_status::unsupported_joint_type:
      return "joint type has no single axis of motion";
    case import_status::axis_not_rational_unit:
      return "joint axis is not a rational unit vector";
    case import_status::rotation_deviation_exceeded:
      return "no exact rotation lies within tolerance of the stated one";
    case import_status::translation_deviation_exceeded:
      return "no exact translation lies within tolerance of the stated one";
    case import_status::no_actuated_joints:
      return "chain has no actuated joints";
    case import_status::chain_invalid:
      return "imported chain failed validation";
  }
  return "unknown";
}

struct import_options {
  // Denominators up to a million recover every decimal literal a URDF is
  // likely to hold, and every multiple of a right angle exactly.
  long max_denominator = 1000000;
  // How far the exact robot may sit from the stated one. The default is loose
  // enough to absorb a truncated pi and tight enough that a real misalignment
  // is refused.
  double rotation_tolerance = 1e-9;
  double translation_tolerance = 1e-9;
};

// What the recovery cost at one joint, kept for every joint so that a model can
// be audited rather than merely accepted.
struct joint_recovery {
  std::string name;
  double rotation_deviation = 0.0;
  double translation_deviation = 0.0;
  // The recovered placement is the stated one to the last bit of a double, so
  // the joint has not moved at all. A joint whose angles were written as
  // truncated decimals is not of this kind: it moves, by the truncation.
  bool unmoved = false;
};

struct import_report {
  import_status status = import_status::ok;
  std::string detail;
  std::vector<joint_recovery> joints;
  double max_rotation_deviation = 0.0;
  double max_translation_deviation = 0.0;

  bool ok() const noexcept { return status == import_status::ok; }
  explicit operator bool() const noexcept { return ok(); }

  // True when the exact chain is the stated one to the last bit throughout, so
  // that exactness cost no motion at all. Usually false, and harmlessly so: a
  // model that writes pi/2 as a decimal is moved onto the right angle it meant,
  // by an amount max_rotation_deviation reports.
  bool indistinguishable_from_file() const noexcept {
    if (!ok()) {
      return false;
    }
    for (const joint_recovery& recovery : joints) {
      if (!recovery.unmoved) {
        return false;
      }
    }
    return true;
  }
};

// Joints from the tip up to the root, reversed into base-to-tip order.
inline import_report collect_joints(const ::urdf::ModelInterface& model,
                                    const std::string& root_link,
                                    const std::string& tip_link,
                                    std::vector<::urdf::JointConstSharedPtr>& out) {
  import_report report;
  out.clear();

  ::urdf::LinkConstSharedPtr link = model.getLink(tip_link);
  if (!link) {
    report.status = import_status::tip_link_not_found;
    report.detail = tip_link;
    return report;
  }
  if (!model.getLink(root_link)) {
    report.status = import_status::root_link_not_found;
    report.detail = root_link;
    return report;
  }

  while (link && link->name != root_link) {
    const ::urdf::JointConstSharedPtr joint = link->parent_joint;
    if (!joint) {
      report.status = import_status::root_not_on_chain;
      report.detail = link->name;
      return report;
    }
    out.push_back(joint);
    link = model.getLink(joint->parent_link_name);
  }

  std::reverse(out.begin(), out.end());
  return report;
}

// Builds the chain, filling the report as it goes. On failure the chain is left
// untouched and the report names the joint responsible.
inline import_report chain_from_model(const ::urdf::ModelInterface& model,
                                      const std::string& root_link,
                                      const std::string& tip_link,
                                      chain<rational>& out,
                                      const import_options& options = {}) {
  std::vector<::urdf::JointConstSharedPtr> joints;
  import_report report = collect_joints(model, root_link, tip_link, joints);
  if (!report.ok()) {
    return report;
  }

  chain<rational> built(model.getName());

  for (const ::urdf::JointConstSharedPtr& source : joints) {
    joint<rational> target;
    target.name = source->name;

    switch (source->type) {
      case ::urdf::Joint::REVOLUTE:
      case ::urdf::Joint::CONTINUOUS:
        target.type = joint_type::revolute;
        break;
      case ::urdf::Joint::PRISMATIC:
        target.type = joint_type::prismatic;
        break;
      case ::urdf::Joint::FIXED:
        target.type = joint_type::fixed;
        break;
      default:
        report.status = import_status::unsupported_joint_type;
        report.detail = source->name;
        return report;
    }

    // Rotation, recovered projectively from the quaternion urdfdom has already
    // built out of the file's roll-pitch-yaw triple.
    const ::urdf::Pose& pose = source->parent_to_joint_origin_transform;
    double qx = 0.0, qy = 0.0, qz = 0.0, qw = 1.0;
    pose.rotation.getQuaternion(qx, qy, qz, qw);
    const rotation_snap rotation = snap_rotation(qx, qy, qz, qw, options.max_denominator);
    if (rotation.deviation_radians > options.rotation_tolerance) {
      report.status = import_status::rotation_deviation_exceeded;
      report.detail = source->name;
      return report;
    }

    const scalar_snap tx = snap_scalar(pose.position.x, options.max_denominator);
    const scalar_snap ty = snap_scalar(pose.position.y, options.max_denominator);
    const scalar_snap tz = snap_scalar(pose.position.z, options.max_denominator);
    const double translation_deviation =
        std::max({tx.deviation, ty.deviation, tz.deviation});
    if (translation_deviation > options.translation_tolerance) {
      report.status = import_status::translation_deviation_exceeded;
      report.detail = source->name;
      return report;
    }

    target.origin = rigid_transform<rational>(
        rotation.exact(), vector3<rational>(tx.value, ty.value, tz.value));

    if (target.is_actuated()) {
      // The axis has to be a rational unit vector, which for the coordinate
      // directions of a real URDF it is. Anything else would rescale the joint
      // variable under the half-angle substitution, so it is refused rather
      // than renormalised.
      const vector3<rational> axis(rationalize(source->axis.x, options.max_denominator),
                                   rationalize(source->axis.y, options.max_denominator),
                                   rationalize(source->axis.z, options.max_denominator));
      if (axis.squared_norm() != rational(1)) {
        report.status = import_status::axis_not_rational_unit;
        report.detail = source->name;
        return report;
      }
      target.axis = axis;

      if (source->type == ::urdf::Joint::REVOLUTE ||
          source->type == ::urdf::Joint::PRISMATIC) {
        if (source->limits) {
          target.has_limits = true;
          target.lower = source->limits->lower;
          target.upper = source->limits->upper;
        }
      }
    }

    joint_recovery recovery;
    recovery.name = source->name;
    recovery.rotation_deviation = rotation.deviation_radians;
    recovery.translation_deviation = translation_deviation;
    recovery.unmoved =
        rotation.round_trips && tx.round_trips && ty.round_trips && tz.round_trips;
    report.joints.push_back(std::move(recovery));

    report.max_rotation_deviation =
        std::max(report.max_rotation_deviation, rotation.deviation_radians);
    report.max_translation_deviation =
        std::max(report.max_translation_deviation, translation_deviation);

    built.add_joint(std::move(target));
  }

  if (built.degrees_of_freedom() == 0) {
    report.status = import_status::no_actuated_joints;
    return report;
  }

  const chain_diagnostic diagnostic = built.validate();
  if (!diagnostic.ok()) {
    report.status = import_status::chain_invalid;
    report.detail = to_string(diagnostic.status);
    return report;
  }

  out = std::move(built);
  return report;
}

// The tip of a model that is a single chain, that is, the one link with no
// children. Convenience for the common case; a branched model needs the tip
// naming which branch is wanted.
inline std::string sole_tip_link(const ::urdf::ModelInterface& model) {
  std::vector<::urdf::LinkSharedPtr> links;
  model.getLinks(links);
  std::string tip;
  for (const ::urdf::LinkSharedPtr& link : links) {
    if (link->child_joints.empty()) {
      if (!tip.empty()) {
        return std::string();
      }
      tip = link->name;
    }
  }
  return tip;
}

}  // namespace urdf_import
}  // namespace varietas

#endif
