// Every configuration at once.
//
// The sweep demonstration answers the forward question: it drives the joints
// and shows that the tool goes where the exactly recovered chain says it goes.
// This one answers the question the library exists for. A target is moved
// through the workspace, the generated solver is asked for the configurations
// that reach it, and *all* of them are drawn together.
//
// What that shows is the thing a numerical solver cannot produce. An iteration
// started somewhere returns the branch it happened to fall into; the quotient
// dimension says how many there are, and the eigenvalue method returns them
// all, so the picture can be complete rather than representative. When the
// target leaves the workspace the arms disappear rather than lunging at a
// point they cannot reach, and as it approaches the boundary the elbow-up and
// elbow-down solutions converge on the straight arm in view.
//
// The header this includes was written by varietas during this build, from the
// URDF, by the path `urdf_codegen --decouple` takes. Nothing from the library
// is linked at solve time: branch_ik::solve is ordinary C++ over Eigen, which
// is the claim the emitter makes.

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <memory>
#include <string>
#include <vector>

#include <rclcpp/rclcpp.hpp>
#include <tf2_ros/static_transform_broadcaster.h>
#include <sensor_msgs/msg/joint_state.hpp>
#include <std_msgs/msg/color_rgba.hpp>
#include <visualization_msgs/msg/marker.hpp>
#include <visualization_msgs/msg/marker_array.hpp>

#include <urdf/model.h>

#include "varietas/ik/parametric_ik.hpp"
#include "varietas/kinematics/evaluate.hpp"
#include "varietas/urdf/urdf_chain.hpp"

#include "branch_ik.hpp"

using namespace std::chrono_literals;

namespace {

using solver = varietas_generated::branch_ik;

struct rgb {
  double r, g, b;
};

// One colour per branch, equally saturated on purpose: no branch is the answer
// and the others alternatives. They are all solutions, and drawing one boldly
// would say otherwise.
const rgb kBranchColours[4] = {
    {0.35, 0.70, 1.00},
    {1.00, 0.62, 0.22},
    {0.40, 0.90, 0.50},
    {0.92, 0.45, 0.86},
};

using vec3 = varietas::vector3<double>;

geometry_msgs::msg::Point point_of(const vec3& v) {
  geometry_msgs::msg::Point p;
  p.x = v[0];
  p.y = v[1];
  p.z = v[2];
  return p;
}

double distance(const vec3& a, const vec3& b) {
  const double dx = a[0] - b[0], dy = a[1] - b[1], dz = a[2] - b[2];
  return std::sqrt(dx * dx + dy * dy + dz * dz);
}

// A cylinder marker runs along its own z, so a limb drawn between two points
// needs the rotation carrying z onto the segment. The shortest arc will do:
// axis z x u, angle between them, written in half-angle form. The antiparallel
// case has no unique axis and any perpendicular one is correct.
void span(visualization_msgs::msg::Marker& m, const vec3& a, const vec3& b) {
  const double dx = b[0] - a[0], dy = b[1] - a[1], dz = b[2] - a[2];
  const double len = std::sqrt(dx * dx + dy * dy + dz * dz);
  m.scale.z = len;
  m.pose.position.x = 0.5 * (a[0] + b[0]);
  m.pose.position.y = 0.5 * (a[1] + b[1]);
  m.pose.position.z = 0.5 * (a[2] + b[2]);
  m.pose.orientation.x = m.pose.orientation.y = m.pose.orientation.z = 0.0;
  m.pose.orientation.w = 1.0;
  if (len < 1e-12) {
    return;
  }
  const double ux = dx / len, uy = dy / len, uz = dz / len;
  if (uz < -1.0 + 1e-9) {
    m.pose.orientation.x = 1.0;
    m.pose.orientation.w = 0.0;
    return;
  }
  const double n = std::sqrt(2.0 * (1.0 + uz));
  m.pose.orientation.x = -uy / n;
  m.pose.orientation.y = ux / n;
  m.pose.orientation.z = 0.0;
  m.pose.orientation.w = n / 2.0;
}

}  // namespace

class branch_node : public rclcpp::Node {
  static constexpr std::chrono::milliseconds kTick{33};
  static constexpr int kMaxBranches = static_cast<int>(solver::max_configurations);
  // Marker ids are branch * kLimbStride + segment, so a branch's limbs can be
  // cleared without touching another's.
  static constexpr int kLimbStride = 16;

 public:
  branch_node() : rclcpp::Node("varietas_branches") {
    const std::string urdf_path = declare_parameter<std::string>("urdf", "");
    period_ = declare_parameter<double>("period", 20.0);

    urdf::Model model;
    if (urdf_path.empty() || !model.initFile(urdf_path)) {
      RCLCPP_FATAL(get_logger(), "could not parse urdf '%s'", urdf_path.c_str());
      throw std::runtime_error("urdf");
    }

    const std::string tip = varietas::urdf_import::sole_tip_link(model);
    root_link_ = declare_parameter<std::string>("root_link", model.getRoot()->name);

    varietas::chain<varietas::rational> exact;
    const auto imported =
        varietas::urdf_import::chain_from_model(model, root_link_, tip, exact);
    if (!imported.ok()) {
      RCLCPP_FATAL(get_logger(), "import refused: %s (%s)",
                   varietas::urdf_import::to_string(imported.status),
                   imported.detail.c_str());
      throw std::runtime_error("import");
    }
    if (exact.degrees_of_freedom() != solver::num_joints) {
      RCLCPP_FATAL(get_logger(), "the generated solver is for %zu joints, this arm has %zu",
                   solver::num_joints, exact.degrees_of_freedom());
      throw std::runtime_error("joints");
    }

    robot_ = varietas::chain_cast<double>(exact);
    for (const auto& j : robot_.joints()) {
      if (j.is_actuated()) {
        names_.push_back(j.name);
      }
    }

    solved_name_ = exact.name();
    measure_reach();
    // Clear of the tallest posture rather than at a height chosen by eye, so
    // the label does not end up behind an arm on a model of another size.
    label_height_ = reach_ * 1.06;
    RCLCPP_INFO(get_logger(),
                "%s: %zu joints, reach %.3f m, up to %zu configurations per target",
                exact.name().c_str(), exact.degrees_of_freedom(), reach_,
                solver::max_configurations);

    joint_states_ = create_publisher<sensor_msgs::msg::JointState>("joint_states", 10);
    markers_ = create_publisher<visualization_msgs::msg::MarkerArray>("varietas_markers", 10);

    setup_refused_arm();

    timer_ = create_wall_timer(kTick, [this] { tick(); });
    report_timer_ = create_wall_timer(2s, [this] { report(); });
    start_ = now();
  }

 private:
  // The other half of the story: a real robot, and what varietas says about it.
  //
  // The verdict is computed rather than written down. The iiwa is imported by
  // the same front end, and the same parametric_position_ik the pipeline uses
  // is asked for a position solver over its seven joints. It refuses on the
  // counts, before any Grobner basis is attempted, and the reason it gives is
  // the reason shown beside the arm. If that ever stops being true this
  // demonstration says something else, which is the point of not hardcoding it.
  void setup_refused_arm() {
    const std::string path = declare_parameter<std::string>("refused_urdf", "");
    // Clear of the workspace rather than merely beside it. The solved arm
    // reaches about 2.4 m, so an offset of 1.9 m stood the iiwa inside the
    // set the four coloured arms sweep through, and they passed through it
    // several times a period. The default is now outside that reach, so the
    // two arms read as two robots on a bench rather than as one drawing
    // superimposed on another.
    refused_offset_ = declare_parameter<std::vector<double>>("refused_offset",
                                                             {0.0, -2.9, 0.0});
    if (path.empty()) {
      return;
    }

    urdf::Model model;
    if (!model.initFile(path)) {
      RCLCPP_WARN(get_logger(), "could not parse '%s'; the refused arm is omitted",
                  path.c_str());
      return;
    }
    const std::string tip = varietas::urdf_import::sole_tip_link(model);
    varietas::chain<varietas::rational> exact;
    const auto imported =
        varietas::urdf_import::chain_from_model(model, model.getRoot()->name, tip, exact);
    if (!imported.ok()) {
      RCLCPP_WARN(get_logger(), "refused arm import failed: %s",
                  varietas::urdf_import::to_string(imported.status));
      return;
    }

    refused_root_ = model.getRoot()->name;
    refused_ = varietas::chain_cast<double>(exact);
    for (const auto& j : refused_.joints()) {
      if (j.is_actuated()) {
        refused_names_.push_back(j.name);
      }
    }
    refused_joints_ = exact.degrees_of_freedom();
    refused_name_ = exact.name();

    // Seven joints against three position coordinates. The counting check
    // settles it: P < N cannot cut out a finite solution set, whatever the arm.
    const auto verdict = varietas::ik::parametric_position_ik<7, 3>(
        exact, std::array<std::size_t, 3>{0, 1, 2});
    refused_reason_ = varietas::ik::to_string(verdict.status);
    RCLCPP_INFO(get_logger(), "%s: %zu joints, recovered to %.1e rad; asked for a position "
                              "solver it answers: %s",
                exact.name().c_str(), refused_joints_, imported.max_rotation_deviation,
                refused_reason_.c_str());

    static_tf_ = std::make_shared<tf2_ros::StaticTransformBroadcaster>(this);
    geometry_msgs::msg::TransformStamped placement;
    placement.header.stamp = now();
    placement.header.frame_id = root_link_;
    placement.child_frame_id = refused_root_;
    placement.transform.translation.x = refused_offset_.size() > 0 ? refused_offset_[0] : 0.0;
    placement.transform.translation.y = refused_offset_.size() > 1 ? refused_offset_[1] : -2.9;
    placement.transform.translation.z = refused_offset_.size() > 2 ? refused_offset_[2] : 0.0;
    placement.transform.rotation.w = 1.0;
    static_tf_->sendTransform(placement);

    refused_states_ =
        create_publisher<sensor_msgs::msg::JointState>("refused/joint_states", 10);
  }

  // The refused arm is posed, not solved: nothing here asks it to reach
  // anything, because varietas has already said it cannot be asked. It moves so
  // that it reads as a robot rather than a prop.
  void publish_refused(const rclcpp::Time& stamp, double t) {
    if (refused_names_.empty() || !refused_states_) {
      return;
    }
    sensor_msgs::msg::JointState js;
    js.header.stamp = stamp;
    js.name = refused_names_;
    js.position.resize(refused_names_.size());
    for (std::size_t i = 0; i < refused_names_.size(); ++i) {
      const double f = 1.0 + 0.34 * static_cast<double>(i);
      js.position[i] = 0.85 * std::sin(2.0 * M_PI * f * t / (period_ * 1.7) + 0.6 * i);
    }
    refused_states_->publish(js);
  }

  // The furthest the tool actually gets, found by sampling rather than by
  // adding link lengths up: the sum is only an upper bound once a joint origin
  // carries a rotation, and the trajectory has to cross the boundary exactly,
  // not approximately, for the arms to vanish where the picture says they do.
  void measure_reach() {
    std::vector<double> q(robot_.degrees_of_freedom(), 0.0);
    double best = 0.0;
    const int steps = 21;
    for (int a = 0; a < steps; ++a) {
      for (int b = 0; b < steps; ++b) {
        for (int c = 0; c < steps; ++c) {
          q[0] = -M_PI + 2.0 * M_PI * a / (steps - 1);
          if (q.size() > 1) q[1] = -M_PI + 2.0 * M_PI * b / (steps - 1);
          if (q.size() > 2) q[2] = -M_PI + 2.0 * M_PI * c / (steps - 1);
          // By value. `translation()` hands back a reference into the
          // transform, and the transform here is a temporary: binding a
          // const reference to it does not extend its life, because what is
          // bound is the returned reference rather than the temporary
          // itself. The sampled reach came back 0.000 m that way, which the
          // trajectory below multiplies into a target that never leaves the
          // origin, so the arms neither moved nor vanished and every tick
          // reported the same four configurations.
          const vec3 p = varietas::forward_kinematics(robot_, q).translation();
          best = std::max(best, std::sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]));
        }
      }
    }
    reach_ = best;
  }

  // A closed path that deliberately leaves the workspace. The radius pulses
  // twice per revolution between about a third of the reach and a little past
  // it, so each revolution crosses the boundary four times.
  //
  // The overshoot is small on purpose. What is worth watching is the approach
  // rather than the vanishing: as the target nears the boundary the arm
  // straightens and the elbow-up and elbow-down solutions converge, so the
  // four arms close into two before they go.
  //
  // What the count does on the way out is the part worth having, and it is
  // not what was written here before. It steps four, two, zero rather than
  // four straight to zero, and two is not a measure-zero accident that the
  // sampling misses: it is held over a band of radii, about a quarter of all
  // ticks. The offset shoulder is the reason. The four postures are two
  // solutions of the reduced two-joint problem taken twice over, once with
  // the base facing the target and once with it turned half a revolution
  // away, and the offset puts the shoulder on the near side of the base axis
  // in the first family and on the far side in the second. The two families
  // therefore have different reach -- they differ by twice the offset -- and
  // between the two radii only the facing family arrives. So the annulus
  // where exactly two configurations exist is genuinely thick, the arms go in
  // pairs rather than all at once, and the trail is coloured by the count
  // below so that the two boundaries can be seen on the curve itself.
  //
  // The elevation is held positive. It swung through zero once, which sent
  // the target and its trail below the floor for part of every revolution:
  // the curve left the frame through the grid, and an arm reaching under the
  // bench it is bolted to is not a picture of anything. The band is set so
  // the lowest point of the path clears the plinth instead.
  std::array<double, 3> target_at(double t) const {
    const double u = 2.0 * M_PI * t / period_;
    const double radius = reach_ * (0.66 + 0.36 * std::sin(2.0 * u));
    const double elevation = 0.36 * std::sin(3.0 * u) + 0.40;
    const double c = std::cos(elevation);
    return {radius * c * std::cos(u), radius * c * std::sin(u), radius * std::sin(elevation)};
  }

  void tick() {
    const double t = (now() - start_).seconds();
    const auto target = target_at(t);

    solver::status state{};
    std::array<double, solver::max_configurations * solver::num_joints> out{};
    const int found = solver::solve(target.data(), out.data(), kMaxBranches, &state);

    last_count_ = found;
    last_state_ = state;
    ++ticks_;
    if (found < 0) {
      ++refusals_;
    } else {
      ++seen_[static_cast<std::size_t>(found)];
    }

    trail_.push_back({target, found});
    const std::size_t capacity =
        static_cast<std::size_t>(period_ * 1000.0 / static_cast<double>(kTick.count()));
    while (trail_.size() > capacity) {
      trail_.erase(trail_.begin());
    }

    publish(target, out, found, state);
    publish_refused(now(), t);
  }

  void publish(const std::array<double, 3>& target,
               const std::array<double, solver::max_configurations * solver::num_joints>& out,
               int found, solver::status /*state*/) {
    visualization_msgs::msg::MarkerArray array;
    const auto stamp = now();
    double worst = 0.0;

    for (int b = 0; b < kMaxBranches; ++b) {
      if (b >= found) {
        for (int seg = 0; seg < kLimbStride; ++seg) {
          array.markers.push_back(erased("limb", b * kLimbStride + seg, stamp));
        }
        array.markers.push_back(erased("joints", b, stamp));
        continue;
      }

      const std::vector<double> q(
          out.begin() + static_cast<std::size_t>(b) * solver::num_joints,
          out.begin() + static_cast<std::size_t>(b + 1) * solver::num_joints);
      const auto frames = varietas::link_frames(robot_, q);
      const auto reached = frames.back().translation();

      // The configuration is drawn from the same forward kinematics that checks
      // it, so a branch that does not reach the target would be visibly short
      // rather than quietly wrong.
      worst = std::max(worst, std::sqrt(std::pow(reached[0] - target[0], 2) +
                                        std::pow(reached[1] - target[1], 2) +
                                        std::pow(reached[2] - target[2], 2)));

      // Coincident frames carry no limb. A chain with joints at a shared origin
      // produces several of them, and drawing zero-length cylinders leaves
      // artefacts at the base.
      std::vector<vec3> joints;
      for (const auto& f : frames) {
        if (joints.empty() || distance(joints.back(), f.translation()) > 1e-9) {
          joints.push_back(f.translation());
        }
      }

      const rgb& colour = kBranchColours[b % 4];
      for (std::size_t seg = 0; seg + 1 < joints.size() && seg < kLimbStride; ++seg) {
        auto limb = base_marker("limb", b * kLimbStride + static_cast<int>(seg), stamp);
        limb.type = visualization_msgs::msg::Marker::CYLINDER;
        limb.scale.x = limb.scale.y = 0.055;
        limb.color.r = colour.r;
        limb.color.g = colour.g;
        limb.color.b = colour.b;
        limb.color.a = 0.93;
        span(limb, joints[seg], joints[seg + 1]);
        array.markers.push_back(limb);
      }
      for (std::size_t seg = joints.size() ? joints.size() - 1 : 0; seg < kLimbStride; ++seg) {
        array.markers.push_back(erased("limb", b * kLimbStride + static_cast<int>(seg), stamp));
      }

      auto knuckles = base_marker("joints", b, stamp);
      knuckles.type = visualization_msgs::msg::Marker::SPHERE_LIST;
      knuckles.scale.x = knuckles.scale.y = knuckles.scale.z = 0.085;
      knuckles.color.r = colour.r;
      knuckles.color.g = colour.g;
      knuckles.color.b = colour.b;
      knuckles.color.a = 1.0;
      for (const auto& j : joints) {
        knuckles.points.push_back(point_of(j));
      }
      array.markers.push_back(knuckles);
    }
    worst_residual_ = std::max(worst_residual_, worst);

    // A plinth at the origin and the axis the swept joint turns about. Without
    // them the four arms read as a cage floating over a grid: they share a
    // base, and the picture has to say where it is.
    auto plinth = base_marker("base", 0, stamp);
    plinth.type = visualization_msgs::msg::Marker::CYLINDER;
    plinth.scale.x = plinth.scale.y = 0.40;
    plinth.scale.z = 0.10;
    plinth.pose.position.z = 0.05;
    plinth.color.r = plinth.color.g = plinth.color.b = 0.38;
    plinth.color.a = 1.0;
    array.markers.push_back(plinth);

    auto axis = base_marker("base", 1, stamp);
    axis.type = visualization_msgs::msg::Marker::LINE_LIST;
    axis.scale.x = 0.006;
    axis.color.r = axis.color.g = axis.color.b = 0.55;
    axis.color.a = 0.5;
    geometry_msgs::msg::Point lo, hi;
    lo.x = lo.y = 0.0; lo.z = 0.0;
    hi = lo; hi.z = reach_ * 0.75;
    axis.points.push_back(lo);
    axis.points.push_back(hi);
    array.markers.push_back(axis);

    // The target. Its colour is the verdict: reached, out of reach, or refused
    // because the pose sits on the locus the parametric basis cannot describe.
    auto goal = base_marker("target", 0, stamp);
    goal.type = visualization_msgs::msg::Marker::SPHERE;
    goal.scale.x = goal.scale.y = goal.scale.z = 0.13;
    goal.pose.position.x = target[0];
    goal.pose.position.y = target[1];
    goal.pose.position.z = target[2];
    if (found < 0) {
      goal.color.r = 0.92; goal.color.g = 0.22; goal.color.b = 0.22;
    } else if (found == 0) {
      goal.color.r = 0.96; goal.color.g = 0.76; goal.color.b = 0.18;
    } else {
      goal.color.r = 0.97; goal.color.g = 0.97; goal.color.b = 0.97;
    }
    goal.color.a = 1.0;
    array.markers.push_back(goal);

    // The trail, coloured by how many configurations reached each point of it.
    //
    // A flat grey curve said only where the target had been, which the moving
    // sphere already says. Coloured by the count it carries the thing the
    // demonstration is about: the curve is one closed loop, and the places
    // where it changes colour are where the target crossed a boundary of the
    // reachable set. Those crossings are the same instants at which arms
    // appear and vanish, so the still frame explains the motion.
    //
    // It is a LINE_LIST rather than a LINE_STRIP because a strip takes one
    // colour for the whole curve; a list takes a colour per vertex, and each
    // segment is emitted as its two endpoints. The count belongs to the
    // segment rather than to either end, so both endpoints are given the
    // count of the later sample, which puts the change of colour exactly at
    // the tick where the count changed.
    auto path = base_marker("trail", 0, stamp);
    path.type = visualization_msgs::msg::Marker::LINE_LIST;
    path.scale.x = 0.012;
    path.color.a = 1.0;
    for (std::size_t i = 1; i < trail_.size(); ++i) {
      const std_msgs::msg::ColorRGBA c = trail_colour(trail_[i].found);
      path.points.push_back(point_of_array(trail_[i - 1].at));
      path.colors.push_back(c);
      path.points.push_back(point_of_array(trail_[i].at));
      path.colors.push_back(c);
    }
    array.markers.push_back(path);

    array.markers.push_back(verdict_label(found, stamp));
    if (!refused_reason_.empty()) {
      array.markers.push_back(refused_label(stamp));
    }

    markers_->publish(array);

    if (found > 0) {
      sensor_msgs::msg::JointState js;
      js.header.stamp = stamp;
      js.name = names_;
      js.position.assign(out.begin(), out.begin() + solver::num_joints);
      joint_states_->publish(js);
    }
  }

  // The same three verdicts the target sphere carries, so that a reader who
  // has learned the sphere's colours can read the curve without a key: white
  // where every configuration was found, amber where none was, and the count
  // between them where the target is reachable by one family of postures and
  // not the other.
  static std_msgs::msg::ColorRGBA trail_colour(int found) {
    std_msgs::msg::ColorRGBA c;
    c.a = 0.85f;
    if (found < 0) {
      c.r = 0.92f; c.g = 0.22f; c.b = 0.22f;
    } else if (found == 0) {
      c.r = 0.96f; c.g = 0.76f; c.b = 0.18f;
    } else if (found < kMaxBranches) {
      c.r = 0.55f; c.g = 0.80f; c.b = 0.95f;
    } else {
      c.r = 0.93f; c.g = 0.93f; c.b = 0.95f;
    }
    return c;
  }

  // Breaks a sentence onto lines of at most `columns` characters and sets the
  // spaces in underscores.
  //
  // Two constraints of the renderer, neither of them negotiable from here.
  // The refusal the library returns is a sentence rather than a token, and a
  // TEXT_VIEW_FACING marker lays a line out in world units: at a legible
  // height, seventy characters on one line is metres wide and runs off both
  // sides of the frame, so it has to be wrapped. And rviz_rendering's
  // MovableText advances a space by a width of its own choosing rather than
  // the font's, which comes out at some eight characters: a line of ordinary
  // prose arrives with its words scattered across the frame and unreadable.
  // The marker message has no field for it -- MovableText::setSpaceWidth is
  // not reachable through visualization_msgs -- so the text cannot contain a
  // space. Underscores are what is left, and they are not out of place here,
  // since both arms are named in the picture by their URDF identifiers, which
  // are written that way already.
  //
  // What matters is that the wrapping and the substitution are typographic
  // only: every word of the library's own sentence survives them in order, so
  // the label cannot come to say something the library does not.
  static std::string wrap(const std::string& text, std::size_t columns) {
    std::string out;
    std::size_t line = 0;
    for (std::size_t i = 0; i < text.size();) {
      std::size_t end = text.find(' ', i);
      if (end == std::string::npos) {
        end = text.size();
      }
      const std::size_t word = end - i;
      if (line != 0 && line + 1 + word > columns) {
        out += '\n';
        line = 0;
      } else if (line != 0) {
        out += ' ';
        ++line;
      }
      out.append(text, i, word);
      line += word;
      i = end + 1;
    }
    for (char& c : out) {
      if (c == ' ') {
        c = '_';
      }
    }
    return out;
  }

  static geometry_msgs::msg::Point point_of_array(const std::array<double, 3>& a) {
    geometry_msgs::msg::Point p;
    p.x = a[0];
    p.y = a[1];
    p.z = a[2];
    return p;
  }

  // What the arms are, said in the frame rather than in the log.
  //
  // The count is live because it is the reading: an arm that has just vanished
  // is indistinguishable in a still from an arm hidden behind another, and the
  // number settles it. The wording follows the verdict rather than the count
  // alone, since zero configurations found and a solver that declined to
  // answer are different statements and the picture should not merge them.
  visualization_msgs::msg::Marker verdict_label(int found,
                                                const rclcpp::Time& stamp) const {
    auto m = base_marker("label", 0, stamp);
    m.type = visualization_msgs::msg::Marker::TEXT_VIEW_FACING;
    m.pose.position.z = label_height_;
    m.scale.z = 0.17;
    m.color.r = m.color.g = m.color.b = 0.93;
    m.color.a = 0.95;
    std::string line = solved_name_ + "\n" + wrap(std::to_string(solver::num_joints) +
                                                   " joints solved exactly", 30) +
                       "\nconfigurations_drawn:";
    line += (found < 0) ? "declined"
                        : std::to_string(found) + "/" + std::to_string(kMaxBranches);
    m.text = line;
    return m;
  }

  // The refusal, beside the arm it is about. The string is the one
  // parametric_position_ik returned at start-up, not a paraphrase of it, so
  // this label cannot drift away from what the library actually says.
  visualization_msgs::msg::Marker refused_label(const rclcpp::Time& stamp) const {
    auto m = base_marker("refused", 0, stamp);
    m.type = visualization_msgs::msg::Marker::TEXT_VIEW_FACING;
    m.pose.position.x = refused_offset_.size() > 0 ? refused_offset_[0] : 0.0;
    m.pose.position.y = refused_offset_.size() > 1 ? refused_offset_[1] : -2.9;
    m.pose.position.z = 1.52;
    m.scale.z = 0.145;
    m.color.r = 0.96;
    m.color.g = 0.62;
    m.color.b = 0.25;
    m.color.a = 0.95;
    m.text = refused_name_ + "\n" +
             wrap(std::to_string(refused_joints_) + " joints refused:", 26) + "\n" +
             wrap(refused_reason_, 26);
    return m;
  }

  visualization_msgs::msg::Marker base_marker(const std::string& ns, int id,
                                              const rclcpp::Time& stamp) const {
    visualization_msgs::msg::Marker m;
    m.header.frame_id = root_link_;
    m.header.stamp = stamp;
    m.ns = ns;
    m.id = id;
    m.action = visualization_msgs::msg::Marker::ADD;
    m.pose.orientation.w = 1.0;
    return m;
  }

  visualization_msgs::msg::Marker erased(const std::string& ns, int id,
                                         const rclcpp::Time& stamp) const {
    auto m = base_marker(ns, id, stamp);
    m.action = visualization_msgs::msg::Marker::DELETE;
    return m;
  }

  // What the picture cannot show. That every drawn configuration actually
  // reaches the target is a number, and it belongs in the log at the size it
  // really is rather than in a marker nobody can measure by eye.
  void report() {
    std::string histogram;
    for (std::size_t k = 0; k < seen_.size(); ++k) {
      if (seen_[k] == 0) {
        continue;
      }
      histogram += (histogram.empty() ? "" : ", ") + std::to_string(k) + "x" +
                   std::to_string(seen_[k]);
    }
    RCLCPP_INFO(get_logger(),
                "%d now; worst tool error over %ld ticks %.3e m; counts seen [%s]; %ld refused",
                last_count_, ticks_, worst_residual_, histogram.c_str(), refusals_);
  }

  varietas::chain<double> robot_;
  std::vector<std::string> names_;
  std::string root_link_;
  double period_ = 20.0;
  double reach_ = 1.0;

  // A point of the target's path together with the number of configurations
  // that reached it, which is what the curve is coloured by.
  struct waypoint {
    std::array<double, 3> at;
    int found;
  };
  std::vector<waypoint> trail_;
  int last_count_ = 0;
  solver::status last_state_ = solver::status::ok;
  long ticks_ = 0;
  long refusals_ = 0;
  std::array<long, solver::max_configurations + 1> seen_{};
  double worst_residual_ = 0.0;

  varietas::chain<double> refused_;
  std::vector<std::string> refused_names_;
  std::string refused_root_;
  std::string refused_reason_;
  std::string refused_name_;
  std::string solved_name_;
  double label_height_ = 2.0;
  std::vector<double> refused_offset_{0.0, -2.9, 0.0};
  std::size_t refused_joints_ = 0;

  rclcpp::Publisher<sensor_msgs::msg::JointState>::SharedPtr joint_states_;
  rclcpp::Publisher<sensor_msgs::msg::JointState>::SharedPtr refused_states_;
  std::shared_ptr<tf2_ros::StaticTransformBroadcaster> static_tf_;
  rclcpp::Publisher<visualization_msgs::msg::MarkerArray>::SharedPtr markers_;
  rclcpp::TimerBase::SharedPtr timer_;
  rclcpp::TimerBase::SharedPtr report_timer_;
  rclcpp::Time start_;
};

int main(int argc, char** argv) {
  rclcpp::init(argc, argv);
  rclcpp::spin(std::make_shared<branch_node>());
  rclcpp::shutdown();
  return 0;
}
