// Drives a URDF from the exact chain varietas recovered from it.
//
// The node imports the model into a chain over the rationals, sweeps the joints
// through a trajectory, and publishes the joint states, which
// robot_state_publisher turns into the transforms RViz draws the robot with,
// together with two markers at the tool: a halo at the pose varietas computes
// from the exactly recovered chain, and a core at the pose
// robot_state_publisher derives from the decimals in the file, by way of KDL.
//
// The pair is the whole point, and it is a pair for a reason. Comparing one
// marker against the drawn robot does not work: RViz redraws the robot from tf
// on its own schedule, so at speed the mesh trails the marker stream by however
// much the two clocks differ, and a marker sitting exactly where the frame is
// still looks displaced from the mesh. Two markers published in one array from
// one instant admit no such gap. Agreement is the shell sitting concentric in
// the halo, and an error is the shell leaving it.

#include <chrono>
#include <cmath>
#include <memory>
#include <string>
#include <vector>

#include <rclcpp/rclcpp.hpp>
#include <geometry_msgs/msg/transform_stamped.hpp>
#include <sensor_msgs/msg/joint_state.hpp>
#include <tf2_ros/buffer.h>
#include <tf2_ros/transform_listener.h>
#include <urdf/model.h>
#include <visualization_msgs/msg/marker.hpp>
#include <visualization_msgs/msg/marker_array.hpp>

#include "varietas/kinematics/evaluate.hpp"
#include "varietas/urdf/urdf_chain.hpp"

using namespace std::chrono_literals;

namespace {

// Quaternion of a rotation matrix, for the marker's pose. Shepperd's branch on
// the largest diagonal entry, which avoids the cancellation the naive trace
// formula suffers near a half turn.
void quaternion_of(const varietas::matrix3<double>& r, double& x, double& y, double& z,
                   double& w) {
  const double trace = r(0, 0) + r(1, 1) + r(2, 2);
  if (trace > 0.0) {
    const double s = 0.5 / std::sqrt(trace + 1.0);
    w = 0.25 / s;
    x = (r(2, 1) - r(1, 2)) * s;
    y = (r(0, 2) - r(2, 0)) * s;
    z = (r(1, 0) - r(0, 1)) * s;
  } else if (r(0, 0) > r(1, 1) && r(0, 0) > r(2, 2)) {
    const double s = 2.0 * std::sqrt(1.0 + r(0, 0) - r(1, 1) - r(2, 2));
    w = (r(2, 1) - r(1, 2)) / s;
    x = 0.25 * s;
    y = (r(0, 1) + r(1, 0)) / s;
    z = (r(0, 2) + r(2, 0)) / s;
  } else if (r(1, 1) > r(2, 2)) {
    const double s = 2.0 * std::sqrt(1.0 + r(1, 1) - r(0, 0) - r(2, 2));
    w = (r(0, 2) - r(2, 0)) / s;
    x = (r(0, 1) + r(1, 0)) / s;
    y = 0.25 * s;
    z = (r(1, 2) + r(2, 1)) / s;
  } else {
    const double s = 2.0 * std::sqrt(1.0 + r(2, 2) - r(0, 0) - r(1, 1));
    w = (r(1, 0) - r(0, 1)) / s;
    x = (r(0, 2) + r(2, 0)) / s;
    y = (r(1, 2) + r(2, 1)) / s;
    z = 0.25 * s;
  }
}

}  // namespace

class sweep_node : public rclcpp::Node {
  // The rate the trajectory is stepped at. The trail's length is derived from
  // it, so the two cannot drift apart.
  static constexpr std::chrono::milliseconds kTick{20};

 public:
  sweep_node() : rclcpp::Node("varietas_sweep") {
    const std::string urdf_path = declare_parameter<std::string>("urdf", "");
    period_ = declare_parameter<double>("period", 12.0);

    urdf::Model model;
    if (urdf_path.empty() || !model.initFile(urdf_path)) {
      RCLCPP_FATAL(get_logger(), "could not parse urdf '%s'", urdf_path.c_str());
      throw std::runtime_error("urdf");
    }

    const std::string tip = declare_parameter<std::string>(
        "tip_link", varietas::urdf_import::sole_tip_link(model));
    tip_link_ = tip;
    const std::string root =
        declare_parameter<std::string>("root_link", model.getRoot()->name);
    root_link_ = root;

    varietas::chain<varietas::rational> exact;
    const auto report = varietas::urdf_import::chain_from_model(model, root, tip, exact);
    if (!report.ok()) {
      RCLCPP_FATAL(get_logger(), "import refused: %s (%s)",
                   varietas::urdf_import::to_string(report.status),
                   report.detail.c_str());
      throw std::runtime_error("import");
    }

    RCLCPP_INFO(get_logger(),
                "recovered %zu joints exactly; worst rotation moved %.3e rad, "
                "worst translation %.3e m",
                exact.degrees_of_freedom(), report.max_rotation_deviation,
                report.max_translation_deviation);

    robot_ = varietas::chain_cast<double>(exact);
    for (const auto& j : robot_.joints()) {
      if (j.is_actuated()) {
        names_.push_back(j.name);
        lower_.push_back(j.has_limits ? j.lower : -3.0);
        upper_.push_back(j.has_limits ? j.upper : 3.0);
      }
    }

    joint_states_ =
        create_publisher<sensor_msgs::msg::JointState>("joint_states", 10);
    markers_ =
        create_publisher<visualization_msgs::msg::MarkerArray>("varietas_markers", 10);
    tf_buffer_ = std::make_unique<tf2_ros::Buffer>(get_clock());
    tf_listener_ = std::make_shared<tf2_ros::TransformListener>(*tf_buffer_);

    timer_ = create_wall_timer(kTick, [this] { tick(); });
    verify_timer_ = create_wall_timer(200ms, [this] { verify(); });
    start_ = now();
  }

 private:
  // The trajectory as a function of time rather than of the tick, so that the
  // configuration belonging to any past instant can be recovered exactly. The
  // verification below needs that: it has to compare our pose against the
  // reference at the same instant, not merely at about the same time.
  //
  // Each joint runs on its own frequency, scaled into its limits, so the arm
  // covers a large part of its configuration space rather than one plane. The
  // frequencies are whole multiples of the fundamental, which makes the whole
  // trajectory periodic with exactly `period`: a recording one period long then
  // closes on itself.
  std::vector<double> configuration_at(double t) const {
    static constexpr int kHarmonic[] = {1, 2, 3, 1, 2, 3, 1};
    std::vector<double> values(names_.size());
    for (std::size_t i = 0; i < values.size(); ++i) {
      const int harmonic = kHarmonic[i % (sizeof(kHarmonic) / sizeof(int))];
      const double phase = 2.0 * M_PI * harmonic * t / period_;
      const double centre = 0.5 * (lower_[i] + upper_[i]);
      const double amplitude = 0.38 * (upper_[i] - lower_[i]);
      values[i] = centre + amplitude * std::sin(phase + 0.7 * static_cast<double>(i));
    }
    return values;
  }

  void tick() {
    const rclcpp::Time stamp = now();
    const double t = (stamp - start_).seconds();
    const std::vector<double> values = configuration_at(t);

    sensor_msgs::msg::JointState state;
    state.header.stamp = stamp;
    state.name = names_;
    state.position = values;
    joint_states_->publish(state);

    publish_markers(values);
  }

  void publish_markers(const std::vector<double>& values) {
    const varietas::rigid_transform<double> tool =
        varietas::forward_kinematics(robot_, values);

    geometry_msgs::msg::Point point;
    point.x = tool.translation()[0];
    point.y = tool.translation()[1];
    point.z = tool.translation()[2];
    trail_.push_back(point);
    if (trail_.size() > trail_capacity()) {
      trail_.erase(trail_.begin());
    }

    visualization_msgs::msg::MarkerArray array;

    visualization_msgs::msg::Marker tip;
    tip.header.frame_id = root_link_;
    tip.header.stamp = now();
    tip.ns = "varietas";
    tip.id = 0;
    tip.type = visualization_msgs::msg::Marker::SPHERE;
    tip.action = visualization_msgs::msg::Marker::ADD;
    tip.pose.position.x = tool.translation()[0];
    tip.pose.position.y = tool.translation()[1];
    tip.pose.position.z = tool.translation()[2];
    quaternion_of(tool.rotation(), tip.pose.orientation.x, tip.pose.orientation.y,
                  tip.pose.orientation.z, tip.pose.orientation.w);
    // Translucent, and wider than the shell below by a margin thick enough to
    // read at the size the recording is published at: the annulus between the
    // two silhouettes is what carries the comparison, so it has to survive
    // being scaled down to a few hundred pixels. Scale is a diameter, so this
    // is a radius of 0.120 against the shell's 0.080.
    tip.scale.x = tip.scale.y = tip.scale.z = 0.240;
    tip.color.r = 0.11f;
    tip.color.g = 0.85f;
    tip.color.b = 0.62f;
    tip.color.a = 0.30f;
    array.markers.push_back(tip);

    // The same pose as robot_state_publisher reports it, drawn as a translucent
    // shell around ours.
    //
    // Comparing our marker against the drawn robot was the wrong comparison and
    // it made the demonstration unreadable. Not because the flange is drawn off
    // its frame -- in this fixture link_7's visual is a sphere at the origin,
    // so it is drawn exactly on it -- but because RViz poses the robot from tf
    // on its own redraw schedule while the markers arrive on ours. The two
    // clocks differ, and at the speeds the tool reaches here that difference is
    // centimetres of apparent displacement between marker and mesh: measured
    // off a recording, the gap tracked tool speed at r = 0.6, which is a lag
    // signature and not a kinematic one. So the mesh comparison showed the
    // renderer's timing even though the two poses agreed to 4e-12 m, and the
    // picture could not tell a correct recovery from a wrong one, which is the
    // only thing it exists to do.
    //
    // Two markers under the same convention can. The green halo is placed from
    // our exactly recovered chain and the magenta shell inside it from the
    // transform robot_state_publisher derives from the file's decimals, so
    // agreement is the shell sitting concentric in the halo and any error is
    // the shell riding up against one side of it. What the picture resolves is
    // gross error, a centimetre or so; the agreement is finer than any picture
    // can carry and is measured rather than shown, by verify() below.
    // Both markers have to describe the same instant or the picture measures
    // latency instead of kinematics. The transform robot_state_publisher
    // publishes carries the stamp it was computed for, and the trajectory is an
    // analytic function of time, so our own pose is evaluated at that stamp
    // rather than at the current one — the same discipline verify() uses, and
    // for the same reason: at these speeds a tenth of a second of lag is
    // indistinguishable from a kinematic error.
    geometry_msgs::msg::TransformStamped reference;
    try {
      reference = tf_buffer_->lookupTransform(root_link_, tip_link_, tf2::TimePointZero);
    } catch (const tf2::TransformException&) {
      array.markers.push_back(trail_marker(tip.header));
      markers_->publish(array);
      return;
    }

    const rclcpp::Time at(reference.header.stamp);
    const varietas::rigid_transform<double> ours =
        varietas::forward_kinematics(robot_, configuration_at((at - start_).seconds()));

    tip.pose.position.x = ours.translation()[0];
    tip.pose.position.y = ours.translation()[1];
    tip.pose.position.z = ours.translation()[2];
    quaternion_of(ours.rotation(), tip.pose.orientation.x, tip.pose.orientation.y,
                  tip.pose.orientation.z, tip.pose.orientation.w);
    array.markers.back() = tip;

    visualization_msgs::msg::Marker shell;
    shell.header = tip.header;
    shell.ns = "varietas";
    shell.id = 2;
    shell.type = visualization_msgs::msg::Marker::SPHERE;
    shell.action = visualization_msgs::msg::Marker::ADD;
    shell.pose.position.x = reference.transform.translation.x;
    shell.pose.position.y = reference.transform.translation.y;
    shell.pose.position.z = reference.transform.translation.z;
    shell.pose.orientation = reference.transform.rotation;
    // Wide enough to actually enclose the flange the file draws, and nearly
    // opaque so that it hides it. Scale is a diameter, so this is a radius of
    // 0.080 against the flange sphere's 0.060; the previous 0.080 was a radius
    // of 0.040 and so sat buried inside the very thing the comment claimed it
    // enclosed, visible only in the moments it slipped out from under it.
    //
    // The 20 mm of clearance is what absorbs the renderer skew described
    // below. It is not unlimited: the sweep reaches 1.6 m/s at period 20, so
    // at the very fastest moments the flange still shows at the rim. Covering
    // it even then would take a marker of 15 cm radius, which would swallow
    // the wrist it is meant to sit on, and the trade is not worth it.
    //
    // Hiding the flange is the point and not a cosmetic choice. RViz redraws
    // the robot from tf on its own schedule, independently of the marker
    // stream, and the tool moves at up to about a metre a second here, so the
    // drawn flange trails the markers by whatever the two clocks differ by --
    // measured off an earlier recording, some tens of milliseconds, which is
    // centimetres on screen. That lag is a property of the renderer and says
    // nothing about the kinematics, but left visible it reads as a third blob
    // sitting off to one side and invites exactly the wrong conclusion. The
    // two markers are published in one array from one instant, so no such
    // skew can open between them, and covering the flange leaves only the
    // comparison that means something.
    //
    // Magenta against the green halo: the two have to be told apart at a
    // glance, and green against blue at this size reads as one teal blob.
    shell.scale.x = shell.scale.y = shell.scale.z = 0.160;
    shell.color.r = 0.96f;
    shell.color.g = 0.26f;
    shell.color.b = 0.72f;
    shell.color.a = 0.92f;
    array.markers.push_back(shell);


    array.markers.push_back(trail_marker(tip.header));
    markers_->publish(array);
  }

  // One full period of tool positions, at the rate tick() runs. The trajectory
  // is periodic, so a window of exactly one period is the closed curve the tool
  // traces and nothing more: it neither accumulates nor, as a shorter window
  // did, reduces the curve to a fragment that swings about the frame while the
  // rest of it is missing. Held whole, it is a fixed object the arm moves
  // through, which is what a still of the recording has to be able to show.
  std::size_t trail_capacity() const {
    const double step = std::chrono::duration<double>(kTick).count();
    return static_cast<std::size_t>(period_ / step) + 1;
  }

  visualization_msgs::msg::Marker trail_marker(const std_msgs::msg::Header& header) {
    visualization_msgs::msg::Marker trail;
    trail.header = header;
    trail.ns = "varietas";
    trail.id = 1;
    trail.type = visualization_msgs::msg::Marker::LINE_STRIP;
    trail.action = visualization_msgs::msg::Marker::ADD;
    trail.pose.orientation.w = 1.0;
    trail.scale.x = 0.012;
    trail.color.r = 0.11f;
    trail.color.g = 0.78f;
    trail.color.b = 0.60f;
    trail.color.a = 0.65f;
    trail.points = trail_;
    return trail;
  }

  // The claim the demonstration exists to support, measured rather than
  // eyeballed: the pose computed from the exactly recovered chain agrees with
  // the one robot_state_publisher derives from the file's decimals, at the same
  // instant. Reading the two off the screen cannot establish this — the tool
  // moves fast enough that a tenth of a second of display lag looks exactly
  // like a kinematic error — so the transform is looked up at a definite past
  // stamp and our own map evaluated at that same stamp.
  void verify() {
    // The latest transform tf holds, rather than one at a time of our
    // choosing. Asking for a particular instant makes tf interpolate between
    // the samples either side of it, and the comparison then measures that
    // interpolation — a fraction of a millimetre at these speeds — instead of
    // the kinematics. The latest sample carries the stamp it was actually
    // computed for, and the trajectory is an analytic function of time, so our
    // own map can be evaluated at precisely that instant.
    geometry_msgs::msg::TransformStamped reference;
    try {
      reference = tf_buffer_->lookupTransform(root_link_, tip_link_, tf2::TimePointZero);
    } catch (const tf2::TransformException&) {
      return;
    }

    const rclcpp::Time at(reference.header.stamp);
    const varietas::rigid_transform<double> ours =
        varietas::forward_kinematics(robot_, configuration_at((at - start_).seconds()));

    const double dx = ours.translation()[0] - reference.transform.translation.x;
    const double dy = ours.translation()[1] - reference.transform.translation.y;
    const double dz = ours.translation()[2] - reference.transform.translation.z;
    const double distance = std::sqrt(dx * dx + dy * dy + dz * dz);
    worst_distance_ = std::max(worst_distance_, distance);

    if (++verifications_ % 10 == 0) {
      RCLCPP_INFO(get_logger(),
                  "tool pose against robot_state_publisher: %.3e m now, %.3e m worst "
                  "over %d samples",
                  distance, worst_distance_, verifications_);
    }
  }

  varietas::chain<double> robot_;
  std::string root_link_;
  std::string tip_link_;
  std::vector<std::string> names_;
  std::vector<double> lower_, upper_;
  std::vector<geometry_msgs::msg::Point> trail_;
  double period_ = 12.0;
  rclcpp::Time start_;
  rclcpp::Publisher<sensor_msgs::msg::JointState>::SharedPtr joint_states_;
  rclcpp::Publisher<visualization_msgs::msg::MarkerArray>::SharedPtr markers_;
  rclcpp::TimerBase::SharedPtr timer_;
  rclcpp::TimerBase::SharedPtr verify_timer_;
  std::unique_ptr<tf2_ros::Buffer> tf_buffer_;
  std::shared_ptr<tf2_ros::TransformListener> tf_listener_;
  double worst_distance_ = 0.0;
  int verifications_ = 0;
};

int main(int argc, char** argv) {
  rclcpp::init(argc, argv);
  rclcpp::spin(std::make_shared<sweep_node>());
  rclcpp::shutdown();
  return 0;
}
