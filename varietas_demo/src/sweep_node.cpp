// Drives a URDF from the exact chain varietas recovered from it.
//
// The node imports the model into a chain over the rationals, sweeps the joints
// through a trajectory, and publishes two things: the joint states, which
// robot_state_publisher turns into the transforms RViz draws the robot with,
// and a marker at the tool pose computed by varietas' own forward kinematics
// from the exactly recovered chain.
//
// The marker is the whole point. robot_state_publisher poses the robot from the
// decimals in the file, by way of KDL; the marker is placed from the exact
// chain, whose right angles are exact rather than truncated. If the recovery
// had changed the robot, the marker would drift off the tip as the arm moves.
// It does not, at any configuration, which is what the recording shows.

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

    timer_ = create_wall_timer(20ms, [this] { tick(); });
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
    // Large enough to enclose the link the robot_state_publisher draws, and
    // translucent, so that the two poses are seen to coincide rather than one
    // hiding the other.
    tip.scale.x = tip.scale.y = tip.scale.z = 0.16;
    tip.color.r = 0.11f;
    tip.color.g = 0.78f;
    tip.color.b = 0.60f;
    tip.color.a = 0.45f;
    array.markers.push_back(tip);


    // The path the tool has traced, kept to a fixed length so the recording
    // shows a ribbon rather than an accumulating tangle.
    geometry_msgs::msg::Point point;
    point.x = tool.translation()[0];
    point.y = tool.translation()[1];
    point.z = tool.translation()[2];
    trail_.push_back(point);
    if (trail_.size() > 240) {
      trail_.erase(trail_.begin());
    }

    visualization_msgs::msg::Marker trail;
    trail.header = tip.header;
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
    array.markers.push_back(trail);

    markers_->publish(array);
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
