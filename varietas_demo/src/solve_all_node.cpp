// Every configuration, beside the one an iteration finds.
//
// The arm is general_3r, whose three axes are pairwise skew. Nothing about
// it can be decoupled and the symbolic solve over Q(x, y, z) does not finish;
// its solver was reconstructed from fixed poses during this build. Two copies
// of the arm follow the same moving target. The first is driven by that
// generated solver and draws every configuration it returns. The second is
// driven the way a robotics codebase without a closed form would drive it: a
// damped least squares iteration, warm-started from where it was on the
// previous tick, as a tracking controller is.
//
// The contrast is the demonstration. The generated solver returns all the
// configurations there are, the number of them being dim A, and returns none
// when the target leaves the reachable set. The iteration returns one, the
// branch it happens to be following; when that branch ceases to exist it is
// thrown onto another or stalls short of the target, and it has no way to say
// that three others were available all along. The iterating arm is drawn in
// the colour of whichever branch it is on, so a jump shows as a change of
// colour, and the branches keep their colours from tick to tick by continuity
// rather than by the order the eigensolver happened to return them in.
//
// Every tick is also written to a trace file, when one is named, so that a
// recording can carry the numbers as properly typeset text rather than as
// RViz text markers, which cannot draw a space.

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <memory>
#include <numeric>
#include <string>
#include <vector>

#include <Eigen/Dense>

#include <rclcpp/rclcpp.hpp>
#include <std_msgs/msg/color_rgba.hpp>
#include <visualization_msgs/msg/marker.hpp>
#include <visualization_msgs/msg/marker_array.hpp>

#include <urdf/model.h>

#include "varietas/kinematics/evaluate.hpp"
#include "varietas/urdf/urdf_chain.hpp"

#include "general_ik.hpp"

using namespace std::chrono_literals;

namespace {

using solver = varietas_generated::general_ik;
using provenance = varietas_generated::general_ik_provenance;
using vec3 = varietas::vector3<double>;
using config = std::array<double, 3>;

constexpr int kBranches = static_cast<int>(solver::dimension);

struct rgb {
  double r, g, b;
};

// One colour per branch, equally saturated, so that no branch reads as the
// answer and the others as alternatives.
const rgb kColours[4] = {
    {0.35, 0.70, 1.00},
    {1.00, 0.62, 0.22},
    {0.40, 0.90, 0.50},
    {0.92, 0.45, 0.86},
};

double wrapped(double a) {
  a = std::fmod(a + M_PI, 2.0 * M_PI);
  return a < 0 ? a + M_PI : a - M_PI;
}

double joint_distance(const config& a, const config& b) {
  double d = 0.0;
  for (std::size_t i = 0; i < 3; ++i) {
    d = std::max(d, std::abs(wrapped(a[i] - b[i])));
  }
  return d;
}

geometry_msgs::msg::Point point_of(const vec3& v) {
  geometry_msgs::msg::Point p;
  p.x = v[0];
  p.y = v[1];
  p.z = v[2];
  return p;
}

// Rotates a cylinder marker, which runs along its own z, onto the segment a-b.
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
  m.pose.orientation.w = n / 2.0;
}

double median(std::vector<double> v) {
  if (v.empty()) {
    return 0.0;
  }
  std::nth_element(v.begin(), v.begin() + static_cast<std::ptrdiff_t>(v.size() / 2), v.end());
  return v[v.size() / 2];
}

}  // namespace

class solve_all_node : public rclcpp::Node {
  static constexpr std::chrono::milliseconds kTick{33};
  static constexpr int kStride = 8;  // marker ids per drawn arm

 public:
  solve_all_node() : rclcpp::Node("varietas_solve_all") {
    const std::string urdf_path = declare_parameter<std::string>("urdf", "");
    period_ = declare_parameter<double>("period", 24.0);
    offset_ = declare_parameter<std::vector<double>>("iteration_offset", {0.0, -4.4, 0.0});
    labels_ = declare_parameter<bool>("labels", true);
    const std::string trace_path = declare_parameter<std::string>("trace", "");

    urdf::Model model;
    if (urdf_path.empty() || !model.initFile(urdf_path)) {
      RCLCPP_FATAL(get_logger(), "could not parse urdf '%s'", urdf_path.c_str());
      throw std::runtime_error("urdf");
    }
    varietas::chain<varietas::rational> exact;
    const auto imported = varietas::urdf_import::chain_from_model(
        model, model.getRoot()->name, varietas::urdf_import::sole_tip_link(model), exact);
    if (!imported.ok() || exact.degrees_of_freedom() != 3) {
      RCLCPP_FATAL(get_logger(), "the demonstration wants the three-joint general_3r arm");
      throw std::runtime_error("import");
    }
    robot_ = varietas::chain_cast<double>(exact);
    name_ = exact.name();
    measure_reach();

    if (!trace_path.empty()) {
      trace_ = std::fopen(trace_path.c_str(), "w");
      if (trace_ != nullptr) {
        std::fprintf(trace_, "wall,count,iteration_ok,iteration_error,jumps,stalls,"
                             "solve_us,iteration_us,iteration_branch\n");
      }
    }

    RCLCPP_INFO(get_logger(),
                "%s: reach %.2f m; solver reconstructed in %.1f s over %d primes and %d "
                "samples, %d exact checks; decoupling: %s",
                name_.c_str(), reach_, provenance::seconds, provenance::primes,
                provenance::samples, provenance::exact_checks, provenance::decoupling);

    markers_ = create_publisher<visualization_msgs::msg::MarkerArray>("varietas_markers", 10);
    timer_ = create_wall_timer(kTick, [this] { tick(); });
    report_timer_ = create_wall_timer(3s, [this] { report(); });
    start_ = now();
  }

  ~solve_all_node() override {
    if (trace_ != nullptr) {
      std::fclose(trace_);
    }
  }

 private:
  config forward(const config& q) const {
    const std::vector<double> v(q.begin(), q.end());
    const vec3 p = varietas::forward_kinematics(robot_, v).translation();
    return {p[0], p[1], p[2]};
  }

  void measure_reach() {
    double best = 0.0;
    const int steps = 25;
    for (int a = 0; a < steps; ++a) {
      for (int b = 0; b < steps; ++b) {
        for (int c = 0; c < steps; ++c) {
          const config q{-M_PI + 2.0 * M_PI * a / (steps - 1), -M_PI + 2.0 * M_PI * b / (steps - 1),
                         -M_PI + 2.0 * M_PI * c / (steps - 1)};
          const config p = forward(q);
          best = std::max(best, std::sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]));
        }
      }
    }
    reach_ = best;
  }

  // A closed path, planned against the arm's own map of how many
  // configurations exist where. The base turns about z at the origin, so the
  // count depends only on the distance from that axis and the height, and for
  // this arm all four configurations exist in a band about 1.4 to 2.0 m out at
  // a height of about 0.4 m, two on either side of it, and none beyond about
  // 2.4 m. The target circles inside the band and makes one excursion a
  // revolution out through the region of two and briefly past the reach, so
  // that branches die and are born in view and the iteration is led along a
  // branch to the place where it ends. Over a revolution the count is four
  // about 55% of the time, two about 43% and zero about 2%.
  config target_at(double t) const {
    const double u = 2.0 * M_PI * t / period_;
    const double away = std::remainder(u - M_PI, 2.0 * M_PI);
    const double radius = 1.62 + 0.20 * std::sin(2.0 * u) + 1.00 * std::exp(-(away / 0.4) * (away / 0.4));
    const double height = 0.38 + 0.22 * std::sin(3.0 * u);
    return {radius * std::cos(u), radius * std::sin(u), height};
  }

  // Damped least squares, warm-started. A finite-difference Jacobian, so that
  // it needs nothing but the forward map; a damping of 0.05, the order a
  // tracking controller uses, which keeps it from lunging when the target
  // leaves the workspace (at 1e-3 it swung metres from one tick to the next
  // there, which is real behaviour but not a fair picture of the method); and
  // a budget of iterations per tick, as a controller has.
  bool iterate(const config& target, config& q, double& error) const {
    constexpr double damping = 0.05;
    constexpr double step = 1e-7;
    for (int it = 0; it < 60; ++it) {
      const config p = forward(q);
      Eigen::Vector3d e(target[0] - p[0], target[1] - p[1], target[2] - p[2]);
      error = e.norm();
      if (error < 1e-9) {
        return true;
      }
      Eigen::Matrix3d J;
      for (int j = 0; j < 3; ++j) {
        config moved = q;
        moved[static_cast<std::size_t>(j)] += step;
        const config pm = forward(moved);
        for (int i = 0; i < 3; ++i) {
          J(i, j) = (pm[static_cast<std::size_t>(i)] - p[static_cast<std::size_t>(i)]) / step;
        }
      }
      const Eigen::Vector3d dq =
          J.transpose() * (J * J.transpose() + damping * damping * Eigen::Matrix3d::Identity())
                              .ldlt()
                              .solve(e);
      for (int j = 0; j < 3; ++j) {
        q[static_cast<std::size_t>(j)] = wrapped(q[static_cast<std::size_t>(j)] + dq[j]);
      }
    }
    const config p = forward(q);
    error = std::sqrt(std::pow(target[0] - p[0], 2) + std::pow(target[1] - p[1], 2) +
                      std::pow(target[2] - p[2], 2));
    return error < 1e-9;
  }

  // Gives each returned configuration the colour of the nearest one from the
  // previous tick, over all assignments, so that a branch keeps its colour for
  // as long as it exists. A configuration with no predecessor takes a colour
  // no surviving branch is using.
  std::vector<int> colour_by_continuity(const std::vector<config>& now_found) {
    const std::size_t n = now_found.size();
    std::vector<int> colour(n, -1);
    std::vector<int> perm(static_cast<std::size_t>(kBranches));
    std::iota(perm.begin(), perm.end(), 0);
    double best = 1e300;
    std::vector<int> best_perm = perm;
    // perm[i] is the previous slot matched to new configuration i.
    do {
      double cost = 0.0;
      for (std::size_t i = 0; i < n; ++i) {
        const int slot = perm[i];
        cost += previous_alive_[static_cast<std::size_t>(slot)]
                    ? joint_distance(now_found[i], previous_[static_cast<std::size_t>(slot)])
                    : 3.5;  // a birth costs about as much as a large jump
      }
      if (cost < best) {
        best = cost;
        best_perm = perm;
      }
    } while (std::next_permutation(perm.begin(), perm.end()));

    std::array<bool, 4> alive{};
    for (std::size_t i = 0; i < n; ++i) {
      colour[i] = best_perm[i];
      previous_[static_cast<std::size_t>(best_perm[i])] = now_found[i];
      alive[static_cast<std::size_t>(best_perm[i])] = true;
    }
    previous_alive_ = alive;
    return colour;
  }

  void tick() {
    const double t = (now() - start_).seconds();
    const config target = target_at(t);

    // The generated solver, timed alone.
    std::array<double, solver::dimension * solver::num_unknowns> raw{};
    solver::status state{};
    const auto s0 = std::chrono::steady_clock::now();
    const int found = solver::solve(target.data(), raw.data(), kBranches, &state);
    const double solve_us =
        std::chrono::duration<double, std::micro>(std::chrono::steady_clock::now() - s0).count();

    std::vector<config> configurations;
    for (int k = 0; k < std::max(found, 0); ++k) {
      configurations.push_back({2.0 * std::atan(raw[static_cast<std::size_t>(3 * k)]),
                                2.0 * std::atan(raw[static_cast<std::size_t>(3 * k + 1)]),
                                2.0 * std::atan(raw[static_cast<std::size_t>(3 * k + 2)])});
    }
    const std::vector<int> colours = colour_by_continuity(configurations);

    // The iteration, from where it was.
    const config before = iterate_q_;
    double error = 0.0;
    const auto i0 = std::chrono::steady_clock::now();
    const bool ok = iterate(target, iterate_q_, error);
    const double iteration_us =
        std::chrono::duration<double, std::micro>(std::chrono::steady_clock::now() - i0).count();

    // Which branch is it on? The nearest returned configuration, when it has
    // converged onto one.
    int on_branch = -1;
    if (ok) {
      double nearest = 1e300;
      for (std::size_t k = 0; k < configurations.size(); ++k) {
        const double d = joint_distance(iterate_q_, configurations[k]);
        if (d < nearest) {
          nearest = d;
          on_branch = colours[k];
        }
      }
      if (nearest > 1e-4) {
        on_branch = -1;
      }
    }
    if (ok && last_branch_ >= 0 && on_branch >= 0 && on_branch != last_branch_) {
      ++jumps_;
    }
    if (!ok && last_ok_) {
      ++stalls_;
    }
    if (on_branch >= 0) {
      last_branch_ = on_branch;
    }
    last_ok_ = ok;
    (void)before;

    solve_times_.push_back(solve_us);
    iteration_times_.push_back(iteration_us);
    if (solve_times_.size() > 300) {
      solve_times_.erase(solve_times_.begin());
      iteration_times_.erase(iteration_times_.begin());
    }
    ++ticks_;
    counts_[static_cast<std::size_t>(std::max(found, 0))]++;

    trail_.push_back({target, found});
    const auto capacity = static_cast<std::size_t>(period_ * 1000.0 / kTick.count());
    while (trail_.size() > capacity) {
      trail_.erase(trail_.begin());
    }

    if (trace_ != nullptr) {
      const double wall =
          std::chrono::duration<double>(std::chrono::system_clock::now().time_since_epoch())
              .count();
      std::fprintf(trace_, "%.4f,%d,%d,%.3e,%ld,%ld,%.2f,%.2f,%d\n", wall, found, ok ? 1 : 0,
                   error, jumps_, stalls_, median(solve_times_), median(iteration_times_),
                   on_branch);
      std::fflush(trace_);
    }

    publish(target, configurations, colours, ok, on_branch);
  }

  // One arm as cylinders and spheres, offset by `shift`, in marker ids
  // [base, base + kStride).
  void draw_arm(visualization_msgs::msg::MarkerArray& array, const config& q, const vec3& shift,
                const rgb& colour, double alpha, int base, const rclcpp::Time& stamp,
                bool axes) const {
    const std::vector<double> v(q.begin(), q.end());
    const auto frames = varietas::link_frames(robot_, v);
    std::vector<vec3> joints;
    for (const auto& f : frames) {
      const vec3 p = f.translation() + shift;
      if (joints.empty() ||
          std::abs(p[0] - joints.back()[0]) + std::abs(p[1] - joints.back()[1]) +
                  std::abs(p[2] - joints.back()[2]) > 1e-9) {
        joints.push_back(p);
      }
    }
    int used = 0;
    for (std::size_t seg = 0; seg + 1 < joints.size() && used < kStride - 2; ++seg, ++used) {
      auto limb = marker("limb", base + used, stamp);
      limb.type = visualization_msgs::msg::Marker::CYLINDER;
      limb.scale.x = limb.scale.y = 0.06;
      limb.color.r = colour.r;
      limb.color.g = colour.g;
      limb.color.b = colour.b;
      limb.color.a = alpha;
      span(limb, joints[seg], joints[seg + 1]);
      array.markers.push_back(limb);
    }
    for (; used < kStride - 2; ++used) {
      array.markers.push_back(erased("limb", base + used, stamp));
    }
    auto knuckles = marker("knuckle", base, stamp);
    knuckles.type = visualization_msgs::msg::Marker::SPHERE_LIST;
    knuckles.scale.x = knuckles.scale.y = knuckles.scale.z = 0.095;
    knuckles.color.r = colour.r;
    knuckles.color.g = colour.g;
    knuckles.color.b = colour.b;
    knuckles.color.a = 1.0;
    for (const auto& j : joints) {
      knuckles.points.push_back(point_of(j));
    }
    array.markers.push_back(knuckles);

    // The joint axes, which are the point of this arm: none is parallel or
    // perpendicular to another.
    auto lines = marker("axis", base, stamp);
    if (!axes) {
      array.markers.push_back(erased("axis", base, stamp));
      return;
    }
    lines.type = visualization_msgs::msg::Marker::LINE_LIST;
    lines.scale.x = 0.014;
    lines.color.r = lines.color.g = lines.color.b = 0.85;
    lines.color.a = 0.8;
    // link_frames leads with the base frame, so joint i's frame is entry i + 1,
    // and a revolute joint's rotation fixes its own axis, so the axis in the
    // world is that frame's rotation applied to it.
    for (std::size_t i = 0; i < robot_.joints().size(); ++i) {
      const auto& j = robot_.joints()[i];
      if (!j.is_actuated()) {
        continue;
      }
      const auto& f = frames[i + 1];
      const vec3 axis = f.rotation() * j.axis;
      const vec3 centre = f.translation() + shift;
      lines.points.push_back(point_of(centre - 0.28 * axis));
      lines.points.push_back(point_of(centre + 0.28 * axis));
    }
    array.markers.push_back(lines);
  }

  void publish(const config& target, const std::vector<config>& configurations,
               const std::vector<int>& colours, bool ok, int on_branch) {
    visualization_msgs::msg::MarkerArray array;
    const auto stamp = now();
    const vec3 shift(offset_.size() > 0 ? offset_[0] : 0.0, offset_.size() > 1 ? offset_[1] : 0.0,
                     offset_.size() > 2 ? offset_[2] : 0.0);
    const vec3 origin(0.0, 0.0, 0.0);

    // All configurations, in their branch colours.
    std::array<bool, 4> drawn{};
    for (std::size_t k = 0; k < configurations.size(); ++k) {
      const int c = colours[k];
      drawn[static_cast<std::size_t>(c)] = true;
      draw_arm(array, configurations[k], origin, kColours[c], 0.95, c * kStride, stamp, false);
    }
    for (int c = 0; c < kBranches; ++c) {
      if (!drawn[static_cast<std::size_t>(c)]) {
        for (int i = 0; i < kStride - 2; ++i) {
          array.markers.push_back(erased("limb", c * kStride + i, stamp));
        }
        array.markers.push_back(erased("knuckle", c * kStride, stamp));
        array.markers.push_back(erased("axis", c * kStride, stamp));
      }
    }

    // The iterating arm, in the colour of its branch, or grey when it has
    // stalled short of the target.
    const rgb stalled{0.62, 0.62, 0.66};
    const rgb& shade = (ok && on_branch >= 0) ? kColours[on_branch] : stalled;
    draw_arm(array, iterate_q_, shift, shade, ok ? 0.95 : 0.8, 4 * kStride, stamp, true);

    // Plinths, targets and trails for both.
    for (int side = 0; side < 2; ++side) {
      const vec3 at = side == 0 ? origin : shift;
      auto plinth = marker("plinth", side, stamp);
      plinth.type = visualization_msgs::msg::Marker::CYLINDER;
      plinth.scale.x = plinth.scale.y = 0.42;
      plinth.scale.z = 0.08;
      plinth.pose.position.x = at[0];
      plinth.pose.position.y = at[1];
      plinth.pose.position.z = -0.04;
      plinth.color.r = plinth.color.g = plinth.color.b = 0.36;
      plinth.color.a = 1.0;
      array.markers.push_back(plinth);

      auto goal = marker("target", side, stamp);
      goal.type = visualization_msgs::msg::Marker::SPHERE;
      goal.scale.x = goal.scale.y = goal.scale.z = 0.14;
      goal.pose.position.x = target[0] + at[0];
      goal.pose.position.y = target[1] + at[1];
      goal.pose.position.z = target[2] + at[2];
      const bool reached = side == 0 ? !configurations.empty() : ok;
      goal.color.r = reached ? 0.97 : 0.96;
      goal.color.g = reached ? 0.97 : 0.76;
      goal.color.b = reached ? 0.97 : 0.18;
      goal.color.a = 1.0;
      array.markers.push_back(goal);

      auto path = marker("trail", side, stamp);
      path.type = visualization_msgs::msg::Marker::LINE_LIST;
      path.scale.x = 0.013;
      for (std::size_t i = 1; i < trail_.size(); ++i) {
        std_msgs::msg::ColorRGBA c;
        c.a = 0.8f;
        const int n = trail_[i].found;
        if (n <= 0) {
          c.r = 0.96f; c.g = 0.76f; c.b = 0.18f;
        } else if (n < kBranches) {
          c.r = 0.55f; c.g = 0.80f; c.b = 0.95f;
        } else {
          c.r = 0.93f; c.g = 0.93f; c.b = 0.95f;
        }
        geometry_msgs::msg::Point a = point_of(vec3(trail_[i - 1].at[0], trail_[i - 1].at[1],
                                                    trail_[i - 1].at[2]) + at);
        geometry_msgs::msg::Point b =
            point_of(vec3(trail_[i].at[0], trail_[i].at[1], trail_[i].at[2]) + at);
        path.points.push_back(a);
        path.colors.push_back(c);
        path.points.push_back(b);
        path.colors.push_back(c);
      }
      array.markers.push_back(path);
    }

    if (labels_) {
      array.markers.push_back(label(0, origin, "varietas:_every_configuration\n" +
                                                   std::to_string(configurations.size()) + "/" +
                                                   std::to_string(kBranches),
                                    stamp));
      array.markers.push_back(label(1, shift,
                                    std::string("damped_least_squares:_one\n") +
                                        (ok ? "converged" : "stalled"),
                                    stamp));
    }
    markers_->publish(array);
  }

  visualization_msgs::msg::Marker label(int id, const vec3& at, const std::string& text,
                                        const rclcpp::Time& stamp) const {
    auto m = marker("label", id, stamp);
    m.type = visualization_msgs::msg::Marker::TEXT_VIEW_FACING;
    m.pose.position.x = at[0];
    m.pose.position.y = at[1];
    m.pose.position.z = 2.5;
    m.scale.z = 0.18;
    m.color.r = m.color.g = m.color.b = 0.93;
    m.color.a = 0.95;
    m.text = text;
    return m;
  }

  visualization_msgs::msg::Marker marker(const std::string& ns, int id,
                                         const rclcpp::Time& stamp) const {
    visualization_msgs::msg::Marker m;
    m.header.frame_id = "world";
    m.header.stamp = stamp;
    m.ns = ns;
    m.id = id;
    m.action = visualization_msgs::msg::Marker::ADD;
    m.pose.orientation.w = 1.0;
    return m;
  }

  visualization_msgs::msg::Marker erased(const std::string& ns, int id,
                                         const rclcpp::Time& stamp) const {
    auto m = marker(ns, id, stamp);
    m.action = visualization_msgs::msg::Marker::DELETE;
    return m;
  }

  void report() {
    RCLCPP_INFO(get_logger(),
                "generated solver %.1f us median, %ld ticks [0:%ld 2:%ld 4:%ld]; iteration %.1f "
                "us median, %ld branch jumps, %ld stalls",
                median(solve_times_), ticks_, counts_[0], counts_[2], counts_[4],
                median(iteration_times_), jumps_, stalls_);
  }

  varietas::chain<double> robot_;
  std::string name_;
  double period_ = 24.0;
  double reach_ = 1.0;
  std::vector<double> offset_;
  bool labels_ = true;
  std::FILE* trace_ = nullptr;

  std::array<config, 4> previous_{};
  std::array<bool, 4> previous_alive_{};
  config iterate_q_{0.3, 0.4, -0.5};
  int last_branch_ = -1;
  bool last_ok_ = true;
  long jumps_ = 0;
  long stalls_ = 0;
  long ticks_ = 0;
  std::array<long, 5> counts_{};
  std::vector<double> solve_times_;
  std::vector<double> iteration_times_;

  struct waypoint {
    config at;
    int found;
  };
  std::vector<waypoint> trail_;

  rclcpp::Publisher<visualization_msgs::msg::MarkerArray>::SharedPtr markers_;
  rclcpp::TimerBase::SharedPtr timer_;
  rclcpp::TimerBase::SharedPtr report_timer_;
  rclcpp::Time start_;
};

int main(int argc, char** argv) {
  rclcpp::init(argc, argv);
  rclcpp::spin(std::make_shared<solve_all_node>());
  rclcpp::shutdown();
  return 0;
}
