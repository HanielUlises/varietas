# Building

ROS 2 Humble, C++17. Eigen and GMP are the only non-ROS dependencies.

## Dependencies

```sh
sudo apt install libeigen3-dev libgmp-dev
sudo apt install ros-humble-urdf ros-humble-kdl-parser ros-humble-robot-state-publisher ros-humble-rviz2
```

## Build

```sh
source /opt/ros/humble/setup.bash
colcon build
source install/setup.bash
```

A single package and what it depends on:

```sh
colcon build --packages-up-to varietas_kinematics
```

## Test

```sh
colcon test
colcon test-result --all --verbose
```

```sh
colcon test --packages-select varietas_codegen
```

## Packages

| Package | Contents | Depends on |
|---|---|---|
| `varietas_core` | Orders, polynomials, ideals, quotient algebra, solving | Eigen |
| `varietas_codegen` | Exact rationals and rational functions | `varietas_core`, GMP |
| `varietas_kinematics` | Chains, rationalisation, workspace, singularities | `varietas_core` |
| `varietas_ik` | Inverse kinematics over Q(pose), ready to emit | `varietas_kinematics`, `varietas_codegen` |
| `varietas_urdf` | URDF to exact chain, the audit, and `urdf_codegen` | `varietas_ik`, `urdf` |
| `varietas_demo` | RViz demonstration | `varietas_urdf`, `rclcpp` |

`varietas_core`, `varietas_codegen`, `varietas_kinematics` and `varietas_ik` are header-only interface targets.

## Without ROS

The algebra needs no ROS. `varietas_core` needs Eigen; `varietas_codegen` and `varietas_kinematics` additionally need GMP.

```sh
g++ -std=c++17 -O2 example.cpp \
  -Ivarietas_core/include -Ivarietas_codegen/include -Ivarietas_kinematics/include \
  -Ivarietas_ik/include \
  -I/usr/include/eigen3 -lgmpxx -lgmp
```

A header `urdf_codegen` produced needs neither: it includes `<cstddef>` and
`<cstdint>`, plus Eigen when it was emitted with the solver.

```sh
g++ -std=c++17 -O2 consumer.cpp -I. -I/usr/include/eigen3
```

`varietas_urdf` and `varietas_demo` require ROS.

## Running

```sh
ros2 run varietas_urdf urdf_report <file.urdf> [tip_link] [root_link]
ros2 run varietas_urdf urdf_codegen <file.urdf> <output.hpp> [--tip L] [--root L] [--coords xy|xyz]
ros2 launch varietas_demo sweep.launch.py urdf:=<file.urdf> period:=12.0
```
