"""Poses a URDF from the exact chain varietas recovers from it.

robot_state_publisher draws the arm from the file; the sweep node places a
marker at the tool pose computed from the exactly recovered chain. Run with

    ros2 launch varietas_demo sweep.launch.py

and pass urdf:=/path/to/model.urdf to use a model other than the iiwa fixture.
"""

import os

from ament_index_python.packages import get_package_share_directory
from launch import LaunchDescription
from launch.actions import DeclareLaunchArgument, OpaqueFunction
from launch.substitutions import LaunchConfiguration
from launch_ros.actions import Node


def _launch(context, *args, **kwargs):
    urdf_path = LaunchConfiguration("urdf").perform(context)
    with open(urdf_path, "r") as handle:
        robot_description = handle.read()

    rviz_config = os.path.join(
        get_package_share_directory("varietas_demo"), "rviz", "sweep.rviz"
    )

    return [
        Node(
            package="robot_state_publisher",
            executable="robot_state_publisher",
            output="screen",
            parameters=[{"robot_description": robot_description}],
        ),
        Node(
            package="varietas_demo",
            executable="sweep_node",
            output="screen",
            parameters=[
                {
                    "urdf": urdf_path,
                    "period": LaunchConfiguration("period"),
                }
            ],
        ),
        Node(
            package="rviz2",
            executable="rviz2",
            arguments=["-d", rviz_config],
            output="log",
        ),
    ]


def generate_launch_description():
    default_urdf = os.path.join(
        get_package_share_directory("varietas_urdf"), "data", "lbr_iiwa14.urdf"
    )
    return LaunchDescription(
        [
            DeclareLaunchArgument("urdf", default_value=default_urdf),
            DeclareLaunchArgument("period", default_value="12.0"),
            OpaqueFunction(function=_launch),
        ]
    )
