"""Draws every configuration that reaches a moving target.

The solver was generated from the URDF during the build, with the base joint
swept out rather than adjoined. This launch moves a target through the
workspace and draws all the configurations the solver returns for it, so the
branches can be seen merging where the arm loses rank and vanishing where the
target leaves the reachable set. Run with

    ros2 launch varietas_demo branches.launch.py

Two arms stand in the frame. The one that is solved is the anthropomorphic 3R
model in varietas_demo/urdf, which the header emitted during the build is for;
a solver for another arm has to be emitted for it, so `urdf:=` alone will not
do. Beside it stands the KUKA LBR iiwa, posed but not solved, because seven
joints against three position coordinates is a count varietas refuses, and the
refusal is computed by the node rather than written down here.
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
        get_package_share_directory("varietas_demo"), "rviz", "branches.rviz"
    )

    refused_path = LaunchConfiguration("refused_urdf").perform(context)
    with open(refused_path, "r") as handle:
        refused_description = handle.read()

    return [
        Node(
            package="robot_state_publisher",
            executable="robot_state_publisher",
            output="log",
            parameters=[{"robot_description": robot_description}],
        ),
        # The arm varietas declines, posed from its own description. Its link
        # names do not collide with the solved arm's, so no frame prefix is
        # needed; the branch node publishes the static transform that stands it
        # beside the other.
        Node(
            package="robot_state_publisher",
            executable="robot_state_publisher",
            namespace="refused",
            output="log",
            parameters=[{"robot_description": refused_description}],
        ),
        Node(
            package="varietas_demo",
            executable="branch_node",
            output="screen",
            parameters=[
                {
                    "urdf": urdf_path,
                    "refused_urdf": refused_path,
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
        get_package_share_directory("varietas_demo"), "urdf",
        "anthropomorphic_offset_3r.urdf",
    )
    return LaunchDescription(
        [
            DeclareLaunchArgument("urdf", default_value=default_urdf),
            DeclareLaunchArgument(
                "refused_urdf",
                default_value=os.path.join(
                    get_package_share_directory("varietas_urdf"), "data", "lbr_iiwa14.urdf"
                ),
            ),
            DeclareLaunchArgument("period", default_value="24.0"),
            OpaqueFunction(function=_launch),
        ]
    )
