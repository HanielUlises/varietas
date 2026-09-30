"""Every configuration of an arm nothing else can solve, beside one iteration.

The arm is general_3r, whose three axes are pairwise skew. The decoupling
refuses it and the symbolic solve over Q(x, y, z) does not finish. Its solver
was reconstructed from fixed poses during the build. Two copies of the arm
follow one moving target. On the left, every configuration the generated
solver returns, each branch in its own colour. On the right, a damped least
squares iteration warm-started from its previous answer, drawn in the colour
of the branch it happens to be on. Run with

    ros2 launch varietas_demo solve_all.launch.py

`trace:=/path/file.csv` writes every tick to a file, which the recording
script turns into the captions of the published figure.
"""

import os

from ament_index_python.packages import get_package_share_directory
from launch import LaunchDescription
from launch.actions import DeclareLaunchArgument
from launch.substitutions import LaunchConfiguration
from launch_ros.actions import Node


def generate_launch_description():
    share = get_package_share_directory("varietas_demo")
    return LaunchDescription(
        [
            DeclareLaunchArgument("period", default_value="24.0"),
            DeclareLaunchArgument("labels", default_value="true"),
            DeclareLaunchArgument("trace", default_value=""),
            DeclareLaunchArgument("rviz", default_value=os.path.join(share, "rviz", "solve_all.rviz")),
            Node(
                package="varietas_demo",
                executable="solve_all_node",
                output="screen",
                parameters=[
                    {
                        "urdf": os.path.join(share, "urdf", "general_3r.urdf"),
                        "period": LaunchConfiguration("period"),
                        "labels": LaunchConfiguration("labels"),
                        "trace": LaunchConfiguration("trace"),
                    }
                ],
            ),
            Node(
                package="rviz2",
                executable="rviz2",
                arguments=["-d", LaunchConfiguration("rviz")],
                output="log",
            ),
        ]
    )
