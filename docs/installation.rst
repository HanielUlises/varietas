============
Installation
============

varietas targets **ROS 2 Humble** and **C++17**. Eigen and GMP are the only
non-ROS dependencies, and the algebraic half of the library needs no ROS at
all.

Dependencies
============

.. code-block:: sh

   sudo apt install libeigen3-dev libgmp-dev
   sudo apt install ros-humble-urdf ros-humble-kdl-parser \
                    ros-humble-robot-state-publisher ros-humble-rviz2

The ROS packages are required only by :doc:`guide/tools` (``varietas_urdf``)
and by the RViz demonstration (``varietas_demo``).

Building
========

.. code-block:: sh

   source /opt/ros/humble/setup.bash
   colcon build
   source install/setup.bash

A single package and what it depends on:

.. code-block:: sh

   colcon build --packages-up-to varietas_kinematics

Testing
=======

.. code-block:: sh

   colcon test
   colcon test-result --all --verbose

   colcon test --packages-select varietas_codegen

Generated code is checked by being compiled. A program links the emitter during
the build, writes a header, and the test suite ``#include``\ s it, so emitted
text that does not parse is a build failure rather than a test failure.

.. _without-ros:

Without ROS
===========

The algebra needs no ROS. ``varietas_core`` needs Eigen; ``varietas_codegen``
and ``varietas_kinematics`` additionally need GMP.

.. code-block:: sh

   g++ -std=c++17 -O2 example.cpp \
     -Ivarietas_core/include -Ivarietas_codegen/include \
     -Ivarietas_kinematics/include -Ivarietas_ik/include \
     -I/usr/include/eigen3 -lgmpxx -lgmp

A header that ``urdf_codegen`` produced needs neither: it includes
``<cstddef>`` and ``<cstdint>``, plus Eigen when it was emitted with the
solver.

.. code-block:: sh

   g++ -std=c++17 -O2 consumer.cpp -I. -I/usr/include/eigen3

``varietas_urdf`` and ``varietas_demo`` require ROS.

Building this documentation
===========================

.. code-block:: sh

   python3 -m pip install -r docs/requirements.txt
   sphinx-build -b html docs docs/_build/html

The published site at https://hanielulises.github.io/varietas/ is built and
deployed by ``.github/workflows/docs.yml`` on every push to ``main``. That
workflow builds with ``-W``, so a dead cross-reference or a malformed table
fails the check rather than shipping.
