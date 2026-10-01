========================
MoveIt kinematics plugin
========================

``varietas_moveit`` is a kinematics plugin for MoveIt. For a planning group
that is a six-joint arm whose last three axes meet, which is most industrial
arms, it replaces the default solver with one that returns every configuration
reaching the pose. Nothing is generated or compiled for the robot. The plugin
reads the URDF MoveIt already holds when it is loaded.

Using it
========

In the MoveIt configuration's ``kinematics.yaml``, name the plugin for the
group.

.. code-block:: yaml

   manipulator:
     kinematics_solver: varietas_moveit/VarietasKinematicsPlugin

The group must be a chain from the arm's base to its flange or tool, as the
Setup Assistant writes it. Timeouts and search resolution are accepted and
ignored, since the plugin does not search.

What happens when it loads
==========================

The chain between the group's base and tip is recovered exactly from the URDF,
split at its wrist centre (:doc:`spherical_wrist`), and the arm that places the
centre is solved, by the decoupling where its base sweeps a plane and by
reconstruction where it does not. On the industrial arm in
``varietas_urdf/test/data`` this takes about a second and a half and happens
once. The joint ranges come from MoveIt rather than the URDF, so ranges set in
``joint_limits.yaml`` are respected.

A group the decomposition does not apply to is refused at load time, with the
reason in the log, for example that the arm has seven joints or that its wrist
axes do not meet.

What a query does
=================

A query evaluates the solver at the pose. It returns every configuration that
reaches it, each angle moved by whole turns into its joint's range, and a
configuration that no turn brings inside is dropped. Where a joint can reach an
angle two ways, a joint wider than a turn or one without limits, the plugin
picks the representative nearest the seed, so the arm is not sent round a full
turn. The configurations are ordered by distance from the seed.

The seed therefore decides which configuration comes first, and never whether
one is found. ``getPositionIK`` returns the first, the overload that returns
several returns them all, and ``searchPositionIK`` offers them in order to the
solution callback, so that a configuration rejected for a collision gives way
to the next rather than ending the search. Consistency limits are applied the
same way. A pose with no configuration inside the ranges is reported at once.

Against the default solver
==========================

``compare_with_kdl``, built with the package's tests, loads this plugin and
MoveIt's KDL plugin for the same group on the industrial arm and asks both for
the same 2000 reachable poses from the same random seeds, KDL with a timeout of
50 ms.

.. list-table::
   :header-rows: 1
   :widths: 22 18 20 20 20

   * -
     - solved
     - median
     - 99th centile
     - worst position
   * - varietas
     - 100.0%
     - 11.5 µs
     - 20.5 µs
     - :math:`1.5\times10^{-14}` m
   * - KDL
     - 97.8%
     - 1476 µs
     - 52 191 µs
     - :math:`1.0\times10^{-5}` m

KDL iterates from the seed and stops at the first configuration it converges
to, or at the timeout, which is where its 99th centile sits. The plugin
evaluates a closed solution, so it costs the same whatever the seed, and it can
return all eight configurations where KDL returns one.
