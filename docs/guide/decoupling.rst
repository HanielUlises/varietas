==========================
Decoupling the first joint
==========================

Two adjoined parameters is the working limit of the parametric path
(:doc:`../status`), and since the counts force :math:`P=N`, that is to say two
joints. Two joints is less of a restriction than it sounds, because **the joint
that costs the most is often the one that need not be adjoined at all**.

The construction
================

:cpp:func:`varietas::ik::decoupled_position_ik` sweeps the first joint out
instead of solving for it.

If the base turns about a fixed axis and the rest of the arm holds the tool in a
plane containing that axis, then turning the base sweeps that plane around it:
the tool's position is a radius and a height in the plane together with the
angle the plane has been turned through, and that angle is an **arctangent of
the target** rather than an eigenvalue. What is left is a two-joint problem in
two parameters.

The anthropomorphic arm, base yawing about :math:`z` with shoulder and elbow
pitching about :math:`y`, reduces this way in about **45 ms**, against the
fifteen minutes that produced nothing when all three coordinates were adjoined.

.. code-block:: sh

   ros2 run varietas_urdf urdf_codegen arm.urdf arm_ik.hpp --decouple

goes from the URDF to the reduced header in about a second.

Each solution of the reduced problem is a configuration of the arm **twice
over**, facing the target and reversed half a turn away, which recovers the
four postures such an arm is known to have. That is the same count that solving
the whole system over :math:`\Q` at a fixed pose reports, by an entirely
different computation.

Whether an arm admits it
========================

The decomposition is a question about the arm, so it is asked of the arm.

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - ``decoupling_status``
     - Why it was refused
   * - ``wrong_degrees_of_freedom``
     - Fewer than two actuated joints, or a count that disagrees with the one
       the call was instantiated for.
   * - ``first_joint_not_revolute``
     - The first joint has to turn; the whole construction is about the plane
       it sweeps.
   * - ``first_axis_not_a_coordinate_direction``
     - An arbitrary axis can be decoupled too, by rotating the frame first,
       but that rotation is not in general rational and would leave the exact
       field.
   * - ``first_joint_is_displaced``
     - The placement must commute with its own rotation, which for a pure
       translation means lying along the axis. Anything else moves the axis
       itself and there is no fixed plane to sweep.
   * - ``does_not_reduce``
     - The rest of the arm does not hold the tool in a plane through the axis.
       The tool's coordinate along the swept direction must be identically zero
       with the base held still.
   * - ``reduced_problem_refused``
     - The reduced two-joint problem was itself refused; the reason comes from
       there.

A **planar** arm fails the last of these, since its base turns in the plane the
arm already works in and there is nothing to sweep, so it is refused rather than
reduced.

What is generated
=================

A solver for the whole arm, in two parts:

``urdf_ik_reduced``
   The reduced two-joint problem, an ordinary generated header
   (:doc:`generated_headers`).

``urdf_ik``
   The wrapper, carrying the arctangent, the pairing of each reduced solution
   with the two turns of the plane, and the inversion of the half-angle
   substitution.

.. code-block:: cpp

   using solver = varietas_generated::urdf_ik;   // num_joints, max_configurations

   const double target[3] = {0.4, 0.1, 0.3};
   double out[solver::max_configurations * solver::num_joints];

   solver::status state{};
   const int found = solver::solve(target, out, solver::max_configurations, &state);
   // out + k * num_joints is one configuration, the swept joint first.

.. important::

   The wrapper returns **joint angles in radians**, not the ring's variables,
   and deliberately so: the reversed family puts the base angle near :math:`\pi`
   routinely, and :math:`\tan(q/2)` is unbounded there. Returning the angle
   costs one arctangent per joint and is defined everywhere.
