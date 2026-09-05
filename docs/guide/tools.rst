==================
Command line tools
==================

All three live in ``varietas_urdf`` and are run with ``ros2 run``.

.. _urdf_report:

``urdf_report``: the audit
==========================

.. code-block:: sh

   ros2 run varietas_urdf urdf_report <file.urdf> [tip_link] [root_link]

Reads the model, recovers the chain over :math:`\Q`, and prints joint by joint
which placements were recovered without moving, which were moved onto the exact
geometry the decimals approximated, and by how far. Nothing is solved.

Run this first on any model you have not used before. A chain that cannot be
recovered is refused here, with the reason and the joint responsible, rather
than approximated silently. See :doc:`../theory/exactness`.

.. _urdf_codegen:

``urdf_codegen``: URDF to header
=================================

.. code-block:: sh

   ros2 run varietas_urdf urdf_codegen <file.urdf> <output.hpp> \
       [--tip L] [--root L] [--coords xy|xz|yz|xyz] \
       [--name N] [--namespace NS] [--decouple] [--matrices-only]

.. list-table::
   :header-rows: 1
   :widths: 26 74

   * - Option
     - Meaning
   * - ``--tip L``
     - Tip link. Defaults to the sole leaf; required if the model branches.
   * - ``--root L``
     - Root link. Defaults to the model's root.
   * - ``--coords xy``
     - Which position coordinates are adjoined as parameters. A coordinate left
       out is checked to be **identically zero** before its equation is dropped
       Dropping an equation that constrains the joints would silently answer
       about a larger variety.
   * - ``--name N``
     - Name of the generated struct. Default ``urdf_ik``.
   * - ``--namespace NS``
     - Namespace of the generated header. Default ``varietas_generated``.
   * - ``--decouple``
     - Sweep the base joint out instead of adjoining it: three joints, and no
       ``--coords``. See :doc:`decoupling`.
   * - ``--matrices-only``
     - Emit the action matrices and variable coordinates but not ``solve``, so
       the header needs only ``<cstddef>`` and ``<cstdint>``. The default
       runtime is ``eigen``.

The counts settle most of it before any Gröbner basis is attempted; see
:doc:`../status` for why :math:`P=N` is the only arrangement that can produce a
solver at all.

.. _urdf_solve:

``urdf_solve``: one pose, exactly
==================================

.. code-block:: sh

   ros2 run varietas_urdf urdf_solve <file.urdf> --xyz X Y Z \
       [--rpy R P Y | --quat X Y Z W] \
       [--tip L] [--root L] [--position-only] [--denominator D]

.. list-table::
   :header-rows: 1
   :widths: 26 74

   * - Option
     - Meaning
   * - ``--xyz X Y Z``
     - The target position. Required.
   * - ``--rpy R P Y``
     - Target orientation as roll–pitch–yaw.
   * - ``--quat X Y Z W``
     - Target orientation as a quaternion.
   * - ``--position-only``
     - Solve for the position and leave the orientation free.
   * - ``--denominator D``
     - Denominator bound for snapping the target onto a nearby exact pose.
       Default :math:`10^6`.

Up to five joints in practice. See :doc:`pose_ik`.

The demonstration
=================

.. code-block:: sh

   ros2 launch varietas_demo sweep.launch.py urdf:=<file.urdf> period:=12.0

Sweeps a target and drives the chain with the recovered kinematics in RViz,
publishing the tool pose computed from the exact chain alongside the one
``robot_state_publisher`` derives from the file.

Test models
===========

``varietas_urdf/test/data`` ships three: ``planar_2r.urdf``,
``anthropomorphic_3r.urdf`` and ``lbr_iiwa14.urdf``. The first two are solvable
by their respective paths; the iiwa is there because the interesting thing it
does is get **refused**, promptly and with a reason (:doc:`../status`).
