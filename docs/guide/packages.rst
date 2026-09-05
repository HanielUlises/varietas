============
The packages
============

.. list-table::
   :header-rows: 1
   :widths: 24 46 30

   * - Package
     - Role
     - Depends on
   * - ``varietas_core``
     - Monomials, polynomials, orders, ideals, quotient algebra, solving
     - Eigen
   * - ``varietas_codegen``
     - Exact rationals over GMP; the offline field; code emission
     - ``varietas_core``, GMP
   * - ``varietas_kinematics``
     - Chains, rationalisation, workspace, singularities
     - ``varietas_core``
   * - ``varietas_ik``
     - Inverse kinematics posed over :math:`\Q(\p)`, ready to emit
     - ``varietas_kinematics``, ``varietas_codegen``
   * - ``varietas_urdf``
     - URDF → exact chain over :math:`\Q`, with an audit; the command line tools
     - ``varietas_ik``, ``urdf``
   * - ``varietas_demo``
     - RViz demonstration of the recovered chain
     - ``varietas_urdf``, ``rclcpp``

``varietas_core``, ``varietas_codegen``, ``varietas_kinematics`` and
``varietas_ik`` are **header-only interface targets**. Only ``varietas_urdf``
and ``varietas_demo`` require ROS.

Templated on the field
======================

``varietas_core`` is templated on the coefficient type so that the same code
runs over ``double`` at runtime and over an exact rational type in the offline
generator. Nothing in it changes between the two, and the two are exercised by
the same test bodies instantiated twice. See :doc:`../theory/exactness`.

The layering that matters
=========================

``varietas_kinematics`` **knows nothing of URDF**, which is the point. The
rationalisation, the emitter, the implicitization and the singular locus are all
written against :cpp:class:`varietas::chain`, so a URDF front end arrives later
as a translation into it rather than as a dependency of the pipeline, and
hand-written chains and DH tables enter by the same door.
