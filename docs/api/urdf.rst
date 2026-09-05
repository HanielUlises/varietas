==============
varietas_urdf
==============

.. cpp:namespace:: varietas::urdf_import

A robot description turned into a chain over :math:`\Q`, or refused with the
joint responsible named. Namespace ``varietas::urdf_import``. This is the one
package of the pipeline that requires ROS.

Rational approximation
======================

``varietas/urdf/rational_approximation.hpp``

.. cpp:function:: rational limit_denominator(const rational& value, const mpz_class& bound)

   The best rational approximation with denominator at most ``bound``, by
   continued fractions.

.. cpp:function:: rational rationalize(double value, long max_denominator = 1000000)

   ``0.1575`` comes back as :math:`63/400` — the decimal the author wrote, not
   the binary approximation the file stores.

.. cpp:function:: double nearest_double(const rational& value)

   Rounds to nearest. GMP's ``get_d`` truncates towards zero, which is a whole
   unit in the last place and always in the same direction.

.. cpp:function:: double rotation_distance(const matrix3<double>& a, const matrix3<double>& b)

   The chordal distance :math:`\|A-B\|_F = 2\sqrt{2}\sin(\theta/2)`. The
   textbook :math:`\arccos((\operatorname{tr} R - 1)/2)` is useless precisely
   where it matters: near the identity the trace differs from three by the
   *square* of the angle, so a deviation of :math:`10^{-12}` vanishes into
   rounding and the recovery is reported as exact when it is not.

.. cpp:struct:: scalar_snap

   The exact value and how far it moved.

.. cpp:function:: scalar_snap snap_scalar(double value, long max_denominator = 1000000)

.. cpp:struct:: rotation_snap

   The exact rotation and its chordal deviation.

.. cpp:function:: rotation_snap snap_rotation(double qx, double qy, double qz, double qw, \
                                              long max_denominator = 1000000)

   The projective search: divide the quaternion through by its largest entry,
   then approximate each of the four by a rational of bounded denominator. A
   quarter turn becomes :math:`(1,1,0,0)` exactly. See
   :doc:`../theory/exactness`.

Import
======

``varietas/urdf/urdf_chain.hpp``

.. cpp:enum:: import_status

   ``ok``, ``tip_link_not_found``, ``root_link_not_found``,
   ``root_not_on_chain``, ``unsupported_joint_type`` (floating and planar joints
   move in more than one degree of freedom and have no single axis to substitute
   for), ``axis_not_rational_unit``, ``rotation_deviation_exceeded``,
   ``translation_deviation_exceeded``, ``no_actuated_joints``,
   ``chain_invalid`` — the last meaning the importer, not the file, is at fault.

.. cpp:struct:: import_options

   ``max_denominator`` (default :math:`10^6`: enough for every decimal literal a
   URDF is likely to hold, and every multiple of a right angle exactly),
   ``rotation_tolerance`` and ``translation_tolerance`` (default
   :math:`10^{-9}`: loose enough to absorb a truncated :math:`\pi`, tight enough
   that a real misalignment is refused).

.. cpp:struct:: joint_recovery

   ``name``, ``rotation_deviation``, ``translation_deviation``, and ``unmoved``
   — the recovered placement is the stated one to the last bit. A joint whose
   angles were written as truncated decimals is *not* of this kind: it moves, by
   the truncation.

.. cpp:struct:: import_report

   ``status``, ``detail``, ``joints``, ``max_rotation_deviation``,
   ``max_translation_deviation``, ``ok()``, and
   ``indistinguishable_from_file()``.

.. cpp:function:: import_report collect_joints(const ::urdf::ModelInterface& model, ...)

.. cpp:function:: import_report chain_from_model(const ::urdf::ModelInterface& model, \
                                                 const std::string& root, \
                                                 const std::string& tip, \
                                                 chain<rational>& out, \
                                                 const import_options& options = {})

   The whole front end: walk the chain, snap every placement, audit what it
   cost, validate the result.

.. cpp:function:: std::string sole_tip_link(const ::urdf::ModelInterface& model)

   Empty if the model branches, which is when ``--tip`` becomes required.

Executables
===========

``urdf_report``, ``urdf_codegen`` and ``urdf_solve`` — see :doc:`../guide/tools`.
