========================
Building a chain by hand
========================

Everything downstream of the URDF front end is written against
:cpp:class:`varietas::chain`, so a chain typed out in C++, read from a DH table,
or produced by any other means enters the pipeline by the same door.

The pieces
==========

A chain is a sequence of joints, each carrying its type, its axis, and the fixed
placement of its frame as an element of :math:`SE(3)` over the coefficient
field, closed by a tool transform.

.. code-block:: cpp

   #include <varietas/kinematics/chain.hpp>
   #include <varietas/codegen/rational.hpp>

   using varietas::rational;
   using varietas::make_rational;
   using vec = varietas::vector3<rational>;
   using tf  = varietas::rigid_transform<rational>;

   varietas::chain<rational> arm("planar_2r");

   arm.add_joint(varietas::revolute_joint<rational>(
       "shoulder", vec::unit(2), tf::identity()));

   arm.add_joint(varietas::revolute_joint<rational>(
       "elbow", vec::unit(2),
       tf::translation_only(vec(make_rational(1), 0, 0))));

   arm.set_tool(tf::translation_only(vec(make_rational(1), 0, 0)));

:cpp:func:`varietas::revolute_joint`, :cpp:func:`varietas::prismatic_joint` and
:cpp:func:`varietas::fixed_joint` are the convenience constructors for the joints
a URDF actually contains, where the axis is a coordinate direction and the origin
a pure translation.

Rotations that are rational by construction
===========================================

Do **not** assemble a joint origin from roll–pitch–yaw angles and round the
result. Use either

* :cpp:func:`varietas::rotation_from_quaternion` — the Euler–Rodrigues matrix
  divided by the squared norm, exactly orthogonal of determinant one for any
  nonzero rational quaternion, no square root taken; or
* :cpp:func:`varietas::rotation_about_axis` — Rodrigues' formula in a
  cosine–sine pair rather than an angle.

A quarter turn about :math:`x` is ``rotation_from_quaternion(1, 1, 0, 0)``: a
quaternion is homogeneous, so the irrational :math:`(\cos 45^\circ, \sin
45^\circ, 0, 0)` and the integral :math:`(1,1,0,0)` are the same rotation. See
:doc:`../theory/exactness`.

Validate before solving
=======================

.. code-block:: cpp

   const auto diagnostic = arm.validate();
   if (!diagnostic.ok()) {
     std::printf("joint %zu: %s (defect %g)\n", diagnostic.joint_index,
                 varietas::to_string(diagnostic.status), diagnostic.defect);
     return 1;
   }

``validate`` checks orthogonality, properness, unit axes and limits, and returns
a :cpp:struct:`varietas::chain_diagnostic` naming **which** joint failed and
**how** — not a bare bool. Over an exact field the tolerance argument is ignored:
the only admissible orthogonality defect is zero.

Two indexings, kept apart
=========================

A fixed joint contributes no variable, so the index of a joint in the model and
the index of its variable in the polynomial ring are different numbers.
``chain::variable_of_joint`` and ``chain::joint_of_variable`` convert between
them, and they are kept explicit because generated code indexes by variable
while a robot model indexes by joint — conflating them is how a solution ends up
applied to the wrong axis.

``chain::degrees_of_freedom()`` is the :math:`N` the rest of the library is
templated on.

Fold the fixed joints first
===========================

``chain::fold_fixed_joints()`` returns the same chain with runs of fixed joints
multiplied into the origin of the joint that follows them, and any trailing run
folded into the tool frame. The two chains have identical forward kinematics by
associativity in :math:`SE(3)`, and the folded one is what the rationalisation
should be handed: every fixed joint left in place is a matrix product repeated in
every generator.

Then rationalise
================

.. code-block:: cpp

   #include <varietas/kinematics/rationalize.hpp>

   const auto folded = arm.fold_fixed_joints();
   const auto generators =
       varietas::position_ideal_generators<2, varietas::grevlex>(folded, target);
   const auto basis = varietas::groebner_basis(generators);
   const auto answer = varietas::solve_zero_dimensional(basis);

Or hand the chain straight to :cpp:func:`varietas::ik::parametric_position_ik` and
:cpp:func:`varietas::codegen::emit` to get a header instead of an answer
(:doc:`generated_headers`).

Working over ``double``
=======================

Every construction is templated on the field, and
:cpp:func:`varietas::chain_cast` / ``varietas::ik::cast_chain`` convert a
chain between them. :cpp:func:`varietas::forward_kinematics` and
:cpp:func:`varietas::link_frames` evaluate a ``chain<double>`` numerically,
which is what the tests compare the exact path against.
