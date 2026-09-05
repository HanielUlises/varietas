============
varietas_ik
============

.. cpp:namespace:: varietas::ik

The connection between the two halves: it poses the inverse kinematics of a
chain over :math:`\Q(\p)` and returns exactly the
:cpp:struct:`varietas::codegen::parametric_solution` that ``emit`` consumes.
Namespace ``varietas::ik``.

Parametric inverse kinematics
=============================

``varietas/ik/parametric_ik.hpp``

.. cpp:function:: template<std::size_t N, std::size_t P> \
                  parametric_ik_result<N, P> \
                  parametric_position_ik(const chain<rational>& robot, \
                                         const std::array<std::size_t, P>& coordinates, \
                                         buchberger_statistics* statistics = nullptr)

   Poses the inverse kinematics with the target adjoined to the coefficient
   field rather than substituted into the ring, saturates away the half-angle
   loci, and returns the solved system.

.. cpp:struct:: template<std::size_t N, std::size_t P> parametric_ik_result

   ``status``, ``solution``, ``branches`` (:math:`\dim_k A`: the number of
   branches at a general pose, filled in whenever the quotient could be computed
   at all), and ``offending_coordinate``.

.. cpp:enum:: parametric_ik_status

   Most of these are settled by **counting**, before any Gröbner basis is
   attempted.

   .. list-table::
      :header-rows: 1
      :widths: 40 60

      * - Status
        - Why
      * - ``wrong_degrees_of_freedom``
        - The chain does not have the :math:`N` actuated joints the call was
          instantiated for.
      * - ``bad_coordinates``
        - The requested coordinates were not :math:`P` distinct indices below
          three.
      * - ``underdetermined``
        - :math:`P < N`. The residuals are :math:`P` polynomials in :math:`N`
          variables, so by Krull's height theorem every component has dimension
          at least :math:`N-P`, and saturation only *removes* components. The
          result is either the unit ideal or still positive-dimensional; there
          is nothing to emit either way.
      * - ``overdetermined``
        - :math:`P > N`. The positions an :math:`N`-joint arm reaches form a
          variety of dimension at most :math:`N`, and the parameters are
          transcendentals rather than a point of that image — so a general pose
          is out of reach and the ideal is the unit ideal.
      * - ``dropped_coordinate_is_not_identically_zero``
        - Leaving a coordinate out drops its equation, which is legitimate only
          when that equation says nothing — when the tool's coordinate there is
          identically zero as a polynomial, as the :math:`z` of a planar arm is.
          Checked rather than assumed.
      * - ``not_zero_dimensional``
        - The saturated ideal over :math:`\Q(\p)` has no finite solution set at
          a general pose. The usual cause is more joints than the coordinates
          constrain; two joints turning about one axis is another.

   The first two counting cases leave :math:`P = N` as the **only** arrangement
   that can produce a parametric solver — worth knowing before choosing
   ``--coords`` rather than after waiting for a Gröbner basis to say so.

Decoupled inverse kinematics
============================

``varietas/ik/decoupled_ik.hpp``

.. cpp:function:: template<std::size_t N> \
                  decoupled_solution<N> \
                  decoupled_position_ik(const chain<rational>& robot, \
                                        buchberger_statistics* statistics = nullptr)

   Sweeps the first joint out instead of solving for it, leaving a two-joint
   problem in two parameters. See :doc:`../guide/decoupling`.

.. cpp:struct:: sweep_frame

   Which coordinate plays which part once the first axis is known: ``axis``
   (:math:`k`, the direction the joint turns about), ``radial`` (:math:`a`, the
   radius of the reduced problem), ``swept`` (:math:`b`, identically zero on the
   reduced arm), ``reversed``, and ``plane()``.

.. cpp:struct:: template<std::size_t N> decoupled_solution

   ``status``, ``reduced_status``, ``frame``, ``reduced`` (a
   ``parametric_solution<N-1, 2>``, emitted exactly as any other),
   ``first_joint_name``, and ``branches`` — twice the reduced dimension, since
   each reduced solution occurs at the swept angle and again half a turn away
   with the radius negated.

.. cpp:enum:: decoupling_status

   ``ok``, ``wrong_degrees_of_freedom``, ``first_joint_not_revolute``,
   ``first_axis_not_a_coordinate_direction``, ``first_joint_is_displaced``,
   ``does_not_reduce``, ``reduced_problem_refused``. The table in
   :doc:`../guide/decoupling` explains each.

Emitting the decoupled solver
=============================

``varietas/ik/emit_decoupled.hpp``

.. cpp:function:: template<std::size_t N> \
                  std::string emit_decoupled(const decoupled_solution<N>& solution, \
                                             codegen::emit_options options = {})

   Emits the reduced solver as ``<name>_reduced`` and the whole-arm wrapper as
   ``<name>``, the latter through ``emit_options::epilogue``. Requires the
   ``eigen`` runtime, since the wrapper calls the reduced solver.

.. cpp:function:: template<std::size_t N> \
                  std::string decoupled_epilogue(const decoupled_solution<N>& solution, \
                                                 const std::string& reduced_name, \
                                                 const std::string& wrapper_name)

   The wrapper alone, as text.

Casting
=======

``varietas/ik/cast_chain.hpp`` — ``cast_vector3``,
``cast_matrix3``, ``cast_rigid_transform``,
``cast_joint``, ``cast_chain``, between coefficient fields.
