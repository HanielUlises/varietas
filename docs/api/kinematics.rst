====================
varietas_kinematics
====================

.. cpp:namespace:: varietas

The structure the algebra consumes, and the three constructions built on it.
Depends on ``varietas_core`` alone; it knows nothing of URDF.

Geometry
========

``varietas/kinematics/rigid_transform.hpp``

.. cpp:class:: template<class Coeff> vector3

   ``zero()``, ``unit(i)``, ``operator[]``, ``squared_norm()``, ``is_zero()``,
   and the usual arithmetic.

.. cpp:class:: template<class Coeff> matrix3

   ``identity()``, ``from_rows(entries)``, ``operator()(i, j)``,
   ``transpose()``, ``determinant()``, and products with vectors and matrices.

.. cpp:class:: template<class Coeff> rigid_transform

   ``identity()``, ``translation_only(t)``, ``rotation_only(r)``,
   ``rotation()``, ``translation()``, ``apply(p)``, ``inverse()``, and
   composition.

.. cpp:function:: template<class Coeff> \
                  matrix3<Coeff> rotation_from_quaternion(const Coeff& w, const Coeff& x, \
                                                          const Coeff& y, const Coeff& z)

   The Euler–Rodrigues matrix divided by the squared norm: exactly orthogonal of
   determinant one for **any** nonzero rational quaternion, with no square root.
   This is how exact rotations are built.

.. cpp:function:: template<class Coeff> \
                  matrix3<Coeff> rotation_about_axis(const vector3<Coeff>& axis, \
                                                     const Coeff& c, const Coeff& s)

   Rodrigues' formula in a cosine–sine pair rather than an angle.

.. cpp:function:: template<class Coeff> double orthogonality_defect(const matrix3<Coeff>& r)

Chains
======

``varietas/kinematics/chain.hpp``

.. cpp:enum:: joint_type

   ``fixed``, ``revolute``, ``prismatic``. A revolute joint contributes one
   variable :math:`t=\tan(q/2)`, a prismatic joint its displacement directly, a
   fixed joint none.

.. cpp:struct:: template<class Coeff> joint

   ``name``, ``type``, ``axis``, ``origin``, ``has_limits``, ``lower``,
   ``upper``, ``is_actuated()``. Limits are ``double``: they are inequalities,
   they play no part in any ideal, and they are consulted only when solutions
   are filtered down to the reachable ones.

.. cpp:enum:: chain_status

   ``ok``, ``origin_not_orthogonal``, ``origin_improper``, ``axis_degenerate``,
   ``axis_not_unit``, ``limits_inverted``, ``no_actuated_joints``.
   ``to_string`` gives a sentence, not a token.

.. cpp:struct:: chain_diagnostic

   ``status``, ``joint_index``, ``defect``. A named failure with the joint it
   belongs to, in preference to a bare bool.

.. cpp:class:: template<class Coeff> chain

   .. cpp:function:: chain& add_joint(joint<Coeff> j)
   .. cpp:function:: void set_tool(rigid_transform<Coeff> t)
   .. cpp:function:: std::size_t degrees_of_freedom() const noexcept

      The :math:`N` the rest of the library is templated on.

   .. cpp:function:: std::size_t variable_of_joint(std::size_t joint_index) const
   .. cpp:function:: std::size_t joint_of_variable(std::size_t variable_index) const

      Kept explicit: generated code indexes by variable, a robot model by
      joint, and conflating them is how a solution ends up on the wrong axis.

   .. cpp:function:: chain_diagnostic validate(double tolerance = 1e-12) const

      The tolerance is ignored over an exact field, where the only admissible
      defect is zero.

   .. cpp:function:: rigid_transform<Coeff> trailing_transform() const
   .. cpp:function:: chain fold_fixed_joints() const

      Runs of fixed joints multiplied into the joint that follows. Identical
      forward kinematics; far fewer repeated matrix products in every generator.

.. cpp:function:: template<class Coeff> joint<Coeff> revolute_joint(std::string name, const vector3<Coeff>& axis, const rigid_transform<Coeff>& origin)
.. cpp:function:: template<class Coeff> joint<Coeff> prismatic_joint(std::string name, const vector3<Coeff>& axis, const rigid_transform<Coeff>& origin)
.. cpp:function:: template<class Coeff> joint<Coeff> fixed_joint(std::string name, const rigid_transform<Coeff>& origin)

Numerical evaluation
====================

``varietas/kinematics/evaluate.hpp``

.. cpp:function:: template<class To, class From> chain<To> chain_cast(const chain<From>& source)
.. cpp:function:: rigid_transform<double> joint_displacement(const joint<double>& j, double value)
.. cpp:function:: std::vector<rigid_transform<double>> link_frames(const chain<double>& robot, const std::vector<double>& values)
.. cpp:function:: rigid_transform<double> forward_kinematics(const chain<double>& robot, const std::vector<double>& values)

Half-angle rationalisation
==========================

``varietas/kinematics/rationalize.hpp``

.. cpp:class:: template<class Coeff, std::size_t N, class Order> rational_transform

   A rigid transform whose entries are rational functions of the joint
   variables, held as a numerator together with the exponent vector of
   :math:`\prod_i (1+t_i^2)^{e_i}`. The shape is closed under composition, so
   **no common denominator is ever computed and no gcd taken**. See
   :doc:`../theory/rationalisation`.

.. cpp:function:: template<std::size_t N, class Order, class Coeff> \
                  rational_transform<Coeff, N, Order> rational_forward_kinematics(const chain<Coeff>& robot)

.. cpp:function:: template<std::size_t N, class Order, class Coeff> \
                  std::vector<polynomial<Coeff, N, Order>> pose_residuals(...)
.. cpp:function:: template<std::size_t N, class Order, class Coeff> \
                  std::vector<polynomial<Coeff, N, Order>> position_residuals(...)
.. cpp:function:: template<std::size_t N, class Order, class Coeff> \
                  std::vector<polynomial<Coeff, N, Order>> half_angle_multipliers(const chain<Coeff>& robot)
.. cpp:function:: template<std::size_t N, class Order, class Coeff> \
                  std::vector<polynomial<Coeff, N, Order>> kinematic_ideal_generators(...)
.. cpp:function:: template<std::size_t N, class Order, class Coeff> \
                  std::vector<polynomial<Coeff, N, Order>> pose_ideal_generators(...)
.. cpp:function:: template<std::size_t N, class Order, class Coeff> \
                  std::vector<polynomial<Coeff, N, Order>> position_ideal_generators(...)

   The residuals, saturated against the half-angle denominators.

.. cpp:function:: double angle_from_variable(double t)
.. cpp:function:: double variable_from_angle(double q)

   :math:`q = 2\arctan t` and its inverse.

Trigonometric rationalisation
=============================

``varietas/kinematics/trigonometric.hpp``

.. cpp:class:: template<class Coeff, std::size_t N, class Order> trigonometric_transform

   Two variables per joint, :math:`c_i` and :math:`s_i`, and no denominator
   field at all.

.. cpp:function:: template<class Coeff> std::size_t trigonometric_variable_count(const chain<Coeff>& robot)
.. cpp:function:: template<class Coeff> std::vector<std::size_t> cosine_variable_indices(const chain<Coeff>& robot)
.. cpp:function:: template<std::size_t N, class Order, class Coeff> \
                  trigonometric_transform<Coeff, N, Order> trigonometric_forward_kinematics(const chain<Coeff>& robot)
.. cpp:function:: template<std::size_t N, class Order, class Coeff> \
                  std::vector<polynomial<Coeff, N, Order>> circle_relations(const chain<Coeff>& robot)

   :math:`c_i^2+s_i^2-1`, on the same footing as the pose residuals: a genuine
   equation of the problem rather than damage repair.

.. cpp:function:: template<class Coeff> std::vector<Coeff> trigonometric_point(const chain<Coeff>& robot, const std::vector<Coeff>& angles)
.. cpp:function:: double angle_from_cosine_sine(double c, double s)

Workspace
=========

``varietas/kinematics/workspace.hpp``

.. cpp:function:: template<std::size_t V, class Coeff> \
                  std::vector<polynomial<Coeff, 3, grevlex>> workspace_relations(const chain<Coeff>& robot)

   The implicit equations of the Zariski closure of the reachable set, computed
   on the trigonometric ring, a factor of a thousand faster than the half-angle
   route on the torus arm, and provably the same ideal.

.. cpp:function:: template<std::size_t N, class Coeff> \
                  std::vector<polynomial<Coeff, 3, grevlex>> workspace_relations_half_angle(const chain<Coeff>& robot)

   Kept so that the two can be compared rather than trusted.

.. cpp:function:: template<class Coeff> \
                  bool workspace_is_dense(const std::vector<polynomial<Coeff, 3, grevlex>>& relations)

   Reports the full-dimensional case, where the closure is everything and the
   elimination ideal returned only what held identically, rather than leaving
   an empty result to be misread as a failure.

Also ``workspace_residuals``, ``workspace_layout``,
``trigonometric_workspace_generators``, ``trigonometric_workspace_layout`` for
callers wanting the graph rather than its projection.

Singularities
=============

``varietas/kinematics/singular.hpp``

.. cpp:function:: std::vector<std::size_t> position_rows()
.. cpp:function:: std::vector<std::size_t> orientation_rows()
.. cpp:function:: std::vector<std::size_t> pose_rows()
.. cpp:function:: std::vector<std::size_t> planar_pose_rows()

   Which rows of the Jacobian the minors are taken from. An arm is singular
   **for a task**, and asking for the wrong one is how a singularity analysis
   comes to disagree with the robot.

.. cpp:function:: template<std::size_t V, class Order, class Coeff> \
                  dense_matrix<polynomial<Coeff, V, Order>> trigonometric_jacobian(const chain<Coeff>& robot)

   The geometric Jacobian, columns :math:`[\,a_i\times(p-p_i);\,a_i\,]`,
   polynomial in :math:`c` and :math:`s` as it stands. Checked against central
   differences of the numerical forward kinematics, which share no code with it.

.. cpp:function:: template<std::size_t V, class Order, class Coeff> \
                  std::vector<polynomial<Coeff, V, Order>> singular_generators(const chain<Coeff>& robot, const std::vector<std::size_t>& rows)
.. cpp:function:: template<std::size_t V, class Coeff> \
                  std::vector<polynomial<Coeff, V, grevlex>> singular_ideal(const chain<Coeff>& robot, const std::vector<std::size_t>& rows)
.. cpp:function:: template<std::size_t V, class Coeff> \
                  affine_dimension<V> singular_dimension(const chain<Coeff>& robot, const std::vector<std::size_t>& rows)
.. cpp:function:: template<std::size_t V, class Coeff> \
                  std::vector<polynomial<Coeff, 3, grevlex>> singular_workspace_relations(const chain<Coeff>& robot, const std::vector<std::size_t>& rows)

   The singular locus pushed forward into the workspace: the surface the arm
   cannot cross without losing rank.
