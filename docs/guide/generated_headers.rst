=================
Generated headers
=================

``emit`` takes a system solved over :math:`\Q(\p)`, with the pose adjoined to
the coefficient field rather than to the polynomial ring so that one basis
answers every pose instead of one basis per pose, and writes a header.

What the header contains
========================

.. list-table::
   :header-rows: 0
   :widths: 34 66

   * - ``order_id``, ``order_name``
     - The monomial order the basis was computed under, so that a runtime
       assuming a different one is caught rather than silently answering a
       different question.
   * - ``num_unknowns``, ``num_parameters``
     - Sizes, as ``static constexpr``.
   * - ``dimension``
     - :math:`\dim_k A`: the number of solutions counted with multiplicity, and
       the size of every action matrix.
   * - ``one_index``
     - Where the monomial :math:`1` sits in the standard basis. The eigenvalue
       method divides by the eigenvector's component there.
   * - ``action_matrix(v, pose, out)``
     - The action matrix of unknown ``v`` at this pose, column-major.
   * - ``variable_coordinates(v, pose, out)``
     - The coordinates of the normal form of unknown ``v`` in the standard
       basis. A variable is usually **not** a standard monomial, having been
       reduced away, so its value at a point is this combination evaluated
       there rather than a coordinate read off.
   * - ``solve(pose, out, capacity, state, tol)``
     - Present only for the ``eigen`` runtime.

Two runtimes
============

``matrices_only``
   Includes ``<cstddef>`` and ``<cstdint>``, names nothing from this library,
   and leaves the eigenproblem to the caller. Drop it into a project that has
   never heard of Eigen.

``eigen``
   The above plus ``solve``, which builds a separating combination of the
   matrices, decomposes its transpose (left eigenvectors of a multiplication
   operator are the evaluation functionals at the points of the variety), and
   returns the real solutions. This is the default for ``urdf_codegen``;
   ``--matrices-only`` selects the other.

Calling it
==========

.. code-block:: cpp

   #include "arm_ik.hpp"

   using solver = varietas_generated::urdf_ik;

   const double pose[solver::num_parameters] = {0.4, 0.1};
   double out[8 * solver::num_unknowns];

   solver::status state{};
   const int found = solver::solve(pose, out, 8, &state);

``solve`` returns the number of **real** configurations written, or ``-1`` with
``state`` set to why. Real, because a joint angle that is complex is not a
configuration a manipulator can be commanded to; complex points are found and
then discarded.

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - ``status``
     - Meaning
   * - ``ok``
     - The set is complete.
   * - ``bad_pose``
     - A denominator vanished: this pose is off the chart the parametric basis
       describes.
   * - ``eigensolver_failed``
     - Eigen did not converge.
   * - ``deficient``
     - An eigenvector carried no point; the set is **not** certified complete.

The denominator guard
=====================

Every denominator is guarded, so a pose on the locus the parametric basis fails
to describe is **refused rather than answered with infinities**.

The guard does not compare against zero. A denominator that vanishes
mathematically almost never evaluates to ``0.0``; it evaluates to whatever the
cancellation between its terms leaves behind. On the planar 2R arm the pole is
the circle :math:`x^2+y^2+2x=0` (an artefact of the elimination rather than
anything the arm cannot reach) and at :math:`(-1.6, 0.8)`, which is on it
exactly over :math:`\Q`, the same expression in doubles comes to about
:math:`4\times10^{-16}`.

So each denominator is compared against **the sum of the magnitudes of the terms
that produced it**: total cancellation is recognised whether or not it landed on
zero, and a pose merely near the pole is still answered. Distinct denominators
are deduplicated, so the guard costs one comparison per pole rather than one per
matrix entry.

Rationals, not decimals
=======================

Constants are written as the quotient of two exact integer literals rather than
as decimals. A decimal rounds once here and again when the compiler parses it,
and cannot represent :math:`1/3` at all; handing the compiler the same rational
the offline computation held lets it produce the correctly rounded ``double`` in
one step.

Correctness of the emitter
==========================

Generated code is checked by being **compiled**: a program links the emitter
during the build, writes a header, and the test suite ``#include``\ s it, so
emitted text that does not parse is a build failure.

The solutions it returns are required to satisfy the original equations *and* to
agree with :cpp:func:`varietas::solve_zero_dimensional` run on the same system
with the pose substituted beforehand, which are two genuinely different
computations.

The same holds one level up: a pose substituted into the basis computed over
:math:`\Q(\p)` has to give the basis computed over :math:`\Q` with that pose
substituted first, standard monomials included. Two different Buchberger runs,
compared as a test rather than assumed.

``emit`` also checks
:cpp:func:`varietas::codegen::parametric_solution::is_well_formed` before
writing anything, because a header generated from an inconsistent solution
compiles perfectly and answers wrongly, which is the worst failure this code can
have.

Newton steps against the equations as posed
===========================================

A solution can carry the equations it came from, the residuals
:math:`\mathrm{numerator}_k(t) - \mathrm{denominator}(t)\,\mathrm{pose}_k`, and
both solve paths attach them. The header then carries them too, with their
Jacobian differentiated exactly, and ``solve()`` takes up to **two Newton
steps** from each point the eigenvalue method returns, keeping a step only when
it lowers the residual.

The reason is the eigenvectors. The eigenvalue method recovers a point from an
eigenvector of the separating form, and an eigenvector is computed only to
backward stability: near a pose where two solutions give that form nearly the
same value, the eigenvectors lose digits that the matrices never lost, and so do
the points read off them. Newton on the original equations recovers those digits
at the cost of a few polynomial evaluations. A point already at rounding level
comes back unchanged, and near a singular configuration, where the Jacobian
cannot be trusted, a step that would make things worse is not taken.

What it costs
=============

Measured rather than asserted, by ``doc/experiments/solver_cost.cpp``, on the
three-joint arm of :doc:`decoupling` with the header emitted for it:

.. list-table::
   :header-rows: 1
   :widths: 34 22 22 22

   * -
     - median
     - 99th
     - ratio
   * - ``branch_ik::solve``, target reachable
     - 1.5 µs
     - 1.9 µs
     - 1
   * - ``branch_ik::solve``, target out of reach
     - 0.8 µs
     - 1.4 µs
     - 0.5
   * - :cpp:func:`varietas::forward_kinematics`, same arm
     - 0.37 µs
     - 0.8 µs
     - 0.25
   * - damped least squares, one seed
     - 20 µs
     - 350 µs
     - 13

A solve costs about **four forward-kinematics evaluations**, which is the useful
way to hold it, and about two thirds of a million solves a second on one core.
The Newton steps are about 0.6 µs of that; without them a solve cost 0.9 µs,
and was less accurate in the tail. The generated solver is not the expensive
part of any control loop it is likely to sit in.

Generating the header costs **under a tenth of a second, once, offline**. It is
paid by the build, not by the caller.

Against a numerical solver
--------------------------

The comparison is not really about speed, though the speed is not close. A
damped least squares iteration from a random seed takes about thirteen times as
long as one generated solve, converges from only about **91%** of seeds, and
when it does converge returns **one** configuration: whichever one the seed fell
into, with no way to say how many others exist.

Recovering the whole solution set numerically means restarting it. Over three
hundred targets, taking the generated solver's answer as the roll of postures
that exist — certified by :math:`\dim_k A`, which is the point — it took about
**ten seeds per target** to find them all, missed a branch entirely on one
target inside a budget of sixty seeds, and cost some **420 µs per target against
1.5** for the single generated call.

So the generated solver is nearly three hundred times cheaper than the
numerical route for the answer the library actually promises, and unlike it,
returns a count that is a theorem rather than a hope.

Accuracy
--------

Over 71,658 returned configurations, each put back through the forward map and
compared against the target it was asked for:

* median :math:`2.3\times10^{-16}` m, which is the arithmetic's own floor;
* 99th centile :math:`7.5\times10^{-16}` m;
* worst :math:`1.6\times10^{-15}` m, and none above :math:`10^{-12}` m.

Before the Newton steps the worst was :math:`6.5\times10^{-8}` m, with about one
configuration in five thousand above :math:`10^{-12}` m, and the
`denominator guard`_ was the leading suspect. It was not the guard; it was the
eigenvectors, as above.

.. _denominator guard: #the-denominator-guard

The full three-joint header
---------------------------

The same arm has a full solver over :math:`\Q(x,y,z)` as well, written by
``urdf_codegen --reconstruct`` in about two and a half seconds, and
``doc/experiments/full_solver_cost.cpp`` sets the two side by side over twenty
thousand targets:

.. list-table::
   :header-rows: 1
   :widths: 40 30 30

   * -
     - full, :math:`\Q(x,y,z)`
     - decoupled
   * - solve, median
     - 6.2 µs
     - 1.4 µs
   * - residual, 99th centile
     - :math:`8.0\times10^{-16}` m
     - :math:`7.8\times10^{-16}` m
   * - residual, worst
     - :math:`2.0\times10^{-6}` m
     - :math:`1.5\times10^{-15}` m

They return the same configurations on 19,999 targets of 20,000. The one
exception lies within :math:`2\times10^{-4}` of the base axis, where a whole
circle of configurations is about to appear and the Jacobian cannot be trusted.

The full header also refuses a whole plane of targets, and the reason is the
half-angle substitution rather than the solve. A target with :math:`y = 0` is
reached, if at all, with the base at :math:`q_1 = 0` or :math:`q_1 = \pi`, and
:math:`t_1 = \tan(q_1/2)` sends the second to infinity: over :math:`\Q` the
quotient at such a pose has dimension two rather than four, so the parametric
basis has a pole along the plane, and ``solve()`` reports ``bad_pose`` there.
Just off the plane it answers correctly, the configurations with the base near
:math:`\pi` having :math:`t_1` of the order of :math:`1/y`. The decoupled
header recovers the base angle by an arctangent and has no such plane. **Where
an arm decouples, the decoupled header is the one to use**; the full one is for
the arms that do not.

The epilogue hook
=================

:cpp:struct:`varietas::codegen::emit_options` carries an ``epilogue``: verbatim
declarations placed inside the namespace after the struct. A solved system is not
always the whole answer. An arm whose first joint was swept out needs an
arctangent applied to what the header returns, and that arctangent belongs in the
same header as the matrices it accompanies. :doc:`decoupling` is the one user of
it today.
