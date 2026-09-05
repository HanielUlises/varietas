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

The epilogue hook
=================

:cpp:struct:`varietas::codegen::emit_options` carries an ``epilogue``: verbatim
declarations placed inside the namespace after the struct. A solved system is not
always the whole answer. An arm whose first joint was swept out needs an
arctangent applied to what the header returns, and that arctangent belongs in the
same header as the matrices it accompanies. :doc:`decoupling` is the one user of
it today.
