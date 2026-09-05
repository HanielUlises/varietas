==================
The singular locus
==================

A configuration is singular when the differential of the forward kinematics
drops rank there, and the arm loses, instantaneously, the ability to move its
tool in some direction. The standard treatment evaluates the Jacobian at a
configuration and reports its smallest singular value; that is a measurement at
a point, and no number of such measurements says what the singular set *is*, how
many pieces it has, or where they go in the workspace — because those are
questions about a variety.

Rank deficiency is polynomial
=============================

Rank is not a polynomial condition — it is not even a continuous function of the
entries — but rank *deficiency* is: a matrix has rank below :math:`k` exactly
when all its :math:`k\times k` minors vanish.

:cpp:func:`varietas::maximal_minors` computes those minors over the polynomial
ring by Laplace expansion memoised on the set of remaining columns, which is
:math:`2^k k` multiplications rather than :math:`k!` and, more to the point, is
the arrangement that **never divides** — Gaussian elimination needs to invert a
pivot and Bareiss needs exact division, and polynomials offer neither.

Adjoined to the circle relations, the maximal minors generate the ideal of the
singular locus (:cpp:func:`varietas::singular_ideal`).

Why the differential is available
=================================

This is the trigonometric formulation earning its keep a second time. The
differential of a map restricted to a variety is not the Jacobian of the
polynomials that cut the variety out; over the half-angle ring one would
differentiate the cleared numerators and then correct for the denominators
multiplied through — a quotient rule in every entry.

Here the correction is unnecessary, and the reason is exact: along
:math:`c^2+s^2=1` the tangent at :math:`(c,s)` is :math:`(-s,c)`, which is
precisely :math:`(\mathrm{d}c/\mathrm{d}q,\ \mathrm{d}s/\mathrm{d}q)`, so
differentiating the polynomial map along the constraint reproduces
:math:`\partial/\partial q`. The geometric Jacobian — columns
:math:`[\,a_i\times(p-p_i);\,a_i\,]` straight out of the textbook — is already
polynomial in :math:`c` and :math:`s`. It is checked against central differences
of the numerical forward kinematics, which share no code with it.

Singular *for a task*
=====================

Which rows of that Jacobian the minors are taken from is the caller's choice,
because an arm is singular **for a task**, and asking for the wrong one is how a
singularity analysis comes to disagree with the robot. varietas offers
:cpp:func:`varietas::position_rows`, :cpp:func:`varietas::orientation_rows`,
:cpp:func:`varietas::pose_rows` and :cpp:func:`varietas::planar_pose_rows`.

A planar three-link arm is singular *everywhere* for the three-dimensional
position task, correctly: it can never move its tool out of its plane, so the
rank is at most two at every configuration. For the planar pose task it can
span, the same arm is singular on the surface :math:`\sin q_2 = 0` — the elbow
alone, the wrist adding no rank. Both are in the tests, adjacent, because the
pair is the point.

The two-link arm, end to end
============================

.. figure:: ../figures/singular_image.svg
   :width: 60%
   :alt: The image of the singular locus of the planar two-link arm

   The singular image of the planar two-link arm. Elimination returns :math:`z`
   together with :math:`x(x^2+y^2-4)` and :math:`y(x^2+y^2-4)`: the circle of
   radius :math:`l_1+l_2` and the point at :math:`|l_1-l_2|`, which is to say
   the outer boundary of the reachable disc and the hole at its centre. The tint
   is the reachable set, which is semialgebraic and which no ideal cuts out.

Its singular ideal is :math:`\langle c_i^2+s_i^2-1,\ s_2\rangle` — and the minor
itself is :math:`(c_1^2+s_1^2)s_2`, the textbook :math:`l_1l_2\sin q_2`
multiplied by something the determinant has no way of knowing is one, which
reduces to :math:`s_2` modulo the circle relations and not before.

The locus is one dimensional: a circle's worth of configurations with the
shoulder free and the elbow pinned. Splitting along :math:`c_2-1` separates the
straight elbow from the folded one exactly. Eliminating the joint variables
carries it into the workspace as the boundary above — the singular image is the
workspace boundary, computed rather than reasoned about.

Two negative results
====================

**No real singularity, but not the unit ideal.** The torus arm has no real
singularity — its tool is at distance :math:`2+\cos q_2` from the first axis and
that never vanishes — yet its singular ideal is not the unit ideal, because over
:math:`\bar k` the equations solve at :math:`\cos q_2=-2`,
:math:`\sin^2 q_2=-3`. Eliminated into the workspace this reads
:math:`x^2+y^2=0`, :math:`z^2=-3`, which has no real point, and the arm is
therefore nowhere singular; but the ideal could not say so, exactly as it could
not decide unreachability, and for the same reason: **reality is a property of
points, not of ideals**.

.. figure:: ../figures/pinched_torus.svg
   :width: 70%
   :alt: The workspace closure of the arm whose offset equals its tool length

   The same construction on the arm whose offset equals its tool length: the
   inner circle of the torus has closed to a point, and every singular
   configuration maps to that pinch. The elimination returns :math:`x`,
   :math:`y` and :math:`z^2` — not :math:`z`.

**A non-radical image.** The pinched torus, where the offset equals the tool
length, is singular on the circle :math:`c_2=-1`, all of which maps to the
origin, and the elimination returns :math:`x`, :math:`y` and :math:`z^2` — not
:math:`z`. The variety is the single point the Closure Theorem promises, but the
ideal is not radical, and the square is not noise: the arm reaches the pinch
tangentially in :math:`z`, and the elimination ideal has kept the order of
contact that the set alone forgot. Taking a radical would discard it, and
varietas cannot take one anyway.

What is not claimed
===================

A primary decomposition. :cpp:func:`varietas::split_along` is exact and
exhaustive, but the caller chooses :math:`h`, and an algorithm that finds its own
splittings needs multivariate factorisation over :math:`\Q` — a project of its
own, and not yet begun (:doc:`../roadmap`).

What the library offers is the ideal, its dimension, its image in the workspace,
and the ability to separate a branch on a divisor the geometry suggests.
