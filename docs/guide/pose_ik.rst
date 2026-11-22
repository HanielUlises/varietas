==================
One pose at a time
==================

Orientation cannot be adjoined. A general pose is six parameters, since twelve
matrix entries are not independent and a general point of :math:`\A^{12}` is not
a rigid motion at all, and the parametric path reaches a position, three
parameters. So there is no parametric solver for a full pose, and none is
claimed.

**Giving up the parameters is what buys orientation back.** With the target a
constant of :math:`\Q` rather than a parameter of :math:`\Q(\p)`, the twelve
equations of a pose are no harder for Buchberger than the three of a position,
and ``varietas_kinematics`` has posed them since the beginning.
:ref:`urdf_solve` is that path as a command.

What it costs
=============

Buchberger over :math:`\Q` grows steeply with the number of joints, and on a
chain of unit links it had produced nothing after eight minutes at six. Two
things changed that.

**The equations are split at the middle of the chain**
(:cpp:func:`split_pose_residuals`): the first half of the joints on one side,
the inverse of the second half on the other, which is again a chain. The
saturated ideal is the same, and on six joints a basis over :math:`\F_p`
costs a second or two rather than ten or sixty. (Buchberger over :math:`\Q`
does worse with them, so ``--buchberger`` keeps the equations of the whole
chain.)

**The basis over** :math:`\Q` **is lifted from its images modulo primes**
(:cpp:func:`lifted_pose_ideal_generators`), which carry one machine word per
coefficient where Buchberger over :math:`\Q` carries the swell. It is
checked over :math:`\Q` to be a Gröbner basis of an ideal containing the
true one, and is equal to it with high probability; ``--buchberger`` computes
it exactly instead, at the old cost.

On the industrial arm of the tests, cut after its first three, four and five
joints and complete at six, at the pose of a rational configuration:

.. list-table::
   :header-rows: 1
   :widths: 16 42 42

   * - Joints
     - Buchberger over :math:`\Q`, whole chain
     - Lifted, split chain
   * - 3
     - 4 ms
     - 2 ms
   * - 4
     - 0.52 s
     - 5 ms
   * - 5
     - 68 s
     - 29 ms
   * - 6
     - not observed to finish
     - 20 to 30 s, eight solutions

At six joints the time is the images: nineteen primes for coefficients of
262 bits, each image one F4 run of a few seconds, four at a time. A six-joint
arm without a spherical wrist is no harder here than one with: an arm of the
Universal Robots shape takes about sixteen seconds and has its eight
solutions too.

The pose has to be exact
========================

This is not a formality. A target assembled by rounding sines and cosines into
rationals carries denominators near :math:`2^{52}`, and a run over coefficients
that size does not finish. The first attempt at a six-joint pose here failed
that way before it failed for any interesting reason.

So the target is snapped onto a nearby exact pose by the same quaternion
argument the URDF front end uses (:doc:`../theory/exactness`), and **the
deviation is reported rather than hidden**. Right angles survive that snapping
exactly, which is what makes it usable rather than merely defensible. The
``--denominator`` option sets the bound.

Too few joints is not harmless over-specification
=================================================

Asking for a full pose from an arm with too few joints is not a harmless
over-specification. The poses a three-joint arm can strike form a
three-dimensional subvariety of a six-dimensional space, so a target moved even
slightly off it is **unreachable in the strict sense**: the ideal becomes the
unit ideal and the answer is no configuration at all, rather than a
configuration that nearly works. ``urdf_solve`` says so in those words.

Asking such an arm for a **position** instead (``--position-only``) returns the
four postures it genuinely has: the same four the decoupled generated solver
returns, by an entirely different computation, which is the cross-check the two
paths give each other.
