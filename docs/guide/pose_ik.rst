==================
One pose at a time
==================

Orientation cannot be adjoined. A general pose is six parameters — twelve matrix
entries are not independent, and a general point of :math:`\A^{12}` is not a
rigid motion at all — and two adjoined parameters is the working limit, so there
is no parametric solver for a full pose and there is not going to be one at this
size.

**Giving up the parameters is what buys orientation back.** With the target a
constant of :math:`\Q` rather than a parameter of :math:`\Q(\p)`, the twelve
equations of a pose are no harder for Buchberger than the three of a position,
and ``varietas_kinematics`` has posed them since the beginning.
:ref:`urdf_solve` is that path as a command.

What it costs
=============

The growth is steep, and it is the whole story. On a chain of unit links, a full
pose takes

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Joints
     - Time for a full pose
   * - 3
     - about 10 ms
   * - 4
     - about a sixth of a second
   * - 5
     - about 9 seconds
   * - 6
     - nothing after eight minutes

**Five joints is where this stops being interactive, and six is where it
stops.**

The pose has to be exact
========================

This is not a formality. A target assembled by rounding sines and cosines into
rationals carries denominators near :math:`2^{52}`, and a run over coefficients
that size does not finish — the first attempt at a six-joint pose here failed
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
four postures it genuinely has — the same four the decoupled generated solver
returns, by an entirely different computation, which is the cross-check the two
paths give each other.
