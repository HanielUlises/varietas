======
Status
======

Version 0.2.0. What follows is what the library does, what it costs, and what it
refuses, the last being the largest part, and deliberately.

.. list-table::
   :header-rows: 1
   :widths: 62 38

   * - Construction
     - State
   * - Monomials, polynomials, orders
     - complete
   * - Division, Buchberger, reduced bases, membership
     - complete
   * - Dimension, saturation, elimination
     - complete
   * - Quotient algebra, action matrices, spectral solving
     - complete
   * - Exact rationals over GMP
     - complete
   * - Chains, validation, exact geometry
     - complete
   * - URDF front end with audit
     - complete
   * - Half-angle and trigonometric rationalisation
     - complete
   * - Workspace implicitization
     - complete
   * - Singular locus, dimension, workspace image
     - complete
   * - Code emission
     - complete
   * - URDF → header pipeline
     - complete
   * - First-joint decoupling
     - complete
   * - Six joints with a spherical wrist, full pose
     - complete
   * - Joint ranges in generated and runtime solvers
     - complete
   * - MoveIt kinematics plugin
     - complete, six joints with a spherical wrist
   * - Pose inverse kinematics at a given target
     - complete, to five joints
   * - Modular gcd over :math:`\Q(\p)`
     - complete, certified
   * - Parametric solve reconstructed from fixed poses
     - complete, checked at rational poses
   * - Newton refinement in generated solvers
     - complete
   * - Factorisation over :math:`\Q`
     - complete, any number of variables
   * - Decomposition along the factors
     - not begun

The unfinished row is :doc:`roadmap`.

What refusal looks like
=======================

The library refuses more than it accepts, which is the point, and the counts
settle most of it before any Gröbner basis is attempted.

Fewer pose coordinates than unknowns cannot give a finite solution set:
:math:`P` polynomials in :math:`N` variables cut out components of dimension at
least :math:`N-P`, and saturation only removes components. More coordinates than
unknowns cannot give a nonempty one: the positions an :math:`N`-joint arm
reaches form a variety of dimension at most :math:`N`, the parameters are
transcendentals rather than a point of that variety, and a general pose is
simply out of reach, and the ideal is the unit ideal.

So :math:`P=N` is the only arrangement that can produce a solver at all, and
both other cases are rejected by counting. A coordinate left out of the
parameter list is separately checked to be identically zero before its equation
is dropped, since dropping an equation that constrains the joints would silently
answer about a larger variety. What counting cannot catch, two joints turning
about one axis for instance, the quotient dimension does.

Two things that were once wrong
===============================

Worth recording, because the tests that pin them down read oddly without the
history.

**The denominator guard compared against zero.** That is the wrong question in
floating point: a denominator that vanishes mathematically almost never
evaluates to ``0.0``, it evaluates to whatever the cancellation between its
terms leaves behind. On the planar 2R arm the pole is the circle
:math:`x^2+y^2+2x=0`, an artefact of the elimination rather than anything the
arm cannot reach, and at :math:`(-1.6, 0.8)`, which is on it exactly over
:math:`\Q`, the same expression in doubles comes to about
:math:`4\times10^{-16}`. The guard did not fire, the matrices were formed by
dividing by that, and one of the two returned branches did not reach the target,
with neither the count nor the status saying so. The guard now compares each
denominator against the sum of the magnitudes of the terms that produced it.

**Monomial exponents were held in a** ``uint8_t``. Parameter polynomials over a
function field reach that ceiling (degree 254 in one variable was observed in
an underdetermined system), after which a product wrapped to a different
monomial and a polynomial that divided another silently stopped dividing it.
Exponents are sixteen bits now, and the product asserts rather than wraps, which
is the half that matters.

.. _the-cost:

The third thing is scale
========================

It is the one that decided what this pipeline is for. **The cost is not in the
arm but in the number of parameters adjoined**, and in what the solve carries on
the way.

.. list-table::
   :header-rows: 1
   :widths: 40 20 20 20

   * - System
     - Field
     - Symbolic
     - Reconstructed
   * - Planar 2R
     - :math:`\Q(x,y)`
     - ~10 ms
     - ~5 ms
   * - Anthropomorphic 3R, fixed pose
     - :math:`\Q`
     - 4 ms
     - —
   * - Anthropomorphic 3R, decoupled
     - :math:`\Q(r,z)`
     - ~20 ms
     - —
   * - Anthropomorphic 3R
     - :math:`\Q(x,y,z)`
     - ~20 s
     - 0.25 s
   * - Base off its own axis (no sweep)
     - :math:`\Q(x,y,z)`
     - ~35 s
     - 0.25 s
   * - Elbow displaced sideways (no sweep)
     - :math:`\Q(x,y,z)`
     - not in 54 min, 29 GB
     - 0.55 s
   * - Shoulder displaced (the demonstration arm)
     - :math:`\Q(x,y,z)`
     - not after 60 min
     - 2.5 s
   * - Three pairwise skew axes (no sweep)
     - :math:`\Q(x,y,z)`
     - not after 60 min
     - 47 s

The 3R arm at a fixed pose returns a basis of six elements and
:math:`\dim_k A = 4`, exactly the four branches such an arm is known to have.
The arm is not the difficulty.

The first difficulty was cancellation
=====================================

``rational_function`` normalises after every coefficient operation, and with the
subresultant gcd the cancellation was what the run was made of, about **86%**
of the time, at a cost per call growing as about the **fourth power** of the
number of terms in the operands. ``doc/parametric_cost.pdf`` reports the
experiment that settled what to do about it. Repeating the three-parameter
solve over :math:`\F_p(x,y,z)`, where coefficient arithmetic is a machine
multiplication, did not complete either, and the choice of field was worth a
bounded factor of thirteen. Computing the same remainder sequence over several
primes could not have been enough.

The gcd is now Brown's dense modular algorithm, certified by exact division
over :math:`\mathbb{Z}` (:doc:`api/codegen`). On the same operands its cost
grows as about the 1.5th power of the size, and at 455 terms it takes
**under 5 ms against 91 s**. The anthropomorphic arm over :math:`\Q(x,y,z)`, which
produced nothing in fifteen minutes, finishes in about twenty seconds, with the
gcd half of it.

The second was swell
====================

On arms with offsets the symbolic solve still does not finish, and the gcd is
no longer why. With the elbow displaced a quarter of a link sideways the
parameter polynomials Buchberger carries reach **sixteen thousand terms**, the
gcd is a quarter of the time, and the rest is arithmetic on those polynomials.
The action matrices the solve would produce have entries of fourteen terms at most.

The reconstruction (:cpp:func:`varietas::ik::reconstructed_position_ik`) never
forms the intermediate functions. It solves the arm at many poses over prime
fields and recovers each matrix entry from its values; the table above is the
result. It has no final exact test, so it is **checked exactly** at rational
poses instead, and where both finish the two agree entry for entry. The
symbolic solve stays the default of ``urdf_codegen`` because it is certified;
``--reconstruct`` is the route for the arms it cannot reach.

The way to solve a three-joint arm that admits it is still to sweep the base
joint out (:doc:`guide/decoupling`), which is exact, takes tens of
milliseconds, and its header is several times cheaper to call.

Against a real robot
====================

Run against the KUKA iiwa in ``varietas_urdf/test/data``, the pipeline mostly
declines, and the way it declines is the useful part.

The chain recovers exactly, with the audit moving no joint by more than
:math:`5\times10^{-12}` radians, and then has **seven joints**, which is four
more than a tool position can constrain, so ``urdf_codegen`` says so and stops.

Truncating the chain does not rescue it either. Taking the tip to be
``lbr_iiwa_link_3`` gives three joints, but that link's frame is the third
joint's own frame, so the tool sits on the axis that joint turns about and
cannot be moved by it. The position problem is a two-joint problem wearing three
joints, the reduced system comes out positive-dimensional, and the dimension
check catches it in about half a second. **There is no link boundary on this arm
where a well-posed three-joint positioning problem appears.**

That is a fair summary of the present reach of the parametric path, which is a tool
for positioning subsystems of up to three joints, not for a seven-axis
manipulator, and the thing it does well is refuse promptly and say which of the
three reasons applies.

The fixed-pose path
===================

Five joints is where it stops being interactive and six is where it stops; see
:doc:`guide/pose_ik` for the timings.
