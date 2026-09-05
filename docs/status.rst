======
Status
======

Version 0.1.0. What follows is what the library does, what it costs, and what it
refuses — the last being the largest part, and deliberately.

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
   * - Pose inverse kinematics at a given target
     - complete, to five joints
   * - Modular gcd over :math:`\Q(\p)`
     - coprimality test only
   * - Factorisation over :math:`\Q`
     - not begun

The two unfinished rows are :doc:`roadmap`.

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
simply out of reach — the ideal is the unit ideal.

So :math:`P=N` is the only arrangement that can produce a solver at all, and
both other cases are rejected by counting. A coordinate left out of the
parameter list is separately checked to be identically zero before its equation
is dropped, since dropping an equation that constrains the joints would silently
answer about a larger variety. What counting cannot catch — two joints turning
about one axis, say — the quotient dimension does.

Two things that were once wrong
===============================

Worth recording, because the tests that pin them down read oddly without the
history.

**The denominator guard compared against zero.** That is the wrong question in
floating point: a denominator that vanishes mathematically almost never
evaluates to ``0.0``, it evaluates to whatever the cancellation between its
terms leaves behind. On the planar 2R arm the pole is the circle
:math:`x^2+y^2+2x=0` — an artefact of the elimination rather than anything the
arm cannot reach — and at :math:`(-1.6, 0.8)`, which is on it exactly over
:math:`\Q`, the same expression in doubles comes to about
:math:`4\times10^{-16}`. The guard did not fire, the matrices were formed by
dividing by that, and one of the two returned branches did not reach the target,
with neither the count nor the status saying so. The guard now compares each
denominator against the sum of the magnitudes of the terms that produced it.

**Monomial exponents were held in a** ``uint8_t``. Parameter polynomials over a
function field reach that ceiling — degree 254 in one variable was observed in
an underdetermined system — after which a product wrapped to a different
monomial and a polynomial that divided another silently stopped dividing it.
Exponents are sixteen bits now, and the product asserts rather than wraps, which
is the half that matters.

.. _the-cost:

The third thing is scale
========================

It is the one that decides what this pipeline is for. **The cost is not in the
arm but in the number of parameters adjoined.**

.. list-table::
   :header-rows: 1
   :widths: 45 25 30

   * - System
     - Field
     - Time
   * - Planar 2R
     - :math:`\Q(x,y)`
     - ~30 ms
   * - Anthropomorphic 3R, fixed pose
     - :math:`\Q`
     - 4 ms
   * - Anthropomorphic 3R
     - :math:`\Q(x,y,z)`
     - nothing after 15 min
   * - Anthropomorphic 3R, decoupled
     - :math:`\Q(r,z)`
     - ~45 ms

The 3R arm at a fixed pose returns a basis of six elements and
:math:`\dim_k A = 4` — exactly the four branches such an arm is known to have.
The arm is not the difficulty.

The difficulty is cancellation
==============================

``rational_function`` normalises after every coefficient operation, and the
polynomial gcd that cancellation needs is what the run is made of. Sampled over
three minutes:

* about **86%** of the time is spent inside ``polynomial_gcd``;
* the cost per call grows sharply as the parameter polynomials do — the first
  twenty-seven thousand calls take five seconds between them, while a later
  sixteen hundred take sixty-eight.

Cancellation cannot simply be dropped, since without it the entries grow past
any use.

Most of those calls return 1, so normalisation now asks a cheap question first
(the coprimality test described in :doc:`api/codegen`). It is worth being plain
about what that buys, which is not much: on the reduced two-joint problem it
takes about twenty-nine milliseconds down to about twenty-five, and the
three-parameter system still produces no answer. **The bottleneck is real but
the fast path only avoids the cheap half of it.**

``doc/parametric_cost.pdf`` in the repository reports an experiment that narrows
what a proper construction would have to be: repeating the three-parameter solve
over :math:`\F_p(x,y,z)`, where coefficient arithmetic is a machine
multiplication, does not complete either, and an isolated measurement puts the
subresultant gcd at about the **fourth power** of the number of terms in its
operands, with the choice of field worth a bounded factor of thirteen.

Computing the same remainder sequence over several primes is therefore not
enough; only a gcd built on evaluation and interpolation addresses the cost that
was measured.

Until it is, **two adjoined parameters remains the working limit** — and since
:math:`P=N`, that is to say two joints. The way to solve a three-joint arm today
is to sweep one of them out rather than adjoin it (:doc:`guide/decoupling`).

Against a real robot
====================

Run against the KUKA iiwa in ``varietas_urdf/test/data``, the pipeline mostly
declines, and the way it declines is the useful part.

The chain recovers exactly — the audit moves no joint by more than
:math:`5\times10^{-12}` radians — and then has **seven joints**, which is four
more than a tool position can constrain, so ``urdf_codegen`` says so and stops.

Truncating the chain does not rescue it either: taking the tip to be
``lbr_iiwa_link_3`` gives three joints, but that link's frame is the third
joint's own frame, so the tool sits on the axis that joint turns about and
cannot be moved by it. The position problem is a two-joint problem wearing three
joints, the reduced system comes out positive-dimensional, and the dimension
check catches it in about half a second. **There is no link boundary on this arm
where a well-posed three-joint positioning problem appears.**

That is a fair summary of the present reach of the parametric path: it is a tool
for small subsystems, not for a seven-axis manipulator, and the thing it does
well is refuse promptly and say which of the three reasons applies.

The fixed-pose path
===================

Five joints is where it stops being interactive and six is where it stops; see
:doc:`guide/pose_ik` for the timings.
