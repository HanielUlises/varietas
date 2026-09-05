=======
Roadmap
=======

Two constructions are missing, and both are named rather than hedged. Neither is
a matter of tuning; each is a piece of computer algebra that has to be built.

A modular gcd over :math:`\Q(\p)`
=================================

**Why.** About 86% of a parametric solve is spent in
:cpp:func:`varietas::polynomial_gcd`, and the cost per call grows as roughly the
fourth power of the number of terms in its operands (:ref:`the-cost`). This is
what holds the parametric path to two adjoined parameters.

**What is there now.** A subresultant remainder sequence, with a cheap
coprimality test bolted to the front: specialise every parameter but one in a
small prime field, take a univariate gcd on machine integers, and skip the exact
computation when the answer is constant. It is evidence rather than proof, which
is admissible because it only decides whether to *skip* a cancellation — a
skipped one costs size, never correctness.

It buys about 29 ms down to 25 ms on the reduced two-joint problem, and the
three-parameter system still produces no answer.

**What is needed.** A modular gcd proper: **evaluation, interpolation and
rational reconstruction**, rather than a remainder sequence with a filter in
front of it. The measurement in ``doc/parametric_cost.pdf`` is specifically
designed to rule out the cheaper alternative: repeating the three-parameter
solve over :math:`\F_p(x,y,z)`, where coefficient arithmetic is a single machine
multiplication, does not complete either, and the choice of field is worth only
a bounded factor of thirteen. **Computing the same remainder sequence over
several primes is therefore not enough.**

**What it would change.** Three adjoined parameters, which is to say a
three-joint positioning solver without decoupling, and therefore arms that do
not admit the sweep.

Factorisation over :math:`\Q`
=============================

**Why.** :cpp:func:`varietas::split_along` is exact and exhaustive, but the
caller chooses the divisor :math:`h`. It separates a branch the geometry
suggests; it does not *find* the branches.

**What is needed.** Multivariate factorisation over :math:`\Q`. It is what would
turn ``split_along`` into a **decomposition proper** — a singular locus reported
as its irreducible components rather than as one ideal the caller has to probe.

**Status.** Not begun. It is a project of its own.

Not on the roadmap
==================

**A parametric solver for a full pose.** A general pose is six parameters —
twelve matrix entries are not independent, and a general point of :math:`\A^{12}`
is not a rigid motion at all. Two adjoined parameters is the working limit, so
there is no parametric solver for a full pose and there is not going to be one
at this size. The fixed-pose path (:doc:`guide/pose_ik`) is what buys orientation
back, and it stops at five joints.

**Radicals.** varietas cannot take one, and would not want to where it can: the
pinched torus keeps a :math:`z^2` that records the order of contact the set alone
forgets (:doc:`theory/singularities`).

**Semialgebraic reach.** The reachable set of an arm has a boundary, boundaries
are inequalities, and no ideal expresses one. What the library returns is the
Zariski closure, and it says so (:doc:`theory/workspace`).
