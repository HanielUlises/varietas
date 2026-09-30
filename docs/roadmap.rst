=======
Roadmap
=======

What is missing is named rather than hedged. None of it is a matter of tuning;
each is a piece of computer algebra that has to be built.

The modular gcd over :math:`\Q(\p)` that used to head this page is done
(:doc:`api/codegen`), and so is the reconstruction that took over where it
stopped being enough (:doc:`status`).

A certificate for the reconstruction
====================================

**Why.** :cpp:func:`varietas::ik::reconstructed_position_ik` is checked exactly
at rational poses, and a wrong answer that survives those checks would need the
error to vanish at every one of them; but a check at points is evidence, not a
proof, and the symbolic solve it replaces was a proof.

**What is needed.** An exact test of the reconstructed matrices as rational
functions. The candidates are a degree bound on the entries, from which a
finite number of exact point checks becomes a proof by the usual count of the
zeros of a polynomial, or a direct verification that the reconstructed action
matrices commute and annihilate the generators over :math:`\Q(\p)`, which is
arithmetic on functions of tens of terms rather than of thousands.

A second chart for the half-angle substitution
==============================================

**Why.** :math:`t = \tan(q/2)` sends :math:`q = \pi` to infinity, and a full
three-joint solver therefore refuses every target that can only be reached, in
one of its families, with a joint at exactly :math:`\pi`; for a base that yaws
about :math:`z` that is the whole plane :math:`y = 0`
(:doc:`guide/generated_headers`).

**What is needed.** A second solve under the substitution
:math:`t = \tan((q - q_0)/2)` for a fixed rational :math:`q_0`, emitted
alongside the first, with the runtime choosing whichever chart is further from
its poles.

A faster fixed-pose basis over a prime field
============================================

**Why.** Reconstruction spends most of its time in the black box, one
Gröbner basis per sample: about 17 ms each on the arm with its axes in general
position, which is 60% of the 47 s that arm takes.

**What is needed.** A completion algorithm of the F4 family, which reduces
many pairs at once by row reduction over :math:`\F_p`, where a row operation is
a vector of machine words. Since the samples share their shape, the sequence of
reductions from one sample could also be replayed on the next without the
pair selection that found it.

Factorisation over :math:`\Q`
=============================

**Why.** :cpp:func:`varietas::split_along` is exact and exhaustive, but the
caller chooses the divisor :math:`h`. It separates a branch the geometry
suggests; it does not *find* the branches.

**What is needed.** Multivariate factorisation over :math:`\Q`. It is what would
turn ``split_along`` into a **decomposition proper**, a singular locus reported
as its irreducible components rather than as one ideal the caller has to probe.

**Status.** Not begun. It is a project of its own.

Not on the roadmap
==================

**A parametric solver for a full pose.** A general pose is six parameters:
twelve matrix entries are not independent, and a general point of :math:`\A^{12}`
is not a rigid motion at all. The parametric path reaches a position, three
parameters, and a full pose is a different problem rather than a larger one. The fixed-pose path (:doc:`guide/pose_ik`) is what buys orientation
back, and it stops at five joints.

**Radicals.** varietas cannot take one, and would not want to where it can: the
pinched torus keeps a :math:`z^2` that records the order of contact the set alone
forgets (:doc:`theory/singularities`).

**Semialgebraic reach.** The reachable set of an arm has a boundary, boundaries
are inequalities, and no ideal expresses one. What the library returns is the
Zariski closure, and it says so (:doc:`theory/workspace`).
