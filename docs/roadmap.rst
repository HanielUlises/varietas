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

The other half of the lifted basis
==================================

**Why.** Fixed-pose bases over :math:`\Q` are lifted from :math:`\F_p` now
(:doc:`api/codegen`), and that is what brought six joints within reach of
the fixed-pose path. The check over :math:`\Q` proves that the lifted ideal
contains the true one, and not that it is no larger. The reconstruction has
the same gap and closes it with exact checks at rational poses, which are now
lifted bases themselves unless asked otherwise.

**What is needed.** Each element of the lifted basis written over the
generators, with cofactors that can be multiplied out over :math:`\Q`. The
cofactors could be lifted as the basis is, from F4 runs that record them, and
checking them is multiplication rather than reduction. Or a degree bound for
those cofactors, after which membership is one linear system.

Faster images
=============

**Why.** On six joints one image is one to three seconds, and the
coefficients need ten to thirty primes. Images are computed four at a time,
since more contend for memory and gain nothing on the laptop it was measured
on.

**What is needed.** The images of one ideal modulo different primes share
their shape, so the matrices F4 built for the first prime can be rebuilt for
the next with new coefficients and without the pair selection, the symbolic
preprocessing, or the rows that reduced to zero, which on the six-joint arm
are more than half of them. The same trace serves the reconstruction, whose
samples share their shape too. The row operations themselves would gain from
being written for the machine: delayed reduction modulo :math:`p` across
several rows, and rows stored by blocks.

The linear systems of the reconstruction
========================================

**Why.** With the sampling and the checks cheap, the arm with three pairwise
skew axes spends three quarters of its 6.6 s in the linear systems that turn
values into coefficients, a dense elimination per matrix entry.

**What is needed.** The systems for different entries share their matrix
except for one diagonal scaling, and could share most of the elimination. A
sparse interpolation, Zippel's or Ben-Or and Tiwari's, would need fewer
samples as well.

Irreducible components
======================

**Why.** :cpp:func:`varietas::decompose` finds its own splittings now, and a
singular locus comes back as pieces rather than as one ideal the caller has to
probe: on the anthropomorphic arm, the straight elbow, the folded elbow, and the
tool on the base axis. But a piece is only a variety none of whose basis
elements factors over :math:`\Q`, and that is weaker than irreducible in two
ways. A factorisation over an extension of :math:`\Q` is not seen, so
:math:`y^2 - 2x^2` is one piece though it is two planes; and reducibility that
no single basis element shows is not seen either.

**What is needed.** A prime decomposition: the minimal primes of the ideal,
by the method of Gianni, Trager and Zacharias or of Eisenbud, Huneke and
Vasconcelos, each reducing to a zero-dimensional problem in general position
and to factorisation over an algebraic extension of :math:`\Q`. The
factorisation over :math:`\Q` that both need is in place.

**What the factorisation still lacks.** The leading coefficient in the main
variable is imposed whole on both halves of each division rather than shared
out by Wang's method, so the lifting runs to the degree of :math:`l f` rather
than of :math:`f`. The univariate recombination is exponential in the worst
case, which the polynomials of Swinnerton-Dyer reach; van Hoeij's lattice
reduction is the cure.

Not on the roadmap
==================

**A parametric solver for a full pose.** A general pose is six parameters:
twelve matrix entries are not independent, and a general point of :math:`\A^{12}`
is not a rigid motion at all. The parametric path reaches a position, three
parameters, and a full pose is a different problem rather than a larger one. The fixed-pose path (:doc:`guide/pose_ik`) is what buys orientation
back, and it reaches six joints.

**Radicals.** varietas cannot take one, and would not want to where it can: the
pinched torus keeps a :math:`z^2` that records the order of contact the set alone
forgets (:doc:`theory/singularities`).

**Semialgebraic reach.** The reachable set of an arm has a boundary, boundaries
are inequalities, and no ideal expresses one. What the library returns is the
Zariski closure, and it says so (:doc:`theory/workspace`).
