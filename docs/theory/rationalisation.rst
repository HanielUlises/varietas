=================================
Rationalisation and saturation
=================================

The transform that carries the map
==================================

:cpp:class:`varietas::rational_transform` carries a rigid transform whose
entries are rational functions of the joint variables. Under the half-angle
substitution, Rodrigues' formula for a revolute joint about a unit axis reads

.. math::

   R=\frac{1}{1+t^{2}}\Big[(1+t^{2})I+2t\,[u]_\times+2t^{2}[u]_\times^{2}\Big],

a matrix of quadratics over a single denominator.

The property the representation rests on is that this shape is **closed under
composition**: writing :math:`T=(1/D)[R\mid p]`, the product of two such has
rotation :math:`R_aR_b` and translation :math:`R_ap_b+p_aD_b` over
:math:`D_aD_b`, so the numerators multiply and the denominator exponents add.
Nothing else can arise.

A transform is therefore a numerator together with the exponent vector :math:`e`
of :math:`\prod_i(1+t_i^2)^{e_i}`, and **no common denominator ever has to be
computed, no gcd taken**.

What it deliberately is not is ``matrix3<Coeff>``: that type asks
``coefficient_traits`` for inverses, and polynomials are not invertible, so the
substitution needs its own type — which is the type system saying something true
rather than an inconvenience.

Clearing the denominators against a target gives the residuals, and
:cpp:func:`varietas::position_ideal_generators` and
:cpp:func:`varietas::pose_ideal_generators` hand them to saturation. Asking only
*where* the tool is, rather than how it is turned, is the classical problem for
an arm with fewer joints than :math:`SE(3)` has dimensions.

Saturation, and why it is not optional
======================================

:math:`I:h^\infty` is computed by Rabinowitsch's trick: adjoin a variable
:math:`y` and the generator :math:`1-yh`, so that in the enlarged ring the
variety is the part of :math:`\V(I)` where :math:`h\neq 0` with :math:`y`
recording :math:`1/h`, then eliminate :math:`y` and appeal to the Closure
Theorem.

The construction is therefore an elimination and uses the block order machinery
exactly as it stood — its first real application. Two details earn their keep:
the adjoined variable goes first so the block order eliminates it, and the
trailing block is ordered by the caller's own order, so what comes back is
already a Gröbner basis of the saturation under that order rather than something
needing a second Buchberger run. ``embed.hpp`` supplies the variable shifting,
re-sorting terms under the target order rather than copying a term list that was
sorted under a different one.

This is what removes the components the half-angle substitution invents.
Clearing the denominators :math:`1+t_i^2` attaches to the variety the loci where
they vanish, namely :math:`t_i=\pm i` — the images of the configurations
:math:`q_i=\pi` that the substitution sent to infinity. Over the reals they are
invisible, since :math:`1+t^2` never vanishes there, and **that is precisely why
they must be removed symbolically**: they are counted regardless, and the count
is the completeness certificate.

The tests exhibit both ways it goes wrong. In the mild case the spurious locus is
a few extra points and :math:`\dim_k A` reads 3 where the arm has one
configuration. In the severe case — the one that actually arises, where the
denominator divides more than one generator — the spurious locus is a whole
line, the ideal is not zero-dimensional at all, and the Finiteness Theorem
returns no verdict. Saturating is what makes the system *solvable* in the first
place, not merely what makes the count right.

The pipeline on the planar arm
==============================

.. figure:: ../figures/ik_branches.svg
   :width: 55%
   :alt: The two configurations of the inverse kinematics at one pose

   The two elbow configurations at a single prescribed tool position, both
   returned by the spectral solver. The pose is built by evaluating the rational
   map at a rational point, so the problem never leaves the rationals; the
   quotient has dimension two over the field, and that number — not the count of
   what was found — is the certificate that no branch is missing.

On the planar two-link arm the pipeline runs end to end, over :math:`\Q`
throughout. A target built by evaluating the rational map at
:math:`t=(1/2,1/3)` is rational, so nothing leaves the exact field; the
saturated ideal is zero-dimensional with :math:`\dim_k A=2`; the configuration
the target was built from satisfies every generator **exactly**, which over the
rationals is a structural test and not a residual below a threshold; and the
spectral solver returns both elbow branches, each of which puts the tool back on
the target.

The rational map itself is checked against KDL on the iiwa — seven joints, axes
turned by a right angle at every link, which is the composition the closure
property has to survive.

Unreachable two ways
====================

Two boundary cases are worth separating, because conflating them is how a solver
comes to lie.

A point **off the plane** of the arm is unreachable in the strong sense: the
tool has :math:`z=0` identically, a residual reads :math:`0=1`, the variety is
empty over :math:`\bar k` and the ideal is the whole ring — which is what the
weak Nullstellensatz detects.

A point **in the plane but beyond reach** is not that at all. The ideal is not
the unit ideal, and :math:`\dim_k A` still reads 2, because over :math:`\bar k`
the equations solve perfectly well with :math:`\cos q_2 = 23/2` and an angle
that is not real.

Unreachability is a statement about the **real** points of the variety, not
about the ideal, and it appears where the solver separates the real solutions
from the complex ones and nowhere earlier.
