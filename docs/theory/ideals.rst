=====================
Ideals and dimension
=====================

Division and Gröbner bases
==========================

varietas implements the multivariate division algorithm with quotients and
remainder (:cpp:func:`varietas::divide`, :cpp:func:`varietas::normal_form`) and
Buchberger's algorithm with the normal selection strategy and both Buchberger
criteria (:cpp:func:`varietas::buchberger`). The number of critical pairs each
criterion discarded is reported in :cpp:struct:`varietas::buchberger_statistics`
rather than being hidden — the two counts are the difference between a run that
finishes and one that does not.

:cpp:func:`varietas::groebner_basis` minimalises and fully reduces, giving the
*unique* reduced Gröbner basis for the order. Uniqueness is what makes the
cross-checks in the test suite possible: two different computations that ought
to produce the same ideal can be compared element by element.

The order lives in the type
===========================

Polynomials carry their monomial order as a template parameter, so two
polynomials written under different orders cannot be added. The available
orders are

* :cpp:class:`varietas::lex` — lexicographic,
* :cpp:class:`varietas::grlex` — graded lexicographic,
* :cpp:class:`varietas::grevlex` — graded reverse lexicographic,
* ``block_order<Split, First, Second>`` — for elimination,
* ``weight_order<Weights, TieBreak>``.

Each records an :cpp:enum:`varietas::order_id`, which a generated header stores
so that a mismatch between the order used offline and the one assumed at runtime
is caught rather than silently producing a wrong solution set.

Membership, emptiness, dimension
================================

:cpp:func:`varietas::is_member` decides ideal membership by normal form.
:cpp:func:`varietas::is_unit_ideal` detects the unit ideal — that is, an
**empty variety**: a target no configuration reaches.

:cpp:func:`varietas::ideal_dimension` is :math:`\dim \V(I)` read off the leading
terms. A subset of the variables is independent modulo :math:`I` when no leading
monomial is supported inside it, and the dimension is the size of the largest
such subset.

The empty variety falls out of the same statement rather than needing a case:
the unit ideal's leading monomial :math:`1` is supported inside every subset, so
it is reported as **empty** rather than as dimension zero — which is what a
*point* would be. That distinction is not pedantry; it is the difference between
"this target is unreachable" and "this target has exactly one solution".

Elimination
===========

By the Closure Theorem, the elimination ideal :math:`I \cap k[x_{s},\dots,x_n]`
generates the ideal of the Zariski closure of the projection of :math:`\V(I)`.
:cpp:func:`varietas::eliminated_generators` and ``ideal::eliminate(split)``
return the generators supported on the tail variables; the result generates the
elimination ideal whenever the order is an elimination order for that split,
which ``block_order`` is by construction.

Elimination is how the workspace is obtained (:doc:`workspace`) and how the
singular locus is pushed into it (:doc:`singularities`).

Saturation
==========

:cpp:func:`varietas::saturate` computes :math:`I : h^{\infty}` by the Rabinowitsch
trick — adjoin a variable :math:`u` and the relation :math:`1 - u h`, then
eliminate :math:`u`. This is what removes the spurious components the half-angle
substitution introduces at :math:`q_i = \pi`.

The same machinery gives an exact and exhaustive splitting for *any* :math:`h`:

.. math::

   \V(I) \;=\; \V(I : h^{\infty}) \;\cup\; \V(I + \langle h \rangle),

returned by :cpp:func:`varietas::split_along` as
:cpp:struct:`varietas::ideal_splitting`. It separates the branch where
:math:`h` does not vanish from the branch where it does; it is not a primary
decomposition, and turning it into one is what :doc:`../roadmap` calls
factorisation over :math:`\Q`.
