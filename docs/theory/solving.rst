==============================================
Solving: the quotient algebra and its spectrum
==============================================

Finiteness
==========

Let :math:`\prec` be a monomial order, :math:`G=\mathrm{GB}_{\prec}(I)` the
reduced Gröbner basis, and :math:`A=k[\boldsymbol{t}]/I` the quotient algebra.

By the Finiteness Theorem, :math:`\dim_k A<\infty` if and only if :math:`\V(I)`
is finite, a condition decidable from the leading terms of :math:`G` alone:
the quotient is finite dimensional exactly when some leading monomial is a pure
power of each variable. In that case :math:`\dim_k A` **bounds the number of
solutions counted with multiplicity**, so the basis itself certifies that no
branch of the inverse kinematics has been omitted.

This is the completeness claim that separates varietas from a numerical solver.
A Newton iteration returns a solution; the quotient dimension says how many
there are.

Standard monomials
==================

:cpp:func:`varietas::standard_monomials` computes the monomials not in the
initial ideal, a basis of :math:`A` by Macaulay's theorem, together with the
finiteness verdict, as :cpp:struct:`varietas::quotient_basis`.

Action matrices
===============

Multiplication by a fixed element is a linear operator on :math:`A`, and it is
assembled as a dense matrix in the standard monomial basis by
:cpp:func:`varietas::action_matrix` (or, for a single variable,
:cpp:func:`varietas::variable_action_matrix`).

The eigenvalue method
=====================

For a generic linear form :math:`\ell \in A`, the points of :math:`\V(I)` are
recovered from the spectrum of :math:`M_\ell : A \to A`. By the theorem of
Stetter and Möller the **left** eigenvectors of :math:`M_\ell` are the
evaluation functionals at those points: such a functional sends each standard
monomial :math:`m` to :math:`m(p)`. Left eigenvectors of :math:`M` are right
eigenvectors of :math:`M^{\mathsf T}`, which is what is actually decomposed.

Each coordinate is then obtained by applying the functional to the normal form
of :math:`t_i`, divided by what the functional gives to :math:`1`. That division
is what makes it an evaluation rather than a multiple of one.

.. note::

   A single variable will not do as the operator. If two distinct solutions give
   it the same value, the eigenvalue is repeated, the eigenvector is an arbitrary
   element of a plane, and it carries no point. varietas uses a fixed linear
   combination of the action matrices, which separates the solutions except on a
   proper algebraic subset of the coefficient space, and being fixed rather than
   random, two runs on the same pose return the same points in the same order.

Failure is named
================

:cpp:func:`varietas::solve_zero_dimensional` returns a
:cpp:struct:`varietas::solution_set` whose ``status`` is one of

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - ``solve_status``
     - Meaning
   * - ``ok``
     - The point set is complete.
   * - ``empty_variety``
     - The basis is :math:`\{1\}`: no configuration reaches this target.
   * - ``positive_dimensional``
     - Infinitely many solutions; no action matrix exists. Usually a
       redundant arm, or two joints turning about one axis.
   * - ``eigensolver_failed``
     - Eigen did not converge.
   * - ``deficient_eigenstructure``
     - Some eigenvector carried no evaluation functional. The points recovered
       are still returned, but the set is **not** certified complete and the
       caller must treat it as a failure of genericity.

Complex and real points are reported separately, with
``solution_set::real_points(tolerance)`` filtering to the ones a manipulator can
actually be commanded to, and :cpp:func:`varietas::residual` certifies any
point against the original generators.
