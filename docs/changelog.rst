=========
Changelog
=========

0.2.0
=====

Three parameters. The parametric path now solves a tool position for three
joints, which is every arm the counting admits. See ``doc/parametric_cost.pdf``.

.. rubric:: Exactness

* ``modular_gcd``: Brown's dense modular gcd over :math:`\Q`, certified by
  exact division over :math:`\mathbb{Z}`, returning the cofactors that division
  computes. ``rational_function`` normalises with it, and the prime-field
  coprimality heuristic in front of the old gcd is removed: a constant image
  modulo one prime is a proof. The anthropomorphic arm over
  :math:`\Q(x,y,z)`, which produced nothing in fifteen minutes, solves in about
  twenty seconds.
* ``residue``: a word-sized prime field the whole library can compute over.
* ``reconstruction::reconstruct``: rational functions recovered from their
  values modulo primes.

.. rubric:: Inverse kinematics

* ``reconstructed_position_ik``: the parametric solve recovered from fixed-pose
  solves over prime fields, checked exactly at rational poses. Arms with
  offsets, on which the symbolic solve does not finish, take seconds.
* ``urdf_codegen --reconstruct``.

.. rubric:: Code generation

* A solution can carry the equations as posed; the generated ``solve()`` then
  takes up to two Newton steps from each point, keeping a step only when it
  lowers the residual. The worst residual of the decoupled demonstration
  solver goes from :math:`6.5\times10^{-8}` m to :math:`1.6\times10^{-15}` m.

0.1.0
=====

First tagged version. Everything in :doc:`status` marked *complete* is in it.

.. rubric:: Algebra

* Monomials with sixteen-bit exponents and a cached degree; sparse polynomials
  carrying their monomial order in the type.
* Lexicographic, graded lexicographic, graded reverse lexicographic, block and
  weighted orders, each with an ``order_id`` recorded in generated code.
* Multivariate division; Buchberger with the normal selection strategy and both
  criteria, reporting what each discarded; minimalisation and reduction to the
  unique reduced basis; membership; the unit ideal.
* Dimension from the leading terms, with the empty variety reported separately
  from dimension zero.
* Elimination under block orders; saturation by Rabinowitsch's trick, as a
  single elimination returning a basis under the caller's own order; exhaustive
  splitting along a chosen divisor.
* Standard monomials and the finiteness verdict; action matrices; the spectral
  solver of Stetter and Möller, with every failure mode named.
* Maximal minors by memoised Laplace expansion, the arrangement that never
  divides.

.. rubric:: Exactness

* ``varietas::rational`` over GMP, with the ``coefficient_traits``
  specialisation the core consults.
* ``rational_function<P>``: the pose adjoined to the coefficient field, with a
  prime-field coprimality test in front of the gcd.
* ``polynomial::prune`` rejected at compile time over an exact field.

.. rubric:: Kinematics

* Chains of revolute, prismatic and fixed joints over any coefficient field,
  with validation that names the joint and the defect.
* Rotations rational by construction, from quaternions or from a cosine–sine
  pair.
* Both rationalisations, half-angle and trigonometric, carried alongside each
  other, and required by test to agree.
* Workspace implicitization on the trigonometric ring, with the half-angle route
  kept for comparison.
* The singular locus: geometric Jacobian, maximal minors, dimension, and the
  image in the workspace, per task.

.. rubric:: Pipeline

* URDF front end recovering the chain exactly over :math:`\Q` by a projective
  quaternion search, with a per-joint audit.
* ``parametric_position_ik`` over :math:`\Q(\p)`; ``decoupled_position_ik``
  sweeping the base joint out.
* ``emit``: self-contained headers in two runtimes, with a denominator guard
  that recognises cancellation rather than comparing against zero.
* ``urdf_report``, ``urdf_codegen``, ``urdf_solve``, and two RViz
  demonstrations: one driving a model from the chain recovered from it, one
  drawing every configuration the generated solver returns for a moving
  target, against a solver emitted from the URDF during the build.

.. rubric:: Known limits

Two adjoined parameters; five joints for a fixed pose. See :doc:`status`.
