=========
Changelog
=========

0.2.0
=====

Three parameters. The parametric path now solves a tool position for three
joints, which is every arm the counting admits. See ``doc/parametric_cost.pdf``.

.. rubric:: Exactness

* ``modular_gcd``, Brown's dense modular gcd over :math:`\Q`, certified by
  exact division over :math:`\mathbb{Z}`, returning the cofactors that division
  computes. ``rational_function`` normalises with it, and the prime-field
  coprimality heuristic in front of the old gcd is removed, since a constant image
  modulo one prime is a proof. The anthropomorphic arm over
  :math:`\Q(x,y,z)`, which produced nothing in fifteen minutes, solves in about
  twenty seconds.
* ``residue``, a word-sized prime field the whole library can compute over.
* ``reconstruction::reconstruct``, which recovers rational functions from their
  values modulo primes.

.. rubric:: Factorisation

* ``squarefree_decomposition``, Yun's algorithm over :math:`\Q` in any number
  of variables, on the modular gcd and its cofactors.
* ``factor_univariate``, Zassenhaus's algorithm: Cantor–Zassenhaus modulo the
  best of five primes, quadratic Hensel lifting down a factor tree, and a
  recombination whose acceptance test is exact rather than a trial. It agrees
  with SymPy on three hundred random products, and finds the polynomials of
  Swinnerton-Dyer irreducible.
* ``factor``, in any number of variables. The univariate image at a point is
  lifted by Hensel's lemma with the leading coefficient imposed on both halves
  of each division into two groups, and every factor is accepted only when it
  divides, exactly over :math:`\mathbb{Z}`. A point that splits a true factor
  falls back to power series and recombination. It agrees with SymPy on six
  hundred random products in two and three variables, and factors products of
  about 1800 terms in six variables in a second or two.

.. rubric:: Inverse kinematics

* ``reconstructed_position_ik``, the parametric solve recovered from fixed-pose
  solves over prime fields, checked exactly at rational poses. Arms with
  offsets, on which the symbolic solve does not finish, take seconds.
* ``urdf_codegen --reconstruct``.

.. rubric:: Code generation

* A solution can carry the equations as posed; the generated ``solve()`` then
  takes up to two Newton steps from each point, keeping a step only when it
  lowers the residual. The worst residual of the decoupled demonstration
  solver goes from :math:`6.5\times10^{-8}` m to :math:`1.6\times10^{-15}` m.

.. rubric:: Six joints

* ``decompose_spherical_wrist`` and ``emit_spherical_wrist``. A six-joint arm
  whose last three axes meet is split at its wrist centre, decided exactly
  over :math:`\Q`. The arm that places the centre is decoupled or
  reconstructed, and the wrist is split by the subproblems of Paden and
  Kahan, so the generated header solves the full pose and returns up to eight
  configurations. ``urdf_codegen --wrist``.

.. rubric:: Joint ranges and MoveIt

* The decoupled and six-joint headers carry the joints' ranges and a
  ``solve_within_limits`` that keeps the configurations every joint can reach,
  each angle moved by whole turns into its range.
* ``runtime.hpp`` evaluates a solution in process, with the same arithmetic
  the generated headers carry, for an arm known only at run time.
* ``varietas_moveit``, a MoveIt kinematics plugin. For a six-joint arm whose
  last three axes meet it returns every configuration reaching the pose, inside
  the ranges MoveIt knows and nearest the seed first. On the industrial test
  arm it solves every reachable pose in a median of 11.5 µs, where the default
  KDL solver solves 97.8% in a median of 1.5 ms.

.. rubric:: Continuous integration

* Every push builds all six packages and runs every test in a ROS 2 Humble
  container, with the assertions left on in an optimised build.
* Two dependencies that rosdep could not resolve, ``libgmp-dev`` and
  ``orocos_kdl``, are declared by their rosdep keys.

.. rubric:: Demonstrations

* ``solve_all.launch.py``, a third demonstration. An arm with three pairwise
  skew axes, solved in full by reconstruction during the build, drawn with
  every configuration beside a warm-started damped least squares iteration.
  The node can write a trace of every tick, from which the published
  recording is captioned.

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
