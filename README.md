# varietas
Algebraic kinematics for ROS 2. Exact inverse kinematics, workspace implicitization, and singularity variety decomposition via Gröbner bases and elimination theory.
Treats robot kinematics as what it is _a system of polynomial equations over a finite-dimensional quotient ring_ and solves it with the tools of computational algebraic geometry rather than with numerical iteration or structural pattern matching.

Given a URDF, the library rationalizes the forward kinematics map under the tangent half-angle substitution, computes a Gröbner basis of the resulting ideal offline, and emits header-only C++ that resolves the inverse kinematics at runtime by eigendecomposition of a fixed-size action matrix. Solutions are complete: the shape of the ideal certifies that no branch has been missed, and the same ideal yields the implicit equations of the reachable workspace and a decomposition of the singular locus into its irreducible components.

Where IKFast pattern-matches known kinematic structures and fails silently outside them, varietas either produces a solution set with a proof of completeness or reports precisely which hypothesis of the Extension Theorem was violated.

## State of the implementation

`varietas_core` is header-only and depends on Eigen alone. It is templated on the coefficient field so that the same code runs over `double` at runtime and over an exact rational type in the offline generator.

**Representation.** Exponent vectors with a cached total degree (`monomial`), sparse polynomials whose terms are kept sorted and reduced with the monomial order carried in the type (`polynomial`), so that two polynomials written under different orders cannot be combined.

**Monomial orders.** Lexicographic, graded lexicographic, graded reverse lexicographic, block orders `block_order<Split, First, Second>` for elimination, and weighted orders `weight_order<Weights, TieBreak>`. Each order records an `order_id`, which a generated header stores so that a mismatch between the order used offline and the one assumed at runtime is caught rather than silently producing a wrong solution set.

**Ideals.** The multivariate division algorithm with quotients and remainder; Buchberger's algorithm with the normal selection strategy and both Buchberger criteria, reporting how many critical pairs each criterion discarded; minimalisation and full reduction to the unique reduced Gröbner basis; ideal membership; detection of the unit ideal, that is, of an empty variety. The `ideal` class caches the basis and exposes `eliminate(split)`, whose output generates the elimination ideal whenever the order is an elimination order for that split.

**Quotient algebra.** Standard monomials by Macaulay's theorem, together with the finiteness verdict: the quotient is finite dimensional exactly when some leading monomial is a pure power of each variable. Multiplication operators on the quotient are assembled as dense action matrices in that basis.

**Solving.** Zero-dimensional systems are solved by the eigendecomposition of the action matrix of a generic linear form, in the formulation of Stetter and Möller: the left eigenvectors are the evaluation functionals at the points of the variety, and each coordinate is recovered by applying the functional to the normal form of the corresponding variable. Complex and real points are reported separately, a residual is available to certify any point, and every failure mode is named — empty variety, positive dimension, a defective eigenstructure — rather than returning a quietly truncated solution set.

Still to come: the URDF front end, the tangent half-angle rationalisation, exact rational arithmetic, code emission, workspace implicitization and the decomposition of the singular locus.
