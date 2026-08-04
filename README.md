# varietas
Algebraic kinematics for ROS 2. Exact inverse kinematics, workspace implicitization, and singularity variety decomposition via Gröbner bases and elimination theory.
Treats robot kinematics as what it is _a system of polynomial equations over a finite-dimensional quotient ring_ and solves it with the tools of computational algebraic geometry rather than with numerical iteration or structural pattern matching.

Given a URDF, the library rationalizes the forward kinematics map under the tangent half-angle substitution, computes a Gröbner basis of the resulting ideal offline, and emits header-only C++ that resolves the inverse kinematics at runtime by eigendecomposition of a fixed-size action matrix. Solutions are complete: the shape of the ideal certifies that no branch has been missed, and the same ideal yields the implicit equations of the reachable workspace and a decomposition of the singular locus into its irreducible components.

Where IKFast pattern-matches known kinematic structures and fails silently outside them, varietas either produces a solution set with a proof of completeness or reports precisely which hypothesis of the Extension Theorem was violated.
