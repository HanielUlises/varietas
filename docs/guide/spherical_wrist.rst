=============================
Six joints, a spherical wrist
=============================

A tool pose is six parameters, far more than a parametric solve can carry. Most
industrial arms never need it to, because they are built so that their last
three axes meet in a point, the **wrist centre**. Turning those three joints
rotates the tool about that point without moving it, so where the wrist centre
sits depends on the first three joints alone.

That splits the problem in two, which is Pieper's decomposition. The target
pose fixes the wrist centre. Placing it is a three-joint position problem,
which is what the rest of varietas solves, and the rotation the wrist must
then supply is split into its three joint angles in closed form.

.. code-block:: sh

   ros2 run varietas_urdf urdf_codegen arm.urdf arm_ik.hpp --wrist

On the industrial arm in ``varietas_urdf/test/data/industrial_6r.urdf`` this
takes under two seconds and writes a header that solves the full pose.

What is exact and what is closed form
=====================================

Whether the three wrist axes meet is decided over :math:`\Q`, on the chain
recovered exactly from the URDF, and so is the wrist centre. Decimals such as
0.025 and 0.42 arrive as 1/40 and 21/50, so the axes meet exactly where the
drawing says they do rather than missing one another by rounding.

The arm that places the wrist centre goes to the three-joint solvers unchanged.
It is decoupled when its base sweeps a plane (:doc:`decoupling`), which is the
case for the usual industrial geometry, and reconstructed over
:math:`\Q(x,y,z)` when it does not (:cpp:func:`varietas::ik::reconstructed_position_ik`).
Either way every configuration of the arm is returned, and the count is the
quotient dimension.

The rotation left for the wrist is :math:`R_{03}^{T} R\, R_{\mathrm{tool}}^{T}`,
and it is split by the subproblems of Paden and Kahan. The sixth axis is fixed
by its own turn, which leaves two rotations carrying one known vector onto
another, the second subproblem, with two solutions; the sixth angle then
follows from the first subproblem. The wrist axes need only meet. They need not
be orthogonal.

How many configurations
=======================

Up to four for the arm and two for the wrist, so up to eight, which is the most
a six-joint arm with a spherical wrist has. On the industrial arm, twenty
thousand random poses return eight configurations about nine times in ten and
four otherwise, and the configuration each pose came from is always among them.

A wrist whose consecutive axes are not perpendicular reaches only some
orientations. For such a wrist an arm configuration can put the centre in place
and still leave a rotation no setting of the last three joints makes, and that
configuration has no completion, so the count is two, four, six or eight. The
solver discards those configurations rather than returning a wrist that misses
the orientation.

The wrist's singularity
=======================

The commonest wrist rolls, pitches and rolls again, so its first and last axes
are collinear whenever the pitch is zero. Only their sum then matters, and any
split between them is a solution. The decomposition accepts such a wrist, since
its axes still meet, and the solver returns one split at the singularity rather
than refusing the pose. At the singularity itself the configurations it
returns reproduce the pose to rounding. Within about :math:`10^{-9}` rad of it
they lose digits, to about :math:`10^{-9}` in the pose, because the angle is
read off a vector that is nearly parallel to the axis it turns about.

What it refuses
===============

``urdf_codegen --wrist`` refuses an arm that is not six revolute joints, an arm
whose consecutive wrist axes are parallel, and an arm whose last three axes do
not meet, saying which. The KUKA iiwa has seven joints and is refused on the
count. An arm whose wrist axes miss by a centimetre is refused for not having a
centre, which is the right answer, because the position of the tool then
depends on all six joints and the decomposition does not hold.
