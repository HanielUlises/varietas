=========================
Workspace implicitization
=========================

Putting the pose coordinates **into the ring** instead of substituting them
turns the residuals into a description of the graph of the forward kinematics
map; eliminating the joint variables projects that graph onto the pose
coordinates, and by the Closure Theorem the result is the ideal of the Zariski
closure of the reachable set.

The Rabinowitsch variable and the joint variables are eliminated in a single
pass, since computing the saturation in full only to project it afterwards
produces an object nobody wants.

Closure is not reach
====================

The word *closure* carries weight the library should not paper over. The
reachable set of an arm is **semialgebraic**, an annulus, a shell, something
with a boundary, and boundaries are cut out by inequalities that no ideal
expresses. Where the workspace is full dimensional the closure is everything,
and the elimination ideal returns only what held identically.

.. figure:: ../figures/torus_workspace.svg
   :width: 75%
   :alt: The Zariski closure of the workspace of a two-joint arm with
         perpendicular axes

   The workspace closure of the arm with perpendicular axes: the zero set of the
   single quartic :math:`(x^2+y^2+z^2+3)^2 = 16(x^2+y^2)` that elimination
   returned. No torus is parameterised in drawing this; the surface is traced
   from the polynomial itself.

Three cases make the distinction concrete:

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Arm
     - What elimination returns
   * - One joint
     - The circle :math:`z=0,\ x^2+y^2=1`, exactly.
   * - Planar two-joint
     - :math:`z=0` and nothing else, so a point a hundred units away satisfies
       every equation the closure has. The reachable annulus is not an ideal.
   * - Two joints, perpendicular axes
     - The quartic cutting out the torus. The map is not dominant, and the
       closure is a genuine surface.

:cpp:func:`varietas::workspace_is_dense` reports the middle case, rather than
leaving an empty result to be misread as a failure.

Which formulation, and what it cost
===================================

The exercise also measured something about the formulation. Eliminating for the
torus in the half-angle ring takes **seventy seconds**. Presenting the same arm
with two variables per joint and the relation :math:`c^2+s^2=1`, which has no
denominators and therefore needs no saturation at all, eliminates in **sixty
milliseconds** and returns the identical quartic. A factor of a thousand.

The half-angle substitution buys one variable per joint rather than two, which
is what matters for the zero-dimensional inverse kinematics the library is built
around, and pays for it with a denominator, which is what matters here.

.. list-table::
   :header-rows: 1
   :widths: 28 36 36

   * -
     - Half-angle
     - Trigonometric
   * - Variables per joint
     - 1
     - 2
   * - Denominators
     - yes, saturation required
     - none
   * - :math:`q=\pi`
     - sent to infinity
     - ordinary point
   * - Used for
     - inverse kinematics
     - implicitization, singularities

:cpp:func:`varietas::workspace_relations` therefore runs on the trigonometric
ring, and :cpp:func:`varietas::workspace_relations_half_angle` is kept so that
the two can be **compared rather than trusted**: both are run on the torus arm
and their reduced bases are required to be equal, which is what licenses the
migration. The elimination ideal does not depend on the parameterisation it was
computed through; only the cost does.

The trigonometric transform
===========================

Writing :math:`c_i=\cos q_i` and :math:`s_i=\sin q_i` as independent variables
makes Rodrigues' formula polynomial as it stands,

.. math::

   R = I + s\,[u]_\times + (1-c)\,[u]_\times^{2},

so a chain composes as an ordinary product of polynomial matrices and
:cpp:class:`varietas::trigonometric_transform` has no denominator field to
carry.

The trigonometric identity enters as a generator :math:`c_i^2+s_i^2-1` on the
same footing as the pose residuals: a genuine equation of the problem rather
than damage repair. Nothing is cleared, so nothing spurious is attached, so
there is nothing to saturate away: the loci :math:`t_i=\pm i` have no
counterpart here, and :math:`\V(c^2+s^2-1)` is exactly the parameter space of a
revolute joint, including the configuration :math:`q=\pi` that the half-angle
map sends to infinity.
