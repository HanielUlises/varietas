=================================
The problem, stated algebraically
=================================

.. figure:: ../figures/kinematics.png
   :width: 90%
   :alt: Algebraic formulation of the kinematics of a planar serial chain

   Kinematics of a planar serial chain with revolute joints, and its
   formulation as a zero-dimensional polynomial system. Forward kinematics
   evaluates the map; inverse kinematics computes the fibre, which varietas
   realises as the variety of an ideal in a quotient ring of finite dimension.

The map
=======

Let the manipulator have :math:`n` revolute joints with configuration
:math:`\boldsymbol{q}=(q_1,\dots,q_n)\in\mathbb{T}^{n}` and link parameters
:math:`L_1,\dots,L_n`. The forward kinematics is the analytic map

.. math::

   \mathcal{F}:\mathbb{T}^{n}\longrightarrow SE(m),
   \qquad \boldsymbol{q}\longmapsto \p,

with :math:`\p=(x,y,\theta)` for :math:`m=2` and :math:`\p\in SE(3)` in the
spatial case. Inverse kinematics is the determination of the fibre
:math:`\mathcal{F}^{-1}(\p)` for a prescribed pose :math:`\p` — a problem that
is transcendental as stated and therefore not directly amenable to exact
methods.

Rationalisation
===============

Under the tangent half-angle substitution :math:`t_i=\tan(q_i/2)` one has

.. math::

   \cos q_i=\frac{1-t_i^{2}}{1+t_i^{2}},
   \qquad
   \sin q_i=\frac{2t_i}{1+t_i^{2}},

so that, after clearing the denominators :math:`1+t_i^{2}`, the pose equations
become polynomial. Writing :math:`f_x, f_y, f_\theta \in k[t_1,\dots,t_n]` for
the resulting numerators of the pose residuals, the inverse kinematics problem
at :math:`\p` is the affine variety of the ideal

.. math::

   I_{\p}=\big\langle f_x,\;f_y,\;f_\theta\big\rangle
   \;\subseteq\;k[t_1,\dots,t_n],
   \qquad
   \mathcal{F}^{-1}(\p)\;\cong\;\V\!\left(I_{\p}\right),

the correspondence being bijective away from the locus :math:`q_i=\pi`, which
the substitution sends to infinity and which is removed by saturating with
respect to :math:`\prod_i\left(1+t_i^{2}\right)`.

In the library this is
:cpp:func:`varietas::rational_forward_kinematics` and the generator functions
beside it; the saturation is :cpp:func:`varietas::saturate_by_product` with
:cpp:func:`varietas::half_angle_denominator`.

The trigonometric alternative
=============================

The half-angle substitution is not the only rationalisation available. Adjoining
:math:`c_i=\cos q_i` and :math:`s_i=\sin q_i` with the circle relations
:math:`c_i^2+s_i^2=1` gives a polynomial system in :math:`2n` variables with no
excluded locus at all — nothing is sent to infinity, so nothing has to be
saturated back.

The trade is dimension against cleanliness: twice the variables, but every
configuration is represented, including the :math:`q_i=\pi` that the half-angle
chart loses. varietas keeps both.
:cpp:func:`varietas::trigonometric_forward_kinematics` and
:cpp:func:`varietas::circle_relations` build this form, and it is what the
singularity and workspace constructions use, where losing a configuration would
mean losing part of the answer.

Why exactness is the whole point
================================

Each of the statements above — that the quotient is finite dimensional, that a
subset of variables is independent, that a leading monomial is a pure power —
is a statement about *whether a coefficient is zero*. Over floating point those
are not decidable questions: a coefficient that vanishes mathematically almost
never evaluates to ``0.0``, and a Gröbner basis computed from wrong zero tests
is not a Gröbner basis of anything in particular.

So the whole offline path runs over :math:`\Q` with GMP integers underneath,
and the only numerical operation in the library is the eigendecomposition that
reads the configurations off the action matrix at the very end. See
:doc:`exactness`.
