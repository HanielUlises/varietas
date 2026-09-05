==============
varietas_core
==============

.. cpp:namespace:: varietas

Header-only, Eigen only. Templated on the coefficient type throughout.

Configuration
=============

``varietas/core/config.hpp``

.. cpp:type:: std::size_t index_type
.. cpp:type:: double scalar

.. cpp:struct:: template<class Coeff> coefficient_traits

   The field interface every algorithm consults: ``zero``, ``one``, ``is_zero``,
   ``inverse``, ``negate``, ``to_double``, ``from_double``, and the constant
   ``is_exact``. Specialised for ``double``, for
   :cpp:type:`varietas::rational`, and for
   :cpp:class:`varietas::codegen::rational_function`.

   ``is_exact`` is what lets the library refuse, at compile time, an operation
   that is only meaningful over a floating point field —
   ``polynomial::prune`` being the example.

Monomials
=========

``varietas/core/monomial.hpp``

.. cpp:class:: template<std::size_t N> monomial

   An exponent vector with a cached total degree.

   .. cpp:type:: std::uint16_t exponent_type

      Sixteen bits, not eight. Parameter polynomials over a function field
      reach the old ceiling — degree 254 in one variable was observed — after
      which a product wrapped to a *different monomial* and a polynomial that
      divided another silently stopped dividing it. The product now asserts
      rather than wraps.

   .. cpp:function:: static constexpr monomial one() noexcept
   .. cpp:function:: static constexpr monomial variable(std::size_t i, exponent_type power = 1) noexcept
   .. cpp:function:: constexpr degree_type degree() const noexcept
   .. cpp:function:: static constexpr bool divides(const monomial& a, const monomial& b) noexcept
   .. cpp:function:: static constexpr monomial divide(const monomial& b, const monomial& a) noexcept
   .. cpp:function:: static constexpr monomial lcm(const monomial& a, const monomial& b) noexcept
   .. cpp:function:: static constexpr monomial gcd(const monomial& a, const monomial& b) noexcept

Polynomials
===========

``varietas/core/polynomial.hpp``

.. cpp:class:: template<class Coeff, std::size_t N, class Order> polynomial

   A sparse polynomial whose terms are kept sorted and reduced under the
   monomial order **carried in the type**, so that two polynomials written under
   different orders cannot be combined.

   .. cpp:function:: static polynomial constant(const Coeff& c)
   .. cpp:function:: static polynomial variable(std::size_t i, exponent_type power = 1)
   .. cpp:function:: static polynomial from_monomial(const monomial_type& m, const Coeff& c)
   .. cpp:function:: const term& leading_term() const
   .. cpp:function:: const monomial_type& leading_monomial() const
   .. cpp:function:: const Coeff& leading_coefficient() const
   .. cpp:function:: Coeff coefficient_of(const monomial_type& m) const
   .. cpp:function:: polynomial monic() const
   .. cpp:function:: Coeff evaluate(const std::array<Coeff, N>& point) const
   .. cpp:function:: void prune(const Coeff& tolerance)

      Discards negligible terms. ``static_assert``\ s away for an exact
      coefficient type: over :math:`\Q` a small term is not noise, it is the
      ideal.

Monomial orders
===============

``varietas/core/order/``

.. cpp:class:: lex

   Lexicographic.

.. cpp:class:: grlex

   Graded lexicographic.

.. cpp:class:: grevlex

   Graded reverse lexicographic. The default for solving.

.. cpp:class:: template<std::size_t Split, class First, class Second> block_order

   An elimination order for the split: everything supported on the first block
   is greater than anything that is not. This is what makes elimination and
   saturation work.

.. cpp:class:: template<class Weights, class TieBreak> weight_order

.. cpp:enum:: order_id

   ``lex``, ``grlex``, ``grevlex``, ``block``, ``weight``. Recorded in every
   generated header so that a mismatch between the order used offline and the
   one assumed at runtime is caught rather than producing a wrong solution set.
   ``to_string(order_id)`` names it.

Division
========

``varietas/core/ideal/division.hpp``

.. cpp:struct:: template<class Poly> division_result

   Quotients and remainder.

.. cpp:function:: template<class Coeff, std::size_t N, class Order> \
                  division_result<polynomial<Coeff, N, Order>> \
                  divide(const polynomial<Coeff, N, Order>& f, \
                         const std::vector<polynomial<Coeff, N, Order>>& divisors)

.. cpp:function:: template<class Coeff, std::size_t N, class Order> \
                  polynomial<Coeff, N, Order> \
                  normal_form(const polynomial<Coeff, N, Order>& f, \
                              const std::vector<polynomial<Coeff, N, Order>>& divisors)

.. cpp:function:: template<class Coeff, std::size_t N, class Order> \
                  polynomial<Coeff, N, Order> \
                  normal_form_excluding(const polynomial<Coeff, N, Order>& f, \
                                        const std::vector<polynomial<Coeff, N, Order>>& divisors, \
                                        std::size_t excluded)

   Used by the minimalisation pass, where a basis element must be reduced
   against every other one but not against itself.

Gröbner bases
=============

``varietas/core/ideal/buchberger.hpp``

.. cpp:struct:: buchberger_statistics

   How many critical pairs each Buchberger criterion discarded, and how many
   S-polynomials were actually reduced. Reported rather than hidden: the two
   counts are the difference between a run that finishes and one that does not.

.. cpp:function:: template<class Coeff, std::size_t N, class Order> \
                  polynomial<Coeff, N, Order> \
                  s_polynomial(const polynomial<Coeff, N, Order>& f, \
                               const polynomial<Coeff, N, Order>& g)

.. cpp:function:: template<class Coeff, std::size_t N, class Order> \
                  std::vector<polynomial<Coeff, N, Order>> \
                  buchberger(std::vector<polynomial<Coeff, N, Order>> generators, \
                             buchberger_statistics* stats = nullptr)

   Normal selection strategy, both Buchberger criteria.

.. cpp:function:: template<class Coeff, std::size_t N, class Order> \
                  void reduce_groebner_basis(std::vector<polynomial<Coeff, N, Order>>& basis)

.. cpp:function:: template<class Coeff, std::size_t N, class Order> \
                  std::vector<polynomial<Coeff, N, Order>> \
                  groebner_basis(std::vector<polynomial<Coeff, N, Order>> generators)

   Buchberger, minimalised and fully reduced: the **unique** reduced Gröbner
   basis for the order.

.. cpp:function:: template<class Coeff, std::size_t N, class Order> \
                  bool is_member(const polynomial<Coeff, N, Order>& f, \
                                 const std::vector<polynomial<Coeff, N, Order>>& basis)

.. cpp:function:: template<class Coeff, std::size_t N, class Order> \
                  bool is_unit_ideal(const std::vector<polynomial<Coeff, N, Order>>& basis)

   An empty variety, in the language of the weak Nullstellensatz.

Dimension
=========

``varietas/core/ideal/dimension.hpp``

.. cpp:struct:: template<std::size_t N> affine_dimension

   The dimension, the maximal independent set that witnesses it, and a flag for
   the empty variety — which is *not* dimension zero, that being a point.

.. cpp:function:: template<class Coeff, std::size_t N, class Order> \
                  affine_dimension<N> \
                  ideal_dimension(const std::vector<polynomial<Coeff, N, Order>>& basis)

Elimination
===========

``varietas/core/ideal/ideal.hpp``

.. cpp:function:: template<class Coeff, std::size_t N, class Order> \
                  std::vector<polynomial<Coeff, N, Order>> \
                  eliminated_generators(const std::vector<polynomial<Coeff, N, Order>>& basis, \
                                        std::size_t split)

   The basis elements supported on the tail variables. These generate the
   elimination ideal whenever the order is an elimination order for the split.

.. cpp:class:: template<class Coeff, std::size_t N, class Order> ideal

   Caches the reduced basis; exposes ``basis()``, ``eliminate(split)``,
   ``dimension()``, ``contains(f)``, ``is_unit()``.

Saturation
==========

``varietas/core/ideal/saturation.hpp``

.. cpp:function:: template<class Coeff, std::size_t N, class Order> \
                  std::vector<polynomial<Coeff, N, Order>> \
                  saturate(const std::vector<polynomial<Coeff, N, Order>>& generators, \
                           const polynomial<Coeff, N, Order>& h)

   :math:`I : h^{\infty}` by Rabinowitsch's trick, computed as a single
   elimination under a block order whose trailing block is the caller's own
   order — so the result is already a Gröbner basis under that order.

.. cpp:function:: template<class Coeff, std::size_t N, class Order> \
                  std::vector<polynomial<Coeff, N, Order>> \
                  saturate_by_product(const std::vector<polynomial<Coeff, N, Order>>& generators, \
                                      const std::vector<polynomial<Coeff, N, Order>>& factors)

.. cpp:struct:: template<class Coeff, std::size_t N, class Order> ideal_splitting

   The two branches of :math:`\V(I) = \V(I:h^\infty) \cup \V(I + \langle h
   \rangle)`.

.. cpp:function:: template<class Coeff, std::size_t N, class Order> \
                  ideal_splitting<Coeff, N, Order> \
                  split_along(const std::vector<polynomial<Coeff, N, Order>>& generators, \
                              const polynomial<Coeff, N, Order>& h)

   Exact and exhaustive — but the caller chooses :math:`h`. Not a primary
   decomposition; see :doc:`../roadmap`.

.. cpp:function:: template<class Coeff, std::size_t N, class Order> \
                  polynomial<Coeff, N, Order> half_angle_denominator(std::size_t i)

   :math:`1 + t_i^2`.

Quotient algebra
================

``varietas/core/quotient/``

.. cpp:struct:: template<std::size_t N> quotient_basis

   The standard monomials, by Macaulay's theorem, together with
   ``is_zero_dimensional`` — true exactly when some leading monomial is a pure
   power of each variable — and ``dimension()``, which is :math:`\dim_k A`.

.. cpp:function:: template<class Coeff, std::size_t N, class Order> \
                  quotient_basis<N> \
                  standard_monomials(const std::vector<polynomial<Coeff, N, Order>>& basis)

.. cpp:function:: template<class Coeff, std::size_t N, class Order> \
                  action_matrix_type \
                  action_matrix(const polynomial<Coeff, N, Order>& multiplier, \
                                const std::vector<polynomial<Coeff, N, Order>>& basis, \
                                const quotient_basis<N>& quotient)

.. cpp:function:: template<class Coeff, std::size_t N, class Order> \
                  action_matrix_type \
                  variable_action_matrix(std::size_t variable, \
                                         const std::vector<polynomial<Coeff, N, Order>>& basis, \
                                         const quotient_basis<N>& quotient)

Solving
=======

``varietas/core/solve/spectral.hpp``

.. cpp:enum:: solve_status

   ``ok``, ``empty_variety``, ``positive_dimensional``, ``eigensolver_failed``,
   ``deficient_eigenstructure``. See :doc:`../theory/solving` for what each
   means and why none of them is a silently truncated answer.

.. cpp:struct:: template<std::size_t N> solution_set

   ``status``, ``points`` (over :math:`\C`, counted once each),
   ``quotient_dimension``, ``discarded_eigenvectors``; ``ok()`` and
   ``real_points(tolerance = 1e-8)``.

.. cpp:function:: template<class Coeff, std::size_t N, class Order> \
                  solution_set<N> \
                  solve_zero_dimensional(const std::vector<polynomial<Coeff, N, Order>>& basis)

.. cpp:function:: template<class Coeff, std::size_t N, class Order> \
                  double residual(const std::vector<polynomial<Coeff, N, Order>>& generators, \
                                  const std::array<std::complex<double>, N>& point)

   Certifies a point against the original generators, which is a different
   question from certifying it against the basis.

Polynomial gcd
==============

``varietas/core/gcd.hpp``

The subresultant remainder sequence, plus the pieces it needs:
``degree_in``, ``coefficient_in``,
``leading_coefficient_in``, ``divide_exact``,
``pseudo_remainder_in``, ``content_in``,
``main_variable``.

.. cpp:function:: template<class Coeff, std::size_t N, class Order> \
                  polynomial<Coeff, N, Order> \
                  polynomial_gcd(const polynomial<Coeff, N, Order>& a, \
                                 const polynomial<Coeff, N, Order>& b)

   This is where the parametric solve spends about **86% of its time**. See
   :doc:`../status` for the measurement and :doc:`../roadmap` for what would
   change it.

Minors
======

``varietas/core/minors.hpp``

.. cpp:class:: template<class Element> dense_matrix

.. cpp:function:: template<class Element> Element determinant(const dense_matrix<Element>& m)

.. cpp:function:: template<class Element> std::vector<Element> maximal_minors(const dense_matrix<Element>& m)

   Laplace expansion memoised on the set of remaining columns: :math:`2^k k`
   multiplications rather than :math:`k!`, and — the point — an arrangement
   that **never divides**. Gaussian elimination needs to invert a pivot and
   Bareiss needs exact division; polynomials offer neither.

Embedding
=========

``varietas/core/embed.hpp``

``embed_at``, ``project_from``,
``embed_polynomial``, ``project_polynomial``. Variable shifting
between rings of different sizes, re-sorting terms under the target order rather
than copying a term list sorted under a different one. Saturation is the caller.
