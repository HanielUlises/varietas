=================
varietas_codegen
=================

.. cpp:namespace:: varietas

The offline half of the library: the field the Gröbner computation actually runs
over, and the emitter that writes the result out. Namespace ``varietas`` for the
field, ``varietas::codegen`` for the rest.

Exact rationals
===============

``varietas/codegen/rational.hpp``

.. cpp:type:: rational

   Arbitrary-precision rationals backed by GMP, with the
   ``coefficient_traits`` specialisation the core algorithms consult.

.. cpp:function:: rational make_rational(std::int64_t numerator, std::int64_t denominator = 1)
.. cpp:function:: rational rational_from_string(const std::string& text, int base = 10)

The function field
==================

.. cpp:namespace:: varietas::codegen

``varietas/codegen/rational_function.hpp``

.. cpp:class:: template<std::size_t P> rational_function

   An element of :math:`\Q(p_1,\dots,p_P)`: a numerator and a denominator, each
   a polynomial in the :math:`P` pose parameters, normalised after every
   coefficient operation.

   Adjoining the pose to the **coefficient field** rather than to the polynomial
   ring is what lets one basis answer every pose instead of one basis per pose.
   It also specialises ``coefficient_traits``, which is why the same Buchberger
   implementation runs over it unmodified.

   .. rubric:: The cost, and the fast path

   Normalisation needs a polynomial gcd, and that is where the run goes:
   sampled over three minutes, about **86%** of the time is inside
   :cpp:func:`varietas::polynomial_gcd`, with the cost per call growing sharply
   as the parameter polynomials do.

   Most of those calls return 1, so normalisation asks a cheap question first.
   Specialising every parameter but one at fixed values in a small prime field
   leaves two univariate polynomials whose gcd is a Euclidean algorithm on
   machine integers; a constant answer there is strong evidence of coprimality,
   and the exact gcd is computed only when it is not. Every variable is kept in
   turn, since a factor involving only the specialised variables would collapse
   to a constant and be missed.

   The test is **evidence rather than proof**, and that is admissible precisely
   because it only ever decides whether to *skip* a cancellation: a skipped one
   leaves the fraction in higher terms, which costs size and never correctness.

   .. rubric:: What it buys

   Not much: about 29 ms down to about 25 ms on the reduced two-joint problem,
   and the three-parameter system still produces no answer. See
   :doc:`../roadmap`.

The solved system
=================

``varietas/codegen/parametric_solution.hpp``

.. cpp:struct:: template<std::size_t P> parametric_matrix

   An action matrix whose entries are rational functions of the pose.

.. cpp:struct:: template<std::size_t N, std::size_t P> parametric_solution

   Everything the emitter needs about a solved system, and nothing about how it
   was posed, so it can be built by hand in a unit test and the emitter is
   testable without a robot.

   :``order``: the :cpp:enum:`order_id` the basis was computed under.
   :``unknown_names``, ``parameter_names``: for comments and the generated
      signature; empty falls back to ``x0, x1, …``.
   :``quotient``: the standard monomial basis.
   :``action``: one :cpp:struct:`parametric_matrix` per unknown.
   :``one_index``: where the monomial :math:`1` sits. The eigenvalue method
      divides by the eigenvector's component there, and a functional that
      sends :math:`1` to zero is not an evaluation at a point.
   :``variable_coordinates``: row :math:`i` is the normal form of :math:`x_i`
      in the standard basis.

   .. cpp:function:: bool is_well_formed() const

      The invariants ``emit`` would otherwise have to trust. Checked at the top
      of ``emit``, because a header generated from an inconsistent solution
      compiles perfectly and answers wrongly.

.. cpp:function:: template<std::size_t N, std::size_t P, class Order> \
                  varietas::parametric_matrix<P> \
                  parametric_action_matrix(std::size_t variable, \
                                           const std::vector<polynomial<rational_function<P>, N, Order>>& basis, \
                                           const quotient_basis<N>& quotient)

.. cpp:function:: template<std::size_t N, std::size_t P, class Order> \
                  std::vector<std::vector<rational_function<P>>> \
                  parametric_variable_coordinates(const std::vector<polynomial<rational_function<P>, N, Order>>& basis, \
                                                  const quotient_basis<N>& quotient)

Emission
========

``varietas/codegen/emit.hpp``

.. cpp:enum:: runtime_kind

   ``matrices_only``: the header includes only ``<cstddef>`` and ``<cstdint>``
   and names nothing from this library. ``eigen``: the above plus ``solve``.

.. cpp:struct:: emit_options

   ``name``, ``name_space``, ``guard`` (defaulted from the two), ``source_note``
   (free text in the banner), ``runtime``, and ``epilogue``, the last being
   verbatim declarations placed inside the namespace after the struct, written
   out exactly as given.

.. cpp:function:: template<std::size_t N, std::size_t P> \
                  std::string emit(const parametric_solution<N, P>& solution, \
                                                      const emit_options& options = {})

   The header, as a string. See :doc:`../guide/generated_headers` for what it
   contains and how the denominator guard works.
