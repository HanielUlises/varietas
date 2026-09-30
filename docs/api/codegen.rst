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

   .. rubric:: The cost

   Normalisation needs a polynomial gcd after every operation, and with the
   subresultant :cpp:func:`varietas::polynomial_gcd` that was where the run
   went: about **86%** of the time, at a cost growing as the fourth power of the
   operand size. Normalisation now calls
   :cpp:func:`varietas::modular_gcd_with_cofactors`, and every fraction is held
   in lowest terms, not merely usually.

Modular gcd
===========

``varietas/codegen/modular_gcd.hpp``

.. cpp:function:: template<std::size_t N, class Order> \
                  gcd_and_cofactors<polynomial<rational, N, Order>> \
                  modular_gcd_with_cofactors(const polynomial<rational, N, Order>& a, \
                                             const polynomial<rational, N, Order>& b)

   The gcd of two polynomials over :math:`\Q`, monic under ``Order``, together
   with ``a / gcd`` and ``b / gcd``. Brown's dense algorithm: the primitive
   integer operands are reduced modulo word-sized primes, specialised one
   variable at a time down to a Euclidean algorithm on machine words, and
   rebuilt by Newton interpolation; images modulo several primes are joined by
   the Chinese remainder theorem and read back by rational reconstruction.

   An image gcd can be too large but never too small, and too large shows in its
   leading monomial, which is how unlucky primes and evaluation points are
   recognised and discarded. The answer is **certified** by exact division of
   both operands over :math:`\mathbb{Z}`, and the quotients of that division are
   the cofactors returned. A constant image modulo a single prime proves the
   operands coprime, which is how most calls from normalisation end.

   On dense trivariate operands the cost grows as about the 1.5th power of the
   number of terms, against 3.8 for the subresultant sequence; at 455 terms one
   gcd takes under 5 ms rather than 91 s.

.. cpp:function:: template<std::size_t N, class Order> \
                  polynomial<rational, N, Order> \
                  modular_gcd(const polynomial<rational, N, Order>& a, \
                              const polynomial<rational, N, Order>& b)

   The gcd alone.

Reconstruction
==============

``varietas/codegen/prime_field.hpp``, ``varietas/codegen/reconstruct.hpp``

.. cpp:class:: residue

   :math:`\mathbb{Z}/p` for a word-sized prime selected per thread, with the
   ``coefficient_traits`` specialisation that lets Buchberger and the kinematics
   run over it. A rational whose denominator the prime divides has no image,
   and constructing one sets a flag the caller checks.

.. cpp:function:: template<std::size_t P> \
                  bool reconstruction::reconstruct(const black_box<P>& box, \
                                                   const prime_hook& prepare, \
                                                   std::size_t count, \
                                                   std::vector<rational_function<P>>& out, \
                                                   statistics* stats = nullptr, \
                                                   const options& opts = options())

   Recovers ``count`` rational functions of ``P`` parameters from a black box
   that evaluates them at points of :math:`\mathbb{Z}/p`. Per prime: degrees from
   a random line by univariate rational reconstruction under the maximal
   quotient rule, then coefficients from one linear system per function,
   normalised so the leading coefficient of the denominator is one; after the
   first prime only the monomials it found are solved for. Primes are joined as
   in the modular gcd, until two consecutive ones change nothing.

   The result is right with high probability rather than certainly, since there
   is no final exact test on a black box. :doc:`ik` checks it exactly at
   rational poses before using it.

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
