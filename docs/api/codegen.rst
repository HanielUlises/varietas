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
   went, about **86%** of the time, at a cost growing as the fourth power of the
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
   with ``a / gcd`` and ``b / gcd``, by Brown's dense algorithm. The primitive
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

Factorisation
=============

``varietas/codegen/factor.hpp``

.. cpp:struct:: template<class Poly> factorisation

   :math:`f = u\,g_1^{e_1}\cdots g_k^{e_k}`: the ``unit`` :math:`u`, the
   leading coefficient of :math:`f`, and the ``factors`` as
   ``factor_power{base, multiplicity}``, each base monic under the order.

.. cpp:function:: template<std::size_t N, class Order> \
                  factorisation<polynomial<rational, N, Order>> \
                  squarefree_decomposition(const polynomial<rational, N, Order>& f)

   :math:`f = u\,g_1 g_2^2\cdots g_k^k` with the :math:`g_i` squarefree and
   pairwise coprime, one for each multiplicity that occurs, in **any number of
   variables**. The content with respect to the first variable is decomposed
   recursively, and the primitive part by Yun's algorithm, whose gcds are
   :cpp:func:`varietas::modular_gcd_with_cofactors` and whose divisions are the
   cofactors those return.

.. cpp:function:: template<std::size_t N, class Order> \
                  factorisation<polynomial<rational, N, Order>> \
                  factor_univariate(const polynomial<rational, N, Order>& f)

   The factorisation into irreducibles over :math:`\Q` of an :math:`f`
   involving **at most one** of the variables, by Zassenhaus's algorithm. Each
   squarefree part is made a primitive integer polynomial and reduced modulo the
   prime, of the first five that keep it squarefree, at which it has the fewest
   factors; it is factored there by distinct-degree and Cantor–Zassenhaus
   splitting, the factors are lifted by quadratic Hensel steps down a factor tree
   to :math:`p^l > 2B`, :math:`B` being Mignotte's bound scaled by the leading
   coefficient, and recombined.

   The recombination is **exact**: a subset is accepted only when the 1-norms of
   it and its cofactor, read symmetrically modulo :math:`p^l`, multiply to at most
   :math:`B`, which forces their product to equal the polynomial over
   :math:`\mathbb{Z}`. Its cost is exponential in the number of modular factors
   when the polynomial is irreducible but splits finely modulo every prime, as the
   polynomials of Swinnerton-Dyer do; a constant-term divisibility test discards
   most candidates before any product is formed, and van Hoeij's lattice
   reduction, which removes the exponential, is not implemented.

   A polynomial in more than one variable is a precondition violation.

.. cpp:function:: template<std::size_t N, class Order> \
                  factorisation<polynomial<rational, N, Order>> \
                  factor(const polynomial<rational, N, Order>& f)

   The factorisation into irreducibles over :math:`\Q` in **any number of
   variables**, the factors monic and ordered by multiplicity and then degree.

   Each squarefree part has its content in a main variable :math:`x` factored
   with one variable fewer. :math:`x` is a variable in which the leading
   coefficient is constant if there is one, and of least degree otherwise. The
   primitive part is factored at a point: the other variables are moved so that
   a point at which the image keeps its degree and stays squarefree sits at the
   origin, and of the first three such points the one whose image has the fewest
   univariate factors is kept. A single factor there proves irreducibility.

   The univariate factors are split into two groups and both products lifted by
   Hensel's lemma, one total degree at a time, with the whole leading
   coefficient :math:`l` imposed on each, which makes :math:`l f` their product
   and leaves nothing to choose. The primitive part of a lifted group is
   accepted only if it divides :math:`f` over :math:`\mathbb{Z}`, and each group
   is then factored again inside its own factor. When a group fails to divide,
   the point has split a true factor, and the piece falls back to lifting
   :math:`f / l` as power series and recombining subsets by trial division.

   Leading coefficients are imposed rather than shared out as in Wang's
   algorithm, which costs lifting to the degree of :math:`l f` rather than of
   :math:`f`. Counters are in ``factor_counters()``.

Decomposition
=============

``varietas/codegen/decompose.hpp``

.. cpp:function:: template<std::size_t N, class Order> \
                  std::vector<variety_piece<N, Order>> \
                  decompose(const std::vector<polynomial<rational, N, Order>>& generators, \
                            decomposition_statistics* statistics = nullptr)

   :math:`\V(I)` as a union of pieces, each a reduced Gröbner basis under
   ``Order`` none of whose elements factors over :math:`\Q`, together with its
   :cpp:struct:`affine_dimension`, largest first. An empty variety gives no
   pieces.

   It is the factorising Gröbner basis algorithm. When a basis element factors as
   :math:`p_1 \cdots p_k`, the branch for :math:`p_i` is the ideal with
   :math:`p_i` adjoined and carries :math:`p_1, \dots, p_{i-1}` as polynomials
   that do not vanish on what it is responsible for. A branch is dropped when its
   ideal is the unit ideal or contains one of those, and saturated by their
   product before it is reported. A piece whose ideal contains another's is
   dropped at the end.

   The union of the pieces is :math:`\V(I)` **exactly**. The pieces are **not**
   the irreducible components: factors over extensions of :math:`\Q` are not
   found, and reducibility no single basis element shows is not seen. The
   decomposition is of the variety, so repeated factors are taken once.

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
   that evaluates them at points of :math:`\mathbb{Z}/p`. For each prime the
   degrees come from a random line by univariate rational reconstruction under
   the maximal quotient rule, and the coefficients from one linear system per
   function, normalised so the leading coefficient of the denominator is one; after the
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
