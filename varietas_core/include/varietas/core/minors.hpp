#ifndef VARIETAS_CORE_MINORS_HPP
#define VARIETAS_CORE_MINORS_HPP

#include <cstddef>
#include <unordered_map>
#include <utility>
#include <vector>

#include "varietas/core/config.hpp"

namespace varietas {

// Determinants and minors of a matrix whose entries lie in a commutative ring
// rather than a field.
//
// Eigen is not the tool here and no amount of care would make it one. Every
// decomposition it offers divides — Gaussian elimination by a pivot, QR by a
// norm — and the entries this header is written for are polynomials, which are
// not invertible. The same fact the rationalisation ran into with matrix3 and
// answered with rational_transform: a determinant over k[x] must be assembled
// from additions and multiplications alone.
//
// Laplace expansion is the arrangement that never divides. Written naively it
// re-expands the same submatrix once for every path that reaches it, so the
// expansion here is memoised on the set of columns still available: a k by k
// determinant costs 2^k k multiplications rather than k!, which is the
// difference between a six by six being immediate and being a minute. What it
// buys over Bareiss fraction-free elimination is that no exact division is
// required of the entry type at all, and polynomial has none to offer.

// A dense matrix over Element, row major. Deliberately minimal: it exists to be
// handed to the two functions below, and the kinematic Jacobian is the only
// thing that builds one.
template <class Element>
class dense_matrix {
 public:
  dense_matrix() = default;

  dense_matrix(std::size_t rows, std::size_t cols)
      : rows_(rows), cols_(cols), entries_(rows * cols) {}

  std::size_t rows() const noexcept { return rows_; }
  std::size_t cols() const noexcept { return cols_; }

  const Element& operator()(std::size_t i, std::size_t j) const {
    VARIETAS_ASSERT(i < rows_ && j < cols_);
    return entries_[i * cols_ + j];
  }

  Element& operator()(std::size_t i, std::size_t j) {
    VARIETAS_ASSERT(i < rows_ && j < cols_);
    return entries_[i * cols_ + j];
  }

  // The submatrix on the given rows and columns, in the order listed.
  dense_matrix submatrix(const std::vector<std::size_t>& rows,
                         const std::vector<std::size_t>& cols) const {
    dense_matrix r(rows.size(), cols.size());
    for (std::size_t i = 0; i < rows.size(); ++i) {
      for (std::size_t j = 0; j < cols.size(); ++j) {
        r(i, j) = (*this)(rows[i], cols[j]);
      }
    }
    return r;
  }

 private:
  std::size_t rows_ = 0;
  std::size_t cols_ = 0;
  std::vector<Element> entries_;
};

namespace detail {

// Expansion along row `row`, over the columns whose bit is set in `available`.
// The cache is keyed on the column set alone, since the row is determined by
// how many columns have been consumed.
template <class Element>
const Element& memoised_determinant(const dense_matrix<Element>& m, std::size_t row,
                                    unsigned long long available,
                                    std::unordered_map<unsigned long long, Element>& cache) {
  const auto found = cache.find(available);
  if (found != cache.end()) {
    return found->second;
  }

  Element value{};  // the zero of the ring
  bool negate = false;
  for (std::size_t j = 0; j < m.cols(); ++j) {
    const unsigned long long bit = 1ull << j;
    if ((available & bit) == 0) {
      continue;
    }
    const Element& entry = m(row, j);
    if (!(entry == Element{})) {
      // On the last row the cofactor is the identity of the ring, which is
      // never formed: the entry is its own contribution.
      const Element product =
          (row + 1 == m.rows())
              ? entry
              : entry * memoised_determinant(m, row + 1, available & ~bit, cache);
      value = negate ? value - product : value + product;
    }
    negate = !negate;
  }

  return cache.emplace(available, std::move(value)).first->second;
}

}  // namespace detail

// The determinant of a square matrix over a commutative ring.
//
// Element is required to supply nothing but addition, subtraction,
// multiplication, default construction as zero and equality — no identity, no
// inverse, no division — which is exactly what a polynomial ring offers.
template <class Element>
Element determinant(const dense_matrix<Element>& m) {
  VARIETAS_ASSERT(m.rows() == m.cols());
  VARIETAS_ASSERT(m.rows() > 0 && m.rows() < 64);

  std::unordered_map<unsigned long long, Element> cache;
  const unsigned long long all = (1ull << m.cols()) - 1;
  return detail::memoised_determinant(m, 0, all, cache);
}

// Every minor of maximal size, that is of size min(rows, cols). Their common
// zero set is the locus where the matrix drops rank, which is the whole reason
// this header exists: rank is not a polynomial condition, but rank deficiency
// is, and it is cut out by exactly these.
//
// Zero minors are dropped rather than returned, since they generate nothing.
// The count is C(rows, k) * C(cols, k) determinants of size k, so the caller
// choosing which rows are the task is choosing the cost as well.
template <class Element>
std::vector<Element> maximal_minors(const dense_matrix<Element>& m) {
  const std::size_t k = m.rows() < m.cols() ? m.rows() : m.cols();
  std::vector<Element> minors;
  if (k == 0) {
    return minors;
  }

  std::vector<std::size_t> rows(k);
  std::vector<std::size_t> cols(k);

  const auto next_subset = [](std::vector<std::size_t>& subset, std::size_t n) {
    const std::size_t k_local = subset.size();
    std::size_t i = k_local;
    while (i > 0) {
      --i;
      if (subset[i] + (k_local - i) < n) {
        ++subset[i];
        for (std::size_t j = i + 1; j < k_local; ++j) {
          subset[j] = subset[j - 1] + 1;
        }
        return true;
      }
    }
    return false;
  };

  for (std::size_t i = 0; i < k; ++i) {
    rows[i] = i;
  }
  do {
    for (std::size_t i = 0; i < k; ++i) {
      cols[i] = i;
    }
    do {
      Element minor = determinant(m.submatrix(rows, cols));
      if (!(minor == Element{})) {
        minors.push_back(std::move(minor));
      }
    } while (next_subset(cols, m.cols()));
  } while (next_subset(rows, m.rows()));

  return minors;
}

}  // namespace varietas

#endif
