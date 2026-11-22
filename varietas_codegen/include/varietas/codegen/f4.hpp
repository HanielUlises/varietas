#ifndef VARIETAS_CODEGEN_F4_HPP
#define VARIETAS_CODEGEN_F4_HPP

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "varietas/codegen/prime_field.hpp"
#include "varietas/core/config.hpp"
#include "varietas/core/ideal/buchberger.hpp"
#include "varietas/core/monomial.hpp"
#include "varietas/core/polynomial.hpp"

namespace varietas {

// Faugère's F4 over Z/p.
//
// Buchberger reduces one S-polynomial at a time, and each reduction is a walk
// through polynomials term by term, finding a divisor, scaling it, merging it
// in. F4 does the same arithmetic in bulk. It takes every pending pair of the
// lowest sugar degree at once, collects every multiple of a basis element that
// reducing them could call on (the symbolic preprocessing), writes the lot as
// the rows of one sparse matrix whose columns are the monomials in decreasing
// order, and row-reduces it. The rows whose leading monomial no basis element
// divides are the new basis elements, and the rows that vanish are the
// reductions to zero, found without anyone having asked which pairs they were.
//
// Over a prime field below 2^31 a row operation is a loop over machine words,
// which is what makes the matrix form pay: the reconstruction
// (varietas/codegen/reconstruct.hpp) computes one saturated basis per sample,
// hundreds of times per prime, and every one of them is over Z/p.
//
// The pairs are pruned by the criteria of Gebauer and Möller, which subsume
// both of Buchberger's: a pair is dropped when a third element certifies in
// advance that its S-polynomial reduces to zero.
//
// The result is the reduced Gröbner basis, so it is the same object Buchberger
// computes, element for element; nothing downstream can tell which algorithm
// produced it, and the tests check that it cannot.
namespace f4 {

namespace detail {

using word = residue::word;

template <std::size_t N>
struct monomial_hash {
  std::size_t operator()(const monomial<N>& m) const noexcept {
    std::size_t h = 1469598103934665603ull;
    for (std::size_t i = 0; i < N; ++i) {
      h = (h ^ m[i]) * 1099511628211ull;
    }
    return h;
  }
};

// A divisibility mask: for each variable a few bits, bit k set when the
// exponent exceeds k. If a divides b every bit of a's mask is set in b's, so a
// mask with a bit b's lacks rejects a divisor in one instruction, which is
// most of the questions the symbolic preprocessing asks.
template <std::size_t N>
std::uint64_t divisibility_mask(const monomial<N>& m) noexcept {
  constexpr std::size_t per = N == 0 ? 1 : (64 / N < 1 ? 1 : 64 / N);
  constexpr std::size_t bits = per > 8 ? 8 : per;
  std::uint64_t mask = 0;
  for (std::size_t v = 0; v < N && v * bits < 64; ++v) {
    const std::size_t e = m[v] < bits ? m[v] : bits;
    mask |= ((std::uint64_t{1} << e) - 1) << (v * bits);
  }
  return mask;
}

// A basis element, monic, its terms in decreasing order.
template <std::size_t N>
struct element {
  std::vector<monomial<N>> monomials;
  std::vector<word> coefficients;
  unsigned sugar = 0;
  bool redundant = false;

  const monomial<N>& leading() const { return monomials.front(); }
};

template <std::size_t N>
struct critical_pair {
  std::uint32_t i;
  std::uint32_t j;
  monomial<N> lcm;
  unsigned sugar;
};

// A row before it has columns: a basis element times a monomial.
template <std::size_t N>
struct row_spec {
  std::uint32_t element;
  monomial<N> multiplier;
};

template <std::size_t N>
struct row_spec_hash {
  std::size_t operator()(const row_spec<N>& r) const noexcept {
    return monomial_hash<N>()(r.multiplier) * 31u + r.element;
  }
};

template <std::size_t N>
bool operator==(const row_spec<N>& a, const row_spec<N>& b) noexcept {
  return a.element == b.element && a.multiplier == b.multiplier;
}

struct sparse_row {
  std::vector<std::uint32_t> columns;  // increasing, that is decreasing monomials
  std::vector<word> values;
};

template <std::size_t N, class Order>
class engine {
 public:
  using mono = monomial<N>;
  using poly = polynomial<residue, N, Order>;

  explicit engine(buchberger_statistics* statistics)
      : p_(residue::prime()), statistics_(statistics) {}

  // False when the ideal is the whole ring, which is found as soon as a constant
  // turns up and needs no further work.
  bool run(const std::vector<poly>& generators) {
    for (const poly& f : generators) {
      if (f.is_zero()) {
        continue;
      }
      element<N> e;
      const residue inverse = coefficient_traits<residue>::inverse(f.leading_coefficient());
      for (const auto& t : f.terms()) {
        e.monomials.push_back(t.mon);
        e.coefficients.push_back((t.coeff * inverse).value());
      }
      e.sugar = f.degree();
      if (e.leading().is_one()) {
        return false;
      }
      insert(std::move(e));
    }
    while (!pairs_.empty()) {
      if (!round()) {
        return false;
      }
    }
    return true;
  }

  std::vector<poly> basis() const {
    std::vector<poly> out;
    for (const element<N>& e : basis_) {
      if (e.redundant) {
        continue;
      }
      std::vector<typename poly::term> terms;
      terms.reserve(e.monomials.size());
      for (std::size_t k = 0; k < e.monomials.size(); ++k) {
        terms.push_back({e.monomials[k], residue::raw(e.coefficients[k])});
      }
      out.emplace_back(std::move(terms));
    }
    return out;
  }

 private:
  void count(std::size_t buchberger_statistics::*field, std::size_t n = 1) {
    if (statistics_ != nullptr) {
      statistics_->*field += n;
    }
  }

  // The update of Gebauer and Möller, as Becker and Weispfenning give it, for a
  // new element h: the pairs h forms that survive, then the old pairs h makes
  // redundant removed, then the old elements h makes redundant marked.
  void insert(element<N> h) {
    const std::uint32_t hi = static_cast<std::uint32_t>(basis_.size());
    const mono& lh = h.leading();

    std::vector<critical_pair<N>> formed;
    for (std::uint32_t g = 0; g < hi; ++g) {
      if (basis_[g].redundant) {
        continue;
      }
      const element<N>& e = basis_[g];
      const mono l = mono::lcm(lh, e.leading());
      const unsigned sugar = std::max(h.sugar + l.degree() - lh.degree(),
                                      e.sugar + l.degree() - e.leading().degree());
      formed.push_back(critical_pair<N>{g, hi, l, sugar});
    }
    count(&buchberger_statistics::pairs_generated, formed.size());

    const auto coprime = [&](const critical_pair<N>& c) {
      return c.lcm.degree() == lh.degree() + basis_[c.i].leading().degree();
    };

    // A pair (h, g1) goes when another (h, g2) has an lcm properly dividing its
    // own, or an equal one and comes first; a coprime pair is kept for this
    // step, since it may be the one that removes others, and dropped after.
    std::vector<critical_pair<N>> kept;
    for (std::size_t a = 0; a < formed.size(); ++a) {
      const critical_pair<N>& c = formed[a];
      if (coprime(c)) {
        kept.push_back(c);
        continue;
      }
      bool dominated = false;
      for (std::size_t b = 0; b < formed.size() && !dominated; ++b) {
        if (b == a || !mono::divides(formed[b].lcm, c.lcm)) {
          continue;
        }
        dominated = formed[b].lcm != c.lcm || b < a;
      }
      if (dominated) {
        count(&buchberger_statistics::pairs_discarded_chain);
      } else {
        kept.push_back(c);
      }
    }

    // Old pairs whose lcm h's leading monomial divides, and which h would
    // reduce through, are covered by the two pairs h forms with their ends.
    std::vector<critical_pair<N>> old;
    old.reserve(pairs_.size());
    for (const critical_pair<N>& c : pairs_) {
      if (mono::divides(lh, c.lcm) && mono::lcm(basis_[c.i].leading(), lh) != c.lcm &&
          mono::lcm(basis_[c.j].leading(), lh) != c.lcm) {
        count(&buchberger_statistics::pairs_discarded_chain);
        continue;
      }
      old.push_back(c);
    }
    pairs_ = std::move(old);
    for (const critical_pair<N>& c : kept) {
      if (coprime(c)) {
        count(&buchberger_statistics::pairs_discarded_coprime);
      } else {
        pairs_.push_back(c);
      }
    }

    for (element<N>& e : basis_) {
      if (!e.redundant && mono::divides(lh, e.leading())) {
        e.redundant = true;
      }
    }
    active_.erase(std::remove_if(active_.begin(), active_.end(),
                                 [&](const active_element& a) { return basis_[a.index].redundant; }),
                  active_.end());
    active_.push_back(active_element{hi, divisibility_mask(lh), h.monomials.size()});
    basis_.push_back(std::move(h));
  }

  // A basis element whose leading monomial divides m, the shortest of them.
  std::int64_t reducer_for(const mono& m) const {
    const std::uint64_t mask = divisibility_mask(m);
    std::int64_t best = -1;
    std::size_t shortest = 0;
    for (const active_element& a : active_) {
      if ((a.mask & ~mask) != 0 || (best >= 0 && a.length >= shortest) ||
          !mono::divides(basis_[a.index].leading(), m)) {
        continue;
      }
      best = static_cast<std::int64_t>(a.index);
      shortest = a.length;
    }
    return best;
  }

  bool round() {
    // The pairs of least sugar, all of them.
    unsigned lowest = pairs_.front().sugar;
    for (const critical_pair<N>& c : pairs_) {
      lowest = std::min(lowest, c.sugar);
    }
    std::vector<critical_pair<N>> selected;
    std::vector<critical_pair<N>> rest;
    for (const critical_pair<N>& c : pairs_) {
      (c.sugar == lowest ? selected : rest).push_back(c);
    }
    pairs_ = std::move(rest);

    // Rows: both halves of every selected pair, each multiple once.
    std::vector<row_spec<N>> specs;
    std::vector<std::vector<std::uint32_t>> row_columns;  // in first-seen numbering
    std::vector<bool> is_reducer;
    std::unordered_set<row_spec<N>, row_spec_hash<N>> seen;
    std::unordered_map<mono, std::uint32_t, monomial_hash<N>> column_of;
    std::vector<mono> columns;
    std::vector<bool> done;  // a column that already has a row leading there
    std::vector<std::uint32_t> pending;

    const auto add_monomial = [&](const mono& m) -> std::uint32_t {
      const auto [it, fresh] = column_of.emplace(m, static_cast<std::uint32_t>(columns.size()));
      if (fresh) {
        columns.push_back(m);
        done.push_back(false);
        pending.push_back(it->second);
      }
      return it->second;
    };
    const auto add_row = [&](std::uint32_t g, const mono& u, bool reducer) {
      row_spec<N> spec{g, u};
      if (!seen.insert(spec).second) {
        return;
      }
      specs.push_back(spec);
      is_reducer.push_back(reducer);
      std::vector<std::uint32_t> ids;
      ids.reserve(basis_[g].monomials.size());
      for (const mono& m : basis_[g].monomials) {
        ids.push_back(add_monomial(m * u));
      }
      row_columns.push_back(std::move(ids));
    };

    for (const critical_pair<N>& c : selected) {
      add_row(c.i, mono::divide(c.lcm, basis_[c.i].leading()), false);
      add_row(c.j, mono::divide(c.lcm, basis_[c.j].leading()), false);
      done[column_of.at(c.lcm)] = true;
    }

    // Symbolic preprocessing: every monomial some basis element can reduce gets
    // one reducing row, whose own monomials may call for more.
    while (!pending.empty()) {
      const std::uint32_t col = pending.back();
      pending.pop_back();
      if (done[col]) {
        continue;
      }
      done[col] = true;
      const mono m = columns[col];
      const std::int64_t g = reducer_for(m);
      if (g >= 0) {
        const auto gi = static_cast<std::uint32_t>(g);
        add_row(gi, mono::divide(m, basis_[gi].leading()), true);
      }
    }

    // Columns in decreasing order.
    const std::size_t width = columns.size();
    std::vector<std::uint32_t> order(width);
    for (std::uint32_t c = 0; c < width; ++c) {
      order[c] = c;
    }
    std::sort(order.begin(), order.end(), [&](std::uint32_t a, std::uint32_t b) {
      return Order::compare(columns[a], columns[b]) > 0;
    });
    std::vector<std::uint32_t> position(width);
    for (std::uint32_t k = 0; k < width; ++k) {
      position[order[k]] = k;
    }

    // Multiplying by a monomial preserves the order, so each row's columns
    // come out already increasing.
    std::vector<sparse_row> rows(specs.size());
    for (std::size_t r = 0; r < specs.size(); ++r) {
      sparse_row& row = rows[r];
      row.values = basis_[specs[r].element].coefficients;
      row.columns = std::move(row_columns[r]);
      for (std::uint32_t& c : row.columns) {
        c = position[c];
      }
    }

    // A pivot per leading column: every reducer, and for each lcm the
    // shortest of the rows leading there. The others are what get reduced.
    std::vector<std::int64_t> pivot(width, -1);
    for (std::size_t r = 0; r < rows.size(); ++r) {
      if (is_reducer[r]) {
        pivot[rows[r].columns.front()] = static_cast<std::int64_t>(r);
      }
    }
    std::vector<std::size_t> to_reduce;
    for (std::size_t r = 0; r < rows.size(); ++r) {
      if (is_reducer[r]) {
        continue;
      }
      std::int64_t& slot = pivot[rows[r].columns.front()];
      if (slot < 0) {
        slot = static_cast<std::int64_t>(r);
      } else if (rows[r].columns.size() < rows[static_cast<std::size_t>(slot)].columns.size()) {
        to_reduce.push_back(static_cast<std::size_t>(slot));
        slot = static_cast<std::int64_t>(r);
      } else {
        to_reduce.push_back(r);
      }
    }
    count(&buchberger_statistics::pairs_reduced, to_reduce.size());

    // Reduction, one row at a time, in a dense accumulator held below p^2 so
    // that a product can be added without reducing it first.
    const word p = p_;
    const word p2 = p * p;
    std::vector<word> dense(width, 0);
    std::vector<sparse_row> produced;
    std::sort(to_reduce.begin(), to_reduce.end(), [&](std::size_t a, std::size_t b) {
      return rows[a].columns.front() < rows[b].columns.front();
    });
    for (const std::size_t r : to_reduce) {
      const sparse_row& source = rows[r];
      for (std::size_t k = 0; k < source.columns.size(); ++k) {
        dense[source.columns[k]] = source.values[k];
      }
      sparse_row out;
      for (std::uint32_t c = source.columns.front(); c < width; ++c) {
        if (dense[c] == 0) {
          continue;
        }
        const word v = dense[c] % p;
        dense[c] = 0;
        if (v == 0) {
          continue;
        }
        const std::int64_t pr = pivot[c];
        if (pr < 0) {
          out.columns.push_back(c);
          out.values.push_back(v);
          continue;
        }
        const sparse_row& by = static_cast<std::size_t>(pr) < rows.size()
                                   ? rows[static_cast<std::size_t>(pr)]
                                   : produced[static_cast<std::size_t>(pr) - rows.size()];
        const word factor = p - v;  // pivots are monic
        for (std::size_t k = 1; k < by.columns.size(); ++k) {
          word& slot = dense[by.columns[k]];
          slot += factor * by.values[k];
          if (slot >= p2) {
            slot -= p2;
          }
        }
      }
      if (out.columns.empty()) {
        count(&buchberger_statistics::reductions_to_zero);
        continue;
      }
      const word inverse = residue::inverse_of(out.values.front());
      for (word& v : out.values) {
        v = v * inverse % p;
      }
      pivot[out.columns.front()] = static_cast<std::int64_t>(rows.size() + produced.size());
      produced.push_back(std::move(out));
    }

    for (const sparse_row& row : produced) {
      element<N> e;
      e.sugar = lowest;
      e.monomials.reserve(row.columns.size());
      for (const std::uint32_t c : row.columns) {
        e.monomials.push_back(columns[order[c]]);
      }
      e.coefficients = row.values;
      if (e.leading().is_one()) {
        return false;
      }
      insert(std::move(e));
    }
    return true;
  }

  // The elements not made redundant, with what the search for a reducer needs
  // of each, kept together so that the search walks one array.
  struct active_element {
    std::uint32_t index;
    std::uint64_t mask;
    std::size_t length;
  };

  word p_;
  buchberger_statistics* statistics_;
  std::vector<element<N>> basis_;
  std::vector<active_element> active_;
  std::vector<critical_pair<N>> pairs_;
};

}  // namespace detail

// A Gröbner basis of the ideal the generators span, over the prime currently
// selected, not necessarily reduced. {1} when the ideal is the whole ring.
template <std::size_t N, class Order>
std::vector<polynomial<residue, N, Order>> basis(
    const std::vector<polynomial<residue, N, Order>>& generators,
    buchberger_statistics* statistics = nullptr) {
  detail::engine<N, Order> run(statistics);
  if (!run.run(generators)) {
    return {polynomial<residue, N, Order>::constant(residue::raw(1))};
  }
  return run.basis();
}

}  // namespace f4

// Every Gröbner basis over Z/p is computed by F4, through the engine hook in
// varietas/core/ideal/buchberger.hpp.
template <>
struct groebner_engine<residue> {
  template <std::size_t N, class Order>
  static std::vector<polynomial<residue, N, Order>> basis(
      const std::vector<polynomial<residue, N, Order>>& generators,
      buchberger_statistics* statistics) {
    return f4::basis(generators, statistics);
  }
};

}  // namespace varietas

#endif
