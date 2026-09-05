#ifndef VARIETAS_KINEMATICS_RIGID_TRANSFORM_HPP
#define VARIETAS_KINEMATICS_RIGID_TRANSFORM_HPP

#include <array>
#include <cstddef>
#include <utility>

#include "varietas/core/config.hpp"

namespace varietas {

// Fixed geometry of a chain, in the coefficient field the rest of the pipeline
// runs over.
//
// These are deliberately not Eigen types. The constant part of the forward
// kinematics is multiplied into polynomials whose coefficients are exact
// rationals, so the entries of a joint origin have to live in the same field as
// those coefficients; Eigen's numeric traits are built around floating point
// and would have to be taught about the rational type for no benefit at this
// size. Three by three is small enough to write out.
//
// The exactness constraint is not decorative. A rotation matrix assembled from
// roll-pitch-yaw angles has transcendental entries, and rounding them to the
// nearest rational silently perturbs the ideal: the Gröbner basis that comes
// out is then the exact basis of a robot that is not the one on the bench. The
// constructors below therefore build only rotations that are rational by
// construction, and orthogonality is checked rather than assumed.

template <class Coeff>
class vector3 {
 public:
  using traits = coefficient_traits<Coeff>;

  vector3() : e_{traits::zero(), traits::zero(), traits::zero()} {}

  vector3(Coeff x, Coeff y, Coeff z)
      : e_{std::move(x), std::move(y), std::move(z)} {}

  static vector3 zero() { return vector3(); }

  // The i-th standard basis vector, the usual URDF joint axis.
  static vector3 unit(std::size_t i) {
    VARIETAS_ASSERT(i < 3);
    vector3 v;
    v.e_[i] = traits::one();
    return v;
  }

  const Coeff& operator[](std::size_t i) const {
    VARIETAS_ASSERT(i < 3);
    return e_[i];
  }

  Coeff& operator[](std::size_t i) {
    VARIETAS_ASSERT(i < 3);
    return e_[i];
  }

  Coeff squared_norm() const {
    return e_[0] * e_[0] + e_[1] * e_[1] + e_[2] * e_[2];
  }

  bool is_zero() const {
    return traits::is_zero(e_[0]) && traits::is_zero(e_[1]) && traits::is_zero(e_[2]);
  }

  friend vector3 operator+(const vector3& a, const vector3& b) {
    return vector3(a.e_[0] + b.e_[0], a.e_[1] + b.e_[1], a.e_[2] + b.e_[2]);
  }

  friend vector3 operator-(const vector3& a, const vector3& b) {
    return vector3(a.e_[0] - b.e_[0], a.e_[1] - b.e_[1], a.e_[2] - b.e_[2]);
  }

  friend vector3 operator*(const Coeff& c, const vector3& v) {
    return vector3(c * v.e_[0], c * v.e_[1], c * v.e_[2]);
  }

  friend bool operator==(const vector3& a, const vector3& b) {
    for (std::size_t i = 0; i < 3; ++i) {
      if (!traits::is_zero(a.e_[i] - b.e_[i])) {
        return false;
      }
    }
    return true;
  }

  friend bool operator!=(const vector3& a, const vector3& b) { return !(a == b); }

 private:
  std::array<Coeff, 3> e_;
};

template <class Coeff>
class matrix3 {
 public:
  using traits = coefficient_traits<Coeff>;

  matrix3() {
    for (std::size_t i = 0; i < 9; ++i) {
      e_[i] = traits::zero();
    }
  }

  static matrix3 identity() {
    matrix3 m;
    for (std::size_t i = 0; i < 3; ++i) {
      m(i, i) = traits::one();
    }
    return m;
  }

  // Row major, so that (i, j) is the entry in row i and column j.
  static matrix3 from_rows(const std::array<Coeff, 9>& entries) {
    matrix3 m;
    m.e_ = entries;
    return m;
  }

  const Coeff& operator()(std::size_t i, std::size_t j) const {
    VARIETAS_ASSERT(i < 3 && j < 3);
    return e_[3 * i + j];
  }

  Coeff& operator()(std::size_t i, std::size_t j) {
    VARIETAS_ASSERT(i < 3 && j < 3);
    return e_[3 * i + j];
  }

  matrix3 transpose() const {
    matrix3 r;
    for (std::size_t i = 0; i < 3; ++i) {
      for (std::size_t j = 0; j < 3; ++j) {
        r(i, j) = (*this)(j, i);
      }
    }
    return r;
  }

  Coeff determinant() const {
    const matrix3& m = *this;
    return m(0, 0) * (m(1, 1) * m(2, 2) - m(1, 2) * m(2, 1)) -
           m(0, 1) * (m(1, 0) * m(2, 2) - m(1, 2) * m(2, 0)) +
           m(0, 2) * (m(1, 0) * m(2, 1) - m(1, 1) * m(2, 0));
  }

  friend matrix3 operator*(const matrix3& a, const matrix3& b) {
    matrix3 r;
    for (std::size_t i = 0; i < 3; ++i) {
      for (std::size_t j = 0; j < 3; ++j) {
        Coeff sum = traits::zero();
        for (std::size_t k = 0; k < 3; ++k) {
          sum = sum + a(i, k) * b(k, j);
        }
        r(i, j) = sum;
      }
    }
    return r;
  }

  friend vector3<Coeff> operator*(const matrix3& a, const vector3<Coeff>& v) {
    vector3<Coeff> r;
    for (std::size_t i = 0; i < 3; ++i) {
      Coeff sum = traits::zero();
      for (std::size_t k = 0; k < 3; ++k) {
        sum = sum + a(i, k) * v[k];
      }
      r[i] = sum;
    }
    return r;
  }

  friend bool operator==(const matrix3& a, const matrix3& b) {
    for (std::size_t i = 0; i < 9; ++i) {
      if (!traits::is_zero(a.e_[i] - b.e_[i])) {
        return false;
      }
    }
    return true;
  }

  friend bool operator!=(const matrix3& a, const matrix3& b) { return !(a == b); }

 private:
  std::array<Coeff, 9> e_;
};

// How far R is from being a rotation, as the largest absolute entry of
// R^T R - I, reported as a double so that one number serves both fields. Over
// the exact field a valid origin gives exactly zero and validation demands it;
// over double the same quantity is a tolerance to be compared against.
template <class Coeff>
double orthogonality_defect(const matrix3<Coeff>& r) {
  using traits = coefficient_traits<Coeff>;
  const matrix3<Coeff> gram = r.transpose() * r;
  const matrix3<Coeff> identity = matrix3<Coeff>::identity();

  double worst = 0.0;
  for (std::size_t i = 0; i < 3; ++i) {
    for (std::size_t j = 0; j < 3; ++j) {
      double d = traits::to_double(gram(i, j) - identity(i, j));
      if (d < 0.0) {
        d = -d;
      }
      if (d > worst) {
        worst = d;
      }
    }
  }
  return worst;
}

// A rotation of the exact field, from a quaternion with entries in that field.
//
// The Euler-Rodrigues matrix of (w, x, y, z), divided by the squared norm, is a
// rational function of the four entries, so any nonzero rational quaternion
// yields an exactly orthogonal rational rotation of determinant one, with no
// normalisation of the quaternion is required, and no square root is taken.
// This is the constructor a URDF front end will target: an rpy triple is
// admissible exactly when its quaternion can be written rationally.
template <class Coeff>
matrix3<Coeff> rotation_from_quaternion(const Coeff& w, const Coeff& x, const Coeff& y,
                                        const Coeff& z) {
  using traits = coefficient_traits<Coeff>;

  const Coeff norm = w * w + x * x + y * y + z * z;
  VARIETAS_ASSERT(!traits::is_zero(norm));
  const Coeff s = traits::inverse(norm);
  const Coeff two = traits::one() + traits::one();

  matrix3<Coeff> r;
  r(0, 0) = s * (w * w + x * x - y * y - z * z);
  r(0, 1) = s * (two * (x * y - w * z));
  r(0, 2) = s * (two * (x * z + w * y));
  r(1, 0) = s * (two * (x * y + w * z));
  r(1, 1) = s * (w * w - x * x + y * y - z * z);
  r(1, 2) = s * (two * (y * z - w * x));
  r(2, 0) = s * (two * (x * z - w * y));
  r(2, 1) = s * (two * (y * z + w * x));
  r(2, 2) = s * (w * w - x * x - y * y + z * z);
  return r;
}

// Rodrigues' formula about a unit axis, given the cosine and sine of the angle
// directly rather than the angle itself. The caller is responsible for
// c^2 + s^2 = 1 and for a unit axis; both are what chain validation checks.
// Right angles and their multiples, which is most of what appears in a URDF,
// are the case (c, s) in {(0, 1), (-1, 0), (0, -1)}.
template <class Coeff>
matrix3<Coeff> rotation_about_axis(const vector3<Coeff>& axis, const Coeff& c,
                                   const Coeff& s) {
  using traits = coefficient_traits<Coeff>;

  const Coeff d = traits::one() - c;
  const Coeff& ux = axis[0];
  const Coeff& uy = axis[1];
  const Coeff& uz = axis[2];

  matrix3<Coeff> r;
  r(0, 0) = c + d * ux * ux;
  r(0, 1) = d * ux * uy - s * uz;
  r(0, 2) = d * ux * uz + s * uy;
  r(1, 0) = d * uy * ux + s * uz;
  r(1, 1) = c + d * uy * uy;
  r(1, 2) = d * uy * uz - s * ux;
  r(2, 0) = d * uz * ux - s * uy;
  r(2, 1) = d * uz * uy + s * ux;
  r(2, 2) = c + d * uz * uz;
  return r;
}

// An element of SE(3) with entries in the coefficient field: the fixed part of
// a joint's placement, or the tool frame at the tip of the chain.
template <class Coeff>
class rigid_transform {
 public:
  using traits = coefficient_traits<Coeff>;

  rigid_transform()
      : rotation_(matrix3<Coeff>::identity()), translation_(vector3<Coeff>::zero()) {}

  rigid_transform(matrix3<Coeff> rotation, vector3<Coeff> translation)
      : rotation_(std::move(rotation)), translation_(std::move(translation)) {}

  static rigid_transform identity() { return rigid_transform(); }

  static rigid_transform translation_only(vector3<Coeff> t) {
    return rigid_transform(matrix3<Coeff>::identity(), std::move(t));
  }

  static rigid_transform rotation_only(matrix3<Coeff> r) {
    return rigid_transform(std::move(r), vector3<Coeff>::zero());
  }

  const matrix3<Coeff>& rotation() const noexcept { return rotation_; }
  const vector3<Coeff>& translation() const noexcept { return translation_; }

  vector3<Coeff> apply(const vector3<Coeff>& p) const {
    return rotation_ * p + translation_;
  }

  // The inverse in SE(3), which for an orthogonal rotation needs no division
  // and therefore stays in the field. Meaningless if the rotation is not
  // orthogonal, which validation is there to establish.
  rigid_transform inverse() const {
    const matrix3<Coeff> rt = rotation_.transpose();
    const Coeff minus_one = traits::negate(traits::one());
    return rigid_transform(rt, minus_one * (rt * translation_));
  }

  friend rigid_transform operator*(const rigid_transform& a, const rigid_transform& b) {
    return rigid_transform(a.rotation_ * b.rotation_,
                           a.rotation_ * b.translation_ + a.translation_);
  }

  friend bool operator==(const rigid_transform& a, const rigid_transform& b) {
    return a.rotation_ == b.rotation_ && a.translation_ == b.translation_;
  }

  friend bool operator!=(const rigid_transform& a, const rigid_transform& b) {
    return !(a == b);
  }

 private:
  matrix3<Coeff> rotation_;
  vector3<Coeff> translation_;
};

}  // namespace varietas

#endif
