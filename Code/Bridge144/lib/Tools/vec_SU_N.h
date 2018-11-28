/*!
        @file    $Id:: vec_SU_N.h #$

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 00:49:46 #$

        @version $LastChangedRevision: 1561 $
*/

#ifndef VEC_SU_N_INCLUDED
#define VEC_SU_N_INCLUDED

#include <cassert>
#include "Parameters/commonParameters.h"
#include "mat_SU_N.h"

#include "bridge_complex.h"

namespace SU_N {
  class Vec_SU_N {
   private:
    int Nc;
    std::valarray<double> va;
   public:
    Vec_SU_N(int Nci = CommonParameters::Nc(), double r = 0.0)
      : Nc(Nci), va(r, 2 * Nc) // N.B. arg ordering of valarray constructor.
    {
    }

    Vec_SU_N(const Vec_SU_N& v) : Nc(v.nc()), va(v.va)
    {
    }

    Vec_SU_N& operator=(const Vec_SU_N& v)
    {
      assert(Nc == v.nc());
      va = v.va;
      return *this;
    }

    inline double norm() const;
    inline double operator*(const Vec_SU_N&) const;

    inline int nc() const
    { return Nc; }
    Vec_SU_N& dag();
    Vec_SU_N& zero();
    Vec_SU_N& xI();
    Vec_SU_N& operator-();

    Vec_SU_N& operator+=(const Vec_SU_N&);
    Vec_SU_N& operator-=(const Vec_SU_N&);
    Vec_SU_N& operator*=(const double&);
    Vec_SU_N& operator/=(const double&);
    Vec_SU_N& operator*=(const dcomplex&);
    Vec_SU_N& operator/=(const dcomplex&);

    inline int size() const
    { return va.size(); }

    inline double r(const int c) const
    { return va[2 * c]; }
    inline double i(const int c) const
    { return va[2 * c + 1]; }

    inline void set_r(const int c, const double re)
    { va[2 * c] = re; }
    inline void set_i(const int c, const double im)
    { va[2 * c + 1] = im; }
    inline void set(const int c, const double re, const double im)
    {
      va[2 * c]     = re;
      va[2 * c + 1] = im;
    }
  };

  inline double Vec_SU_N::norm() const
  {
    std::valarray<double> tmp = va * va;

    return tmp.sum();
  }


  inline double Vec_SU_N::operator*(const Vec_SU_N& rhs) const
  {
    std::valarray<double> tmp = va * rhs.va;

    return tmp.sum();
  }


  inline Vec_SU_N& Vec_SU_N::dag()
  {
    for (int c = 0; c < Nc; ++c) {
      va[2 * c + 1] = -va[2 * c + 1];
    }
    return *this;
  }


  inline Vec_SU_N& Vec_SU_N::zero()
  {
    va = 0.0;
    return *this;
  }


  inline Vec_SU_N& Vec_SU_N::xI()
  {
    for (unsigned int c = 0; c < va.size() / 2; ++c) {
      double tmp = va[2 * c];
      va[2 * c]     = -va[2 * c + 1];
      va[2 * c + 1] = tmp;
    }
    return *this;
  }


  inline Vec_SU_N& Vec_SU_N::operator-()
  {
    va = -va;
    return *this;
  }


  inline Vec_SU_N& Vec_SU_N::operator+=(const Vec_SU_N& rhs)
  {
    va += rhs.va;
    return *this;
  }


  inline Vec_SU_N& Vec_SU_N::operator-=(const Vec_SU_N& rhs)
  {
    va -= rhs.va;
    return *this;
  }


  inline Vec_SU_N& Vec_SU_N::operator*=(const double& rhs)
  {
    va *= rhs;
    return *this;
  }


  inline Vec_SU_N& Vec_SU_N::operator*=(const dcomplex& rhs)
  {
    std::valarray<double> tmp = va;

    for (unsigned int c = 0; c < va.size() / 2; ++c) {
      va[2 * c]     = (tmp[2 * c] * real(rhs) - tmp[2 * c + 1] * imag(rhs));
      va[2 * c + 1] = (tmp[2 * c] * imag(rhs) + tmp[2 * c + 1] * real(rhs));
    }
    return *this;
  }


  inline Vec_SU_N& Vec_SU_N::operator/=(const double& rhs)
  {
    va /= rhs;
    return *this;
  }


  inline Vec_SU_N& Vec_SU_N::operator/=(const dcomplex& rhs)
  {
    std::valarray<double> tmp = va;

    for (unsigned int c = 0; c < va.size() / 2; ++c) {
      va[2 * c]     = (tmp[2 * c] * real(rhs) + tmp[2 * c + 1] * imag(rhs)) / abs(rhs);
      va[2 * c + 1] = (-tmp[2 * c] * imag(rhs) - tmp[2 * c + 1] * real(rhs)) / abs(rhs);
    }
    return *this;
  }


  inline const Vec_SU_N Ix(const Vec_SU_N& u)
  {
    int      Nc = u.nc();
    Vec_SU_N tmp(Nc);

    for (int c = 0; c < u.size() / 2; ++c) {
      tmp.set(c, -u.i(c), u.r(c));
    }
    return tmp;
  }


  inline const Vec_SU_N operator+(const Vec_SU_N& v1, const Vec_SU_N& v2)
  {
    return Vec_SU_N(v1) += v2;
  }


  inline const Vec_SU_N operator-(const Vec_SU_N& v1, const Vec_SU_N& v2)
  {
    return Vec_SU_N(v1) -= v2;
  }


  inline const Vec_SU_N operator*(const Vec_SU_N& v, const double& r)
  {
    return Vec_SU_N(v) *= r;
  }


  inline const Vec_SU_N operator*(const double& r, const Vec_SU_N& v)
  {
    return Vec_SU_N(v) *= r;
  }


  inline const Vec_SU_N operator/(const Vec_SU_N& v, const double& r)
  {
    return Vec_SU_N(v) /= r;
  }


  inline const Vec_SU_N operator*(const Mat_SU_N& m, const Vec_SU_N& v)
  {
//    int Nc = CommonParameters::Nc();
    int      Nc = v.nc();
    Vec_SU_N tmp(Nc);

    for (int a = 0; a < Nc; ++a) {
      double re = 0.0;
      double im = 0.0;
      for (int b = 0; b < Nc; ++b) {
        re += m.r(a, b) * v.r(b) - m.i(a, b) * v.i(b);
        im += m.r(a, b) * v.i(b) + m.i(a, b) * v.r(b);
      }
      tmp.set(a, re, im);
    }
    return tmp;
  }
}  // namespace SU_N
#endif
