/*!
        @file    $Id:: mat_SU_N.h #$

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 18:35:38 #$

        @version $LastChangedRevision: 1571 $
*/


#ifndef MAT_SU_N_INCLUDED
#define MAT_SU_N_INCLUDED

#include <cassert>
#include <valarray>

#include "decompose_QR_Cmplx.h"
#include "eigen_QR_Cmplx.h"

#include "randomNumbers.h"

class RandomNumbers;

namespace SU_N {
/*!
  SU(N) matrix operations.
  Reunitarization is generalized to general value of Nc.
  Setting random matrix is not for SU(2) which is to be implemented.
                                       [07 May 2014 H.Matsufuru]
*/
  class Mat_SU_N {
   private:
    int m_Nc;
    std::valarray<double> va;
    Mat_SU_N&             (Mat_SU_N::*m_reunit)(); //! pointer to reunitalization.
    Mat_SU_N&             (SU_N::Mat_SU_N::*m_set_random)(RandomNumbers *);
    //!< pointer to random matrix setting.

   public:

    explicit
    Mat_SU_N(int Nc, double r = 0.0) : m_Nc(Nc)
    {
      va.resize(2 * m_Nc * m_Nc, r);

      if (m_Nc == 3) {
        m_reunit     = &Mat_SU_N::reunit_SU3;
        m_set_random = &SU_N::Mat_SU_N::set_random_SU3;
      } else if (m_Nc == 2) {
        m_reunit     = &Mat_SU_N::reunit_SU2;
        m_set_random = &SU_N::Mat_SU_N::set_random_SU2;
      } else {
        m_reunit     = &Mat_SU_N::reunit_general;
        m_set_random = &SU_N::Mat_SU_N::set_random_general;
      }
    }

    Mat_SU_N& dag();
    Mat_SU_N& ht();
    Mat_SU_N& ah();
    Mat_SU_N& at();
    Mat_SU_N& unit();
    Mat_SU_N& zero();
    Mat_SU_N& xI();

    //    Mat_SU_N& reunit();
    Mat_SU_N& reunit() { return (this->*m_reunit)(); }
    Mat_SU_N& reunit_SU3();
    Mat_SU_N& reunit_SU2();
    Mat_SU_N& reunit_general();

    Mat_SU_N& set_random(RandomNumbers *rand)
    {
      return (this->*m_set_random)(rand);
    }

    // H.Matsufuru added (original version for SU(3)) [15 Mar 2012].
    Mat_SU_N& set_random_SU3(RandomNumbers *);
    Mat_SU_N& set_random_SU2(RandomNumbers *);
    Mat_SU_N& set_random_general(RandomNumbers *);

    int nc() const { return m_Nc; }

    Mat_SU_N& operator-();

    Mat_SU_N& operator=(const double&);

    //  Mat_SU_N& operator=(const std::complex<double>&);

    Mat_SU_N& operator+=(const Mat_SU_N&);
    Mat_SU_N& operator+=(const double&);

    //  Mat_SU_N& operator+=(const std::complex<double>&);

    Mat_SU_N& operator-=(const Mat_SU_N&);
    Mat_SU_N& operator-=(const double&);

    //  Mat_SU_N& operator-=(const std::complex<double>&);

    Mat_SU_N& operator*=(const Mat_SU_N&);
    Mat_SU_N& operator*=(const double&);

    //  Mat_SU_N& operator*=(const std::complex<double>&);

    Mat_SU_N& operator/=(const double&);

    //  Mat_SU_N& operator/=(const std::complex<double>&);

    int size() const { return va.size(); }
    double r(int c) const { return va[2 * c]; }
    double i(int c) const { return va[2 * c + 1]; }

    double r(int c1, int c2) const { return r(m_Nc * c1 + c2); }
    double i(int c1, int c2) const { return i(m_Nc * c1 + c2); }

    void set_r(int c, const double& re) { va[2 * c] = re; }
    void set_i(int c, const double& im) { va[2 * c + 1] = im; }

    void set_r(int c1, int c2, const double& re)
    {
      set_r(m_Nc * c1 + c2, re);
    }

    void set_i(int c1, int c2, const double& im)
    {
      set_i(m_Nc * c1 + c2, im);
    }

    void set(int c, double re, const double& im)
    {
      va[2 * c]     = re;
      va[2 * c + 1] = im;
    }

    void set(int c1, int c2, const double& re, const double& im)
    {
      set(m_Nc * c1 + c2, re, im);
    }

    void add(int c, const double& re, const double& im)
    {
      va[2 * c]     += re;
      va[2 * c + 1] += im;
    }

    void add(int c1, int c2, const double& re, const double& im)
    {
      add(m_Nc * c1 + c2, re, im);
    }

    double norm2()     // H.Matsufuru added [15 Mar 2012]
    {
      double t = 0.0;

      for (int i = 0; i < 2 * m_Nc * m_Nc; ++i) { t += va[i] * va[i]; }
      return t;
    }

    void mult_nn(const Mat_SU_N& u1, const Mat_SU_N& u2)
    {
      for (int a = 0; a < m_Nc; ++a) {
        for (int b = 0; b < m_Nc; ++b) {
          va[2 * (b + m_Nc * a)]     = 0.0;
          va[2 * (b + m_Nc * a) + 1] = 0.0;
          for (int c = 0; c < m_Nc; ++c) {
            va[2 * (b + m_Nc * a)] +=
              u1.va[2 * (c + m_Nc * a)] * u2.va[2 * (b + m_Nc * c)]
              - u1.va[2 * (c + m_Nc * a) + 1] * u2.va[2 * (b + m_Nc * c) + 1];
            va[2 * (b + m_Nc * a) + 1] +=
              u1.va[2 * (c + m_Nc * a) + 1] * u2.va[2 * (b + m_Nc * c)]
              + u1.va[2 * (c + m_Nc * a)] * u2.va[2 * (b + m_Nc * c) + 1];
          }
        }
      }
    }

    void multadd_nn(const Mat_SU_N& u1, const Mat_SU_N& u2)
    {
      for (int a = 0; a < m_Nc; ++a) {
        for (int b = 0; b < m_Nc; ++b) {
          for (int c = 0; c < m_Nc; ++c) {
            va[2 * (b + m_Nc * a)] +=
              u1.va[2 * (c + m_Nc * a)] * u2.va[2 * (b + m_Nc * c)]
              - u1.va[2 * (c + m_Nc * a) + 1] * u2.va[2 * (b + m_Nc * c) + 1];
            va[2 * (b + m_Nc * a) + 1] +=
              u1.va[2 * (c + m_Nc * a) + 1] * u2.va[2 * (b + m_Nc * c)]
              + u1.va[2 * (c + m_Nc * a)] * u2.va[2 * (b + m_Nc * c) + 1];
          }
        }
      }
    }

    void mult_nd(const Mat_SU_N& u1, const Mat_SU_N& u2)
    {
      for (int a = 0; a < m_Nc; ++a) {
        for (int b = 0; b < m_Nc; ++b) {
          va[2 * (b + m_Nc * a)]     = 0.0;
          va[2 * (b + m_Nc * a) + 1] = 0.0;
          for (int c = 0; c < m_Nc; ++c) {
            va[2 * (b + m_Nc * a)] +=
              u1.va[2 * (c + m_Nc * a)] * u2.va[2 * (c + m_Nc * b)]
              + u1.va[2 * (c + m_Nc * a) + 1] * u2.va[2 * (c + m_Nc * b) + 1];
            va[2 * (b + m_Nc * a) + 1] +=
              u1.va[2 * (c + m_Nc * a) + 1] * u2.va[2 * (c + m_Nc * b)]
              - u1.va[2 * (c + m_Nc * a)] * u2.va[2 * (c + m_Nc * b) + 1];
          }
        }
      }
    }

    void multadd_nd(const Mat_SU_N& u1, const Mat_SU_N& u2)
    {
      for (int a = 0; a < m_Nc; ++a) {
        for (int b = 0; b < m_Nc; ++b) {
          for (int c = 0; c < m_Nc; ++c) {
            va[2 * (b + m_Nc * a)] +=
              u1.va[2 * (c + m_Nc * a)] * u2.va[2 * (c + m_Nc * b)]
              + u1.va[2 * (c + m_Nc * a) + 1] * u2.va[2 * (c + m_Nc * b) + 1];
            va[2 * (b + m_Nc * a) + 1] +=
              u1.va[2 * (c + m_Nc * a) + 1] * u2.va[2 * (c + m_Nc * b)]
              - u1.va[2 * (c + m_Nc * a)] * u2.va[2 * (c + m_Nc * b) + 1];
          }
        }
      }
    }

    void mult_dn(const Mat_SU_N& u1, const Mat_SU_N& u2)
    {
      for (int a = 0; a < m_Nc; ++a) {
        for (int b = 0; b < m_Nc; ++b) {
          va[2 * (b + m_Nc * a)]     = 0.0;
          va[2 * (b + m_Nc * a) + 1] = 0.0;
          for (int c = 0; c < m_Nc; ++c) {
            va[2 * (b + m_Nc * a)] +=
              u1.va[2 * (a + m_Nc * c)] * u2.va[2 * (b + m_Nc * c)]
              + u1.va[2 * (a + m_Nc * c) + 1] * u2.va[2 * (b + m_Nc * c) + 1];
            va[2 * (b + m_Nc * a) + 1] -=
              u1.va[2 * (a + m_Nc * c) + 1] * u2.va[2 * (b + m_Nc * c)]
              - u1.va[2 * (a + m_Nc * c)] * u2.va[2 * (b + m_Nc * c) + 1];
          }
        }
      }
    }

    void multadd_dn(const Mat_SU_N& u1, const Mat_SU_N& u2)
    {
      for (int a = 0; a < m_Nc; ++a) {
        for (int b = 0; b < m_Nc; ++b) {
          for (int c = 0; c < m_Nc; ++c) {
            va[2 * (b + m_Nc * a)] +=
              u1.va[2 * (a + m_Nc * c)] * u2.va[2 * (b + m_Nc * c)]
              + u1.va[2 * (a + m_Nc * c) + 1] * u2.va[2 * (b + m_Nc * c) + 1];
            va[2 * (b + m_Nc * a) + 1] -=
              u1.va[2 * (a + m_Nc * c) + 1] * u2.va[2 * (b + m_Nc * c)]
              - u1.va[2 * (a + m_Nc * c)] * u2.va[2 * (b + m_Nc * c) + 1];
          }
        }
      }
    }

    void zcopy(double re, double im, const Mat_SU_N& v)
    {
      for (int cc = 0; cc < m_Nc * m_Nc; ++cc) {
        va[2 * cc]     = re * v.va[2 * cc] - im * v.va[2 * cc + 1];
        va[2 * cc + 1] = re * v.va[2 * cc + 1] + im * v.va[2 * cc];
      }
    }

    void zaxpy(double re, double im, const Mat_SU_N& v)
    {
      for (int cc = 0; cc < m_Nc * m_Nc; ++cc) {
        va[2 * cc]     += re * v.va[2 * cc] - im * v.va[2 * cc + 1];
        va[2 * cc + 1] += re * v.va[2 * cc + 1] + im * v.va[2 * cc];
      }
    }
  }; // end of class definition.


  inline Mat_SU_N& Mat_SU_N::dag()
  {
    for (int a = 0; a < m_Nc; ++a) {
      for (int b = a; b < m_Nc; ++b) {
        double re = va[2 * (m_Nc * a + b)];
        double im = va[2 * (m_Nc * a + b) + 1];

        va[2 * (m_Nc * a + b)]     = va[2 * (m_Nc * b + a)];
        va[2 * (m_Nc * a + b) + 1] = -va[2 * (m_Nc * b + a) + 1];

        va[2 * (m_Nc * b + a)]     = re;
        va[2 * (m_Nc * b + a) + 1] = -im;
      }
    }
    return *this;
  }


  inline Mat_SU_N& Mat_SU_N::ht() // hermitian traceless
  {
    for (int a = 0; a < m_Nc; ++a) {
      for (int b = a + 1; b < m_Nc; ++b) {
        double re = va[2 * (m_Nc * a + b)] + va[2 * (m_Nc * b + a)];
        double im = va[2 * (m_Nc * a + b) + 1] - va[2 * (m_Nc * b + a) + 1];

        va[2 * (m_Nc * a + b)]     = 0.5 * re;
        va[2 * (m_Nc * a + b) + 1] = 0.5 * im;

        va[2 * (m_Nc * b + a)]     = 0.5 * re;
        va[2 * (m_Nc * b + a) + 1] = -0.5 * im;
      }
    }
    double tr = 0.0;
    for (int cc = 0; cc < m_Nc; ++cc) {
      tr += va[2 * (m_Nc * cc + cc)];
    }
    tr = tr / m_Nc;
    for (int cc = 0; cc < m_Nc; ++cc) {
      va[2 * (m_Nc * cc + cc)]    -= tr;
      va[2 * (m_Nc * cc + cc) + 1] = 0.0;
    }
    return *this;
  }


  //! antihermitian traceless
  inline Mat_SU_N& Mat_SU_N::at()
  {
    for (int a = 0; a < m_Nc; ++a) {
      for (int b = a + 1; b < m_Nc; ++b) {
        double re = va[2 * (m_Nc * a + b)] - va[2 * (m_Nc * b + a)];
        double im = va[2 * (m_Nc * a + b) + 1] + va[2 * (m_Nc * b + a) + 1];

        va[2 * (m_Nc * a + b)]     = 0.5 * re;
        va[2 * (m_Nc * a + b) + 1] = 0.5 * im;

        va[2 * (m_Nc * b + a)]     = -0.5 * re;
        va[2 * (m_Nc * b + a) + 1] = 0.5 * im;
      }
    }
    double tr = 0.0;
    for (int cc = 0; cc < m_Nc; ++cc) {
      tr += va[2 * (m_Nc * cc + cc) + 1];
    }
    tr = tr / m_Nc;
    for (int cc = 0; cc < m_Nc; ++cc) {
      va[2 * (m_Nc * cc + cc)]      = 0.0;
      va[2 * (m_Nc * cc + cc) + 1] -= tr;
    }
    return *this;
  }


  //! antihermitian
  inline Mat_SU_N& Mat_SU_N::ah()
  {
    for (int a = 0; a < m_Nc; ++a) {
      for (int b = a; b < m_Nc; ++b) {
        double re = va[2 * (m_Nc * a + b)] - va[2 * (m_Nc * b + a)];
        double im = va[2 * (m_Nc * a + b) + 1] + va[2 * (m_Nc * b + a) + 1];
        va[2 * (m_Nc * a + b)]     = 0.5 * re;
        va[2 * (m_Nc * a + b) + 1] = 0.5 * im;
        va[2 * (m_Nc * b + a)]     = -0.5 * re;
        va[2 * (m_Nc * b + a) + 1] = 0.5 * im;
      }
    }
    return *this;
  }


  inline Mat_SU_N& Mat_SU_N::unit()
  {
    va = 0.0;
    for (int c = 0; c < m_Nc; ++c) {
      va[2 * (m_Nc + 1) * c] = 1.0;
    }
    return *this;
  }


  inline Mat_SU_N& Mat_SU_N::zero()
  {
    va = 0.0;
    return *this;
  }


  inline Mat_SU_N& Mat_SU_N::xI()
  {
    for (int c = 0; c < va.size() / 2; ++c) {
      double tmp = va[2 * c];
      va[2 * c]     = -va[2 * c + 1];
      va[2 * c + 1] = tmp;
    }
    return *this;
  }


  inline Mat_SU_N& Mat_SU_N::operator-()
  {
    va = -va;
    return *this;
  }


  inline Mat_SU_N& Mat_SU_N::operator+=(const Mat_SU_N& rhs)
  {
    va += rhs.va;
    return *this;
  }


  inline Mat_SU_N& Mat_SU_N::operator-=(const Mat_SU_N& rhs)
  {
    va -= rhs.va;
    return *this;
  }


  inline Mat_SU_N& Mat_SU_N::operator=(const double& rhs)
  {
    va = rhs;
    return *this;
  }


  inline Mat_SU_N& Mat_SU_N::operator+=(const double& rhs)
  {
    va += rhs;
    return *this;
  }


  inline Mat_SU_N& Mat_SU_N::operator-=(const double& rhs)
  {
    va -= rhs;
    return *this;
  }


  inline Mat_SU_N& Mat_SU_N::operator*=(const Mat_SU_N& rhs)
  {
    std::valarray<double> tmp(2 * m_Nc * m_Nc);

    for (int a = 0; a < m_Nc; ++a) {
      for (int b = 0; b < m_Nc; ++b) {
        tmp[2 * (m_Nc * a + b)]     = 0.0;
        tmp[2 * (m_Nc * a + b) + 1] = 0.0;
        for (int c = 0; c < m_Nc; ++c) {
          tmp[2 * (m_Nc * a + b)] += va[2 * (m_Nc * a + c)] * rhs.va[2 * (m_Nc * c + b)]
                                     - va[2 * (m_Nc * a + c) + 1] * rhs.va[2 * (m_Nc * c + b) + 1];
          tmp[2 * (m_Nc * a + b) + 1] += va[2 * (m_Nc * a + c) + 1] * rhs.va[2 * (m_Nc * c + b)]
                                         + va[2 * (m_Nc * a + c)] * rhs.va[2 * (m_Nc * c + b) + 1];
        }
      }
    }
    va = tmp;
    return *this;
  }


  inline Mat_SU_N& Mat_SU_N::operator*=(const double& rhs)
  {
    va *= rhs;
    return *this;
  }


/*
inline Mat_SU_N& Mat_SU_N::operator*=(const std::complex<double>& rhs){
  std::valarray<double> tmp = va;
  for(int c = 0; c < va.size()/2; ++c){
    va[2*c  ] = (tmp[2*c]*real(rhs) -tmp[2*c+1]*imag(rhs));
    va[2*c+1] = (tmp[2*c]*imag(rhs) +tmp[2*c+1]*real(rhs));
  }
  return *this;
}
*/
  inline Mat_SU_N& Mat_SU_N::operator/=(const double& rhs)
  {
    va /= rhs;
    return *this;
  }


  inline double ReTr(const Mat_SU_N& m)
  {
    int    Nc = m.nc();
    double tr = 0.0;

    for (int c = 0; c < Nc; ++c) {
      tr += m.r(c, c);
    }
    return tr;
  }


  inline double ImTr(const Mat_SU_N& m)
  {
    int    Nc = m.nc();
    double tr = 0.0;

    for (int c = 0; c < Nc; ++c) {
      tr += m.i(c, c);
    }
    return tr;
  }


  inline Mat_SU_N dag(const Mat_SU_N& u)
  {
    int      Nc = u.nc();
    Mat_SU_N tmp(Nc);

    for (int a = 0; a < Nc; a++) {
      for (int b = 0; b < Nc; b++) {
        tmp.set(a, b, u.r(b, a), -u.i(b, a));
      }
    }
    return tmp;
  }


  inline Mat_SU_N xI(const Mat_SU_N& u)
  {
    int      Nc = u.nc();
    Mat_SU_N tmp(Nc);

    for (int c = 0; c < u.size() / 2; ++c) {
      tmp.set(c, -u.i(c), u.r(c));
    }
    return tmp;
  }


  inline Mat_SU_N operator+(const Mat_SU_N& m1, const Mat_SU_N& m2)
  {
    return Mat_SU_N(m1) += m2;
  }


  inline Mat_SU_N operator-(const Mat_SU_N& m1, const Mat_SU_N& m2)
  {
    return Mat_SU_N(m1) -= m2;
  }


  inline Mat_SU_N operator*(const Mat_SU_N& m1, const Mat_SU_N& m2)
  {
    return Mat_SU_N(m1) *= m2;
  }


  inline Mat_SU_N operator*(const Mat_SU_N& m, const double a)
  {
    return Mat_SU_N(m) *= a;
  }


  inline Mat_SU_N operator*(const double a, const Mat_SU_N& m)
  {
    return Mat_SU_N(m) *= a;
  }


  inline Mat_SU_N reunit(const Mat_SU_N& m)
  {
    Mat_SU_N tmp = m;

    tmp.reunit();
    return tmp;
  }


  inline Mat_SU_N mat_exp(const double alpha, const Mat_SU_N& iv, const Mat_SU_N& u, const int Nprec)
  {
    // w = exp(alpha * iv) * u
    //   = (u + alpha * iv * (u + alpha/2 * iv * ( ... (u + alpha/n * iv * u) ... )

    Mat_SU_N p(u), tmp(u);

    tmp = u;

    for (int iprec = 0; iprec < Nprec; ++iprec) {
      double exp_factor = alpha / (Nprec - iprec);

      p    = iv * tmp;
      p   *= exp_factor;
      tmp  = p;
      tmp += u;  // w' = u + alpha/(n-k) alpha iv ( w )
    }

    return tmp.reunit();
  }
} // namespace SU_N
#endif
