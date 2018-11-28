/*!
        @file    $Id:: field_G.h #$

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 18:35:38 #$

        @version $LastChangedRevision: 1571 $
*/


#ifndef FIELD_G_INCLUDED
#define FIELD_G_INCLUDED

#include "field.h"
#include "Tools/randomNumbers.h"

#include "Tools/mat_SU_N.h"
using namespace SU_N;

//! SU(N) gauge field.

/*!
    This class defines SU(N) gauge field, which is used such as
    gauge configuration.
    Original version of this class was written by J.Noaki.
    H.Matsufuru added several functions and modified intefaces
    of several functionality.
                                     [28 Dec 2011 H.Matsufuru]
    (Coding history will be recovered from trac.)
    reunit is added.                 [26 May 2014 Y.Namekawa]
    unique_ptr is introduced to avoid memory leaks
                                     [21 Mar 2015 Y.Namekawa]
 */
class Field_G : public Field
{
 private:
  int m_Nc;     // number of color elements
  int m_Ndf;    // number of components as real values

 public:

  explicit
  Field_G(const int Nvol = CommonParameters::Nvol(), const int Nex = 1) :
    Field(
      2 * CommonParameters::Nc() * CommonParameters::Nc(), Nvol, Nex, COMPLEX
      ),
    m_Nc(CommonParameters::Nc()),
    m_Ndf(2 * CommonParameters::Nc() * CommonParameters::Nc())
  {
    check();
  }

  Field_G clone() const
  {
    return Field_G(nvol(), nex());
  }

  // conversion from Field type

  Field_G(const Field& x) :
    Field(x),
    m_Nc(CommonParameters::Nc()),
    m_Ndf(x.nin())
  {
    assert(m_Ndf == 2 * m_Nc * m_Nc);
    check();
  }

  // assignment
  //Field_G& operator=(const double a) { field = a; return *this; }
  Field_G& operator=(const Field_G& v) { copy(*this, v); return *this; }

  // resize field
  void reset(const int Nvol, const int Nex)
  {
    Field::reset(m_Ndf, Nvol, Nex);
  }

  int nc() const { return m_Nc; }

  // accessors
  double cmp_r(const int cc, const int site, const int mn = 0) const
  {
    return field[myindex(2 * cc, site, mn)];
  }

  double cmp_i(const int cc, const int site, const int mn = 0) const
  {
    return field[myindex(2 * cc + 1, site, mn)];
  }

  void set_r(const int cc, const int site, const int mn, const double re)
  {
    field[myindex(2 * cc, site, mn)] = re;
  }

  void set_i(const int cc, const int site, const int mn, const double im)
  {
    field[myindex(2 * cc + 1, site, mn)] = im;
  }

  void set_ri(const int cc, const int site, const int mn,
              const double re, const double im)
  {
    field[myindex(2 * cc, site, mn)]     = re;
    field[myindex(2 * cc + 1, site, mn)] = im;
  }

  Mat_SU_N mat(const int site, const int mn = 0) const
  {
    Mat_SU_N Tmp(m_Nc);

    for (int cc = 0; cc < m_Nc * m_Nc; ++cc) {
      Tmp.set(cc,
              field[myindex(2 * cc, site, mn)],
              field[myindex(2 * cc + 1, site, mn)]);
    }

    return Tmp;
  }

  Mat_SU_N mat_dag(const int site, const int mn = 0) const
  {
    Mat_SU_N Tmp(m_Nc);

    for (int cc = 0; cc < m_Nc * m_Nc; ++cc) {
      Tmp.set(cc,
              field[myindex(2 * cc, site, mn)],
              field[myindex(2 * cc + 1, site, mn)]);
    }

    return Tmp.dag();
  }

  void mat(Mat_SU_N& Tmp, const int site, const int mn = 0) const
  {
    for (int cc = 0; cc < m_Nc * m_Nc; ++cc) {
      Tmp.set(cc,
              field[myindex(2 * cc, site, mn)],
              field[myindex(2 * cc + 1, site, mn)]);
    }
  }

  void mat_dag(Mat_SU_N& Tmp, const int site, const int mn = 0) const
  {
    for (int c1 = 0; c1 < m_Nc; ++c1) {
      for (int c2 = 0; c2 < m_Nc; ++c2) {
        Tmp.set(c1 + m_Nc * c2,
                field[myindex(2 * (c2 + m_Nc * c1), site, mn)],
                -field[myindex(2 * (c2 + m_Nc * c1) + 1, site, mn)]);
      }
    }
  }

  void set_mat(const int site, const int mn, const Mat_SU_N& U)
  {
    for (int cc = 0; cc < m_Nc * m_Nc; ++cc) {
      field[myindex(2 * cc, site, mn)]     = U.r(cc);
      field[myindex(2 * cc + 1, site, mn)] = U.i(cc);
    }
  }

  void add_mat(const int site, const int mn, const Mat_SU_N& U)
  {
    for (int cc = 0; cc < m_Nc * m_Nc; ++cc) {
      field[myindex(2 * cc, site, mn)]     += U.r(cc);
      field[myindex(2 * cc + 1, site, mn)] += U.i(cc);
    }
  }

  void add_mat(const int site, const int mn, const Mat_SU_N& U, double prf)
  {
    for (int cc = 0; cc < m_Nc * m_Nc; ++cc) {
      field[myindex(2 * cc, site, mn)]     += prf * U.r(cc);
      field[myindex(2 * cc + 1, site, mn)] += prf * U.i(cc);
    }
  }

  void xI()
  {
    for (int i = 0, n = (int)field.size(); i < n; i += 2) {
      double real = field[i];
      field[i]     = -field[i + 1];
      field[i + 1] = real;
    }
  }

  void set_unit();

  void set_random(RandomNumbers *rand);
  void set_random(unique_ptr<RandomNumbers>& rand);

  void reunit();


 private:
  //! check of assumptions for performance implementation.
  void check();
};

//----------------------------------------------------------------
// function style

void mult_Field_Gnn(Field_G& W, const int ex,
                    const Field_G& U1, const int ex1,
                    const Field_G& U2, const int ex2);

void mult_Field_Gdn(Field_G& W, const int ex,
                    const Field_G& U1, const int ex1,
                    const Field_G& U2, const int ex2);

void mult_Field_Gnd(Field_G& W, const int ex,
                    const Field_G& U1, const int ex1,
                    const Field_G& U2, const int ex2);

void mult_Field_Gdd(Field_G& W, const int ex,
                    const Field_G& U1, const int ex1,
                    const Field_G& U2, const int ex2);

void multadd_Field_Gnn(Field_G& W, const int ex,
                       const Field_G& U1, const int ex1,
                       const Field_G& U2, const int ex2,
                       const double ff);

void multadd_Field_Gdn(Field_G& W, const int ex,
                       const Field_G& U1, const int ex1,
                       const Field_G& U2, const int ex2,
                       const double ff);

void multadd_Field_Gnd(Field_G& W, const int ex,
                       const Field_G& U1, const int ex1,
                       const Field_G& U2, const int ex2,
                       const double ff);

void multadd_Field_Gdd(Field_G& W, const int ex,
                       const Field_G& U1, const int ex1,
                       const Field_G& U2, const int ex2,
                       const double ff);

void at_Field_G(Field_G& W, const int ex);
void ah_Field_G(Field_G& W, const int ex);

// W = exp(alpha * iP) * U
//   = (U + alpha * iP * (U + alpha/2 * iP * ( ... (U+ alpha/n * iP * U) ...
void mult_exp_Field_G(Field_G& W, const double alpha, const Field_G& iP, const Field_G& U, const int Nprec);
#endif
