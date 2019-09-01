/*!
        @file    field.h

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/


#ifndef FIELD_INCLUDED
#define FIELD_INCLUDED

#include <iostream>
#include <valarray>
#include <string>
#include <assert.h>

#include "Parameters/commonParameters.h"
#include "Communicator/communicator.h"

#include "bridge_complex.h"
#include "bridge_defs.h"

//! Container of Field-element_type object.

/*!
   This class BAPI defines field-element_type quantity which has three
   size parameters, Nin: on-site degree of freedom,
   Nvol: site d.o.f, Nex: extra d.o.f.
   The detailed structure of these degrees of freedom is not
   defined in this class BAPI but in subclasses.
   Expression template was implemented by J.Noaki.
                                   [28 Dec 2011 H.Matsufuru]
 */
using Element_type::element_type;

class BAPI Field
{
 public:
  //  enum element_type { REAL = 1, COMPLEX = 2 };
  //typedef Element_type::element_type   element_type;
  typedef double               real_t; // [23 Apr 2018 H.M. added]

 protected:
  int m_Nin;           // internal d.o.f.
  int m_Nvol;          // lattice volume
  int m_Nex;           // external d.o.f.
  // Element_type::element_type  m_element_type;
  element_type m_element_type;
  // total degree of freedom is Nin * Nsite * Nex.

  std::valarray<double> field;

  inline
  int myindex(const int jin, const int site, const int jex) const
  {
    return jin + m_Nin * (site + m_Nvol * jex);
  }

  Bridge::VerboseLevel m_vl;

 public:

  Field() :
    m_Nin(0), m_Nvol(0), m_Nex(0), m_element_type(Element_type::COMPLEX),
    field(0),
    m_vl(CommonParameters::Vlevel())
  {
    check();
  }

  Field(const int Nin, const int Nvol, const int Nex,
        const element_type cmpl = Element_type::COMPLEX) :
    m_Nin(Nin), m_Nvol(Nvol), m_Nex(Nex), m_element_type(cmpl),
    field(ntot()),
    m_vl(CommonParameters::Vlevel())
  {
    check();
  }

  Field clone() const
  {
    return Field(m_Nin, m_Nvol, m_Nex, m_element_type);
  }

  void reset(const int Nin, const int Nvol, const int Nex,
             const element_type cmpl = Element_type::COMPLEX)
  {
    if ((m_Nin == Nin) &&
        (m_Nvol == Nvol) &&
        (m_Nex == Nex) &&
        (m_element_type == cmpl)) return;

    m_Nin          = Nin;
    m_Nvol         = Nvol;
    m_Nex          = Nex;
    m_element_type = cmpl;

    field.resize(ntot());
  }

  // assignment
  //Field& operator=(const double a) { field = a; return *this; }
  Field& operator=(const Field& v)
  {
    assert(m_Nin == v.nin());
    assert(m_Nvol == v.nvol());
    assert(m_Nex == v.nex());
    copy(*this, v);
    return *this;
  }

 private:
  void check();

 public:
  int nin() const { return m_Nin; }
  int nvol() const { return m_Nvol; }
  int nex() const { return m_Nex; }
  element_type field_element_type() const { return m_element_type; }

  int ntot() const { return m_Nin * m_Nvol * m_Nex; }
  int size() const { return ntot(); }

  //! checking size parameters. [23 May 2016 H.Matsufuru]
  bool check_size(const int nin, const int nvol, const int nex) const
  {
    bool chk = true;

    if ((m_Nin != nin) || (m_Nvol != nvol) || (m_Nex != nex)) chk = false;
    return chk;
  }

  double cmp(const int jin, const int site, const int jex) const
  {
    return field[myindex(jin, site, jex)];
  }

  double cmp(const int i) const
  {
    return field[i];
  }

  const double *ptr(const int jin, const int site, const int jex) const
  {
    // when using c++03, const_cast is necessary.
    return &(const_cast<std::valarray<double>&>(field)[myindex(jin, site, jex)]);
  }

  double *ptr(const int jin, const int site, const int jex)
  {
    return &field[myindex(jin, site, jex)];
  }

  const double *ptr(const int i) const
  {
    // when using c++03, const_cast is necessary.
    return &(const_cast<std::valarray<double>&>(field)[i]);
  }

  double *ptr(const int i)
  {
    return &field[i];
  }

  void set(const int jin, const int site, const int jex, double v)
  {
    field[myindex(jin, site, jex)] = v;
  }

  void set(const int i, double v)
  {
    field[i] = v;
  }

  void set(double a);

  void add(const int jin, const int site, const int jex, double v)
  {
    field[myindex(jin, site, jex)] += v;
  }

  void add(const int i, double v)
  {
    field[i] += v;
  }

  void setpart_ex(int ex, const Field& w, int exw)
  {
    assert(ex < nex());
    assert(exw < w.nex());
    return copy(*this, ex, w, exw);
  }

  void addpart_ex(int ex, const Field& w, int exw)
  {
    assert(ex < nex());
    assert(exw < w.nex());
    return axpy(*this, ex, 1.0, w, exw);
  }

  void addpart_ex(int ex, const Field& w, int exw, double prf)
  {
    assert(ex < nex());
    assert(exw < w.nex());
    return axpy(*this, ex, prf, w, exw);
  }

  // BLAS-like routine (friend declaration)

  double norm2() const;

  double norm() const { return sqrt(norm2()); }

  friend BAPI
  double dot(const Field& y, const Field& x);

  friend BAPI
  double dot(const Field& y, const int exy, const Field& x, const int exx);

  friend BAPI
  dcomplex dotc(const Field& y, const Field& x);

  friend BAPI
  dcomplex dotc(const Field& y, const int exy, const Field& x, const int exx);

  friend BAPI
  void axpy(Field& y, const double a, const Field& x);

  friend BAPI
  void axpy(Field& y, const int exy, const double a, const Field& x, const int exx);

  friend BAPI
  void axpy(Field& y, const dcomplex a, const Field& x);

  friend BAPI
  void axpy(Field& y, const int exy, const dcomplex a, const Field& x, const int exx);

  friend BAPI
  void scal(Field& x, const double a);

  friend BAPI
  void scal(Field& x, const int exx, const double a);

  friend BAPI
  void scal(Field& x, const dcomplex a);

  friend BAPI
  void scal(Field& x, const int exx, const dcomplex a);

  friend BAPI
  void copy(Field& y, const Field& x);

  friend BAPI
  void copy(Field& y, const int exy, const Field& x, const int exx);

  friend BAPI
  void aypx(const double a, Field& y, const Field& x);

  friend BAPI
  void aypx(const dcomplex a, Field& y, const Field& x);


  //! \brief determines the statistics of the field.
  //! average, maximum value, and deviation is determined
  //! over global lattice. On-site degree of freedom is
  //! sumed over in quadrature, not averaged.
  //! This function works only on single node.
  void stat(double& Fave, double& Fmax, double& Fdev) const;
};

//----------------------------------------------------------------
// BLAS-like routines defined outside the Field class BAPI.

//! dot(y,x) := y^T x
//!   N.B. treat as real vectors.
double BAPI dot(const Field& y, const Field& x);

//! dot(y[j], x[k]) := y[j]^T x[k]
double BAPI dot(const Field& y, const int exy, const Field& x, const int exx);

//! dotc(y, x) := y^dag x
//!   x, y may be real or complex vectors. (must be of same element element_type)
dcomplex BAPI dotc(const Field& y, const Field& x);

//! dotc(y, j, x, k) := y[j]^dag x[k]
dcomplex BAPI dotc(const Field& y, const int exy, const Field& x, const int nexx);

//! axpy(y, a, x):  y := a * x + y
void BAPI axpy(Field& y, const double a, const Field& x);

//! axpy(y, j, a, x, k):  y[j] := a * x[k] + y[j]
void BAPI axpy(Field& y, const int exy, const double a, const Field& x, const int exx);

//! axpy(y, a, x):  y := a * x + y
//!   fails if x or y is real vector
void BAPI axpy(Field& y, const dcomplex a, const Field& x);

//! axpy(y, j, a, x, k):  y[j] := a * x[k] + y[j]
void BAPI axpy(Field& y, const int exy, const dcomplex a, const Field& x, const int exx);

//! scal(x, a): x = a * x
void BAPI scal(Field& x, const double a);

//! scal(x, k, a): x[k] = a * x[k]
void BAPI scal(Field& x, const int exx, const double a);

//! scal(x, a): x = a * x
//!   fails if x is a real vector
void BAPI scal(Field& x, const dcomplex a);

//! scal(x, k, a): x[k] = a * x[k]
void BAPI scal(Field& x, const int exx, const dcomplex a);

//! copy(y, x): y = x
void BAPI copy(Field& y, const Field& x);

//! copy(y, j, x, k): y[j] = x[k]
void BAPI copy(Field& y, const int exy, const Field& x, const int exx);

//! aypx(y, a, x):  y := a * y + x
void BAPI aypx(const double a, Field& y, const Field& x);

//! aypx(y, a, x):  y := a * y + x
//!   fails if x or y is real vector
void BAPI aypx(const dcomplex a, Field& y, const Field& x);

//----------------------------------------------------------------

inline int exchange(int count, Field *recv_buf, Field *send_buf, int idir, int ipm, int tag)
{
  return Communicator::exchange(count, recv_buf->ptr(0), send_buf->ptr(0), idir, ipm, tag);
}


inline int send_1to1(int count, Field *recv_buf, Field *send_buf, int p_to, int p_from, int tag)
{
  return Communicator::send_1to1(count, recv_buf->ptr(0), send_buf->ptr(0), p_to, p_from, tag);
}


//----------------------------------------------------------------
// utility: report statistics.
void BAPI report_field_stat(const Bridge::VerboseLevel vl, const std::string& msg, const Field& f);

//----------------------------------------------------------------
#endif
