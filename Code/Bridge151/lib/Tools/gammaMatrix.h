/*!
        @file    gammaMatrix.h

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/

#ifndef GAMMAMATRIX_INCLUDED
#define GAMMAMATRIX_INCLUDED

//#include <complex>

#include "Parameters/commonParameters.h"

#include "bridge_complex.h"
//typedef std::complex<double> dcomplex;
//typedef std::complex<int> icomplex;

//! Gamma Matrix class BAPI.

/*!
    Gamma matrix is defined.
    This class BAPI just defines the action of gamma matrices.
    Specific form is not given in this class BAPI but outside.
                                        [4 Feb 2012 H.Matsufuru]
    bug fixed:
    definition of index of a product of two gamma matrices in
    mult() and mult_i() was wrong, while there had been no problem
    so far thanks to the property of the gamma matrices.
                                        [9 May 2012 H.Matsufuru]
 */

/*!
 \f$\gamma_{ab}\f$ = m_gmval_i[a], b=m_gmindex[a]  or
 \f$\gamma_{ab}\f$ = m_gmvalue[a], b=m_gmindex[a]
 [8 May 2012 S. Aoki]
 */
class BAPI GammaMatrix {
 public:
  static const std::string class_name;

 private:
  int m_Nd;
  std::vector<int> m_gmindex;        // b=m_gmindex[a] if \gamma_{ab} \not=0 (SA)
  std::vector<icomplex> m_gmval_i;   // \gamma_{ab} in complex<int> (SA)
  std::vector<dcomplex> m_gmvalue;   // \gamma_{ab} in complex<double> (SA)
  std::vector<int> m_gmindex_c;      // =1(0) if m_gmval_i[a] is pure imaginary(real)  (SA)
  std::vector<double> m_gmvalue_r;
  std::vector<double> m_gmvalue_i;

 protected:
  Bridge::VerboseLevel m_vl;
 public:
  GammaMatrix() : m_vl(CommonParameters::Vlevel())
  {
    m_Nd = CommonParameters::Nd();
    m_gmindex.resize(m_Nd);
    m_gmval_i.resize(m_Nd);
    m_gmvalue.resize(m_Nd);
    m_gmindex_c.resize(m_Nd);
    m_gmvalue_r.resize(m_Nd);
    m_gmvalue_i.resize(m_Nd);
  }

  void set(int row, int index, icomplex val_i);

  void set_values(int row);

  void print();

  GammaMatrix mult(GammaMatrix) const;

  GammaMatrix mult_i(GammaMatrix) const;

  GammaMatrix mult(int) const;

  int index(int row) const
  {
    return m_gmindex[row];
  }

  dcomplex value(int row) const
  {
    return m_gmvalue[row];
  }

  int index_c(int row) const
  {
    return m_gmindex_c[row];
  }

  double value_r(int row) const
  {
    return m_gmvalue_r[row];
  }

  double value_i(int row) const
  {
    return m_gmvalue_i[row];
  }
};
#endif
