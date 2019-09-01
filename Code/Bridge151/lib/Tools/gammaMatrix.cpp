/*!
        @file    gammaMatrix.cpp

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/
#include "BridgeLib_Private.h"
#include "gammaMatrix.h"

#include <cassert>

#include "Communicator/communicator.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

const std::string GammaMatrix::class_name = "GammaMatrix";

//====================================================================
void GammaMatrix::set(int row, int index, icomplex val)
{
  assert(row < m_Nd);
  assert(index < m_Nd);

  m_gmindex[row] = index; // m_gmindex[a] = b(=index) if \gamma_{ab} \not=0 (SA)
  m_gmval_i[row] = val;   //  m_gmval_i[a] = \gamma_{ab}(=val) in complex<int> (SA)

  set_values(row);
}


//====================================================================
GammaMatrix GammaMatrix::mult(GammaMatrix gm1) const
{
  // \gamma2_{ab}=\gamma_{ac} * \gamma1_{cb} (SA)

  GammaMatrix gm2;

  for (int row = 0; row < m_Nd; ++row) {
    gm2.m_gmindex[row] = gm1.m_gmindex[m_gmindex[row]];
    gm2.m_gmval_i[row] = m_gmval_i[row]
                         * gm1.m_gmval_i[m_gmindex[row]];
    // \gamma_{ac} * \gamma1_{cb}, c=m_gmindex[a], (a=row) (SA)

    gm2.set_values(row);
  }

  return gm2;
}


//====================================================================
GammaMatrix GammaMatrix::mult_i(GammaMatrix gm1) const
{
  // \gamma2_{ab}=i * \gamma_{ac} * \gamma1_{cb} (SA)

  GammaMatrix gm2;

  for (int row = 0; row < m_Nd; ++row) {
    gm2.m_gmindex[row] = gm1.m_gmindex[m_gmindex[row]];
    gm2.m_gmval_i[row] = icomplex(0, 1) * m_gmval_i[row]
                         * gm1.m_gmval_i[m_gmindex[row]];
    gm2.set_values(row);
  }

  return gm2;
}


//====================================================================
void GammaMatrix::print()
{
  for (int row = 0; row < m_Nd; ++row) {
    if (m_gmindex[row] == 0) {
      vout.general(m_vl, "    (%2d,%2d)\n",
                   real(m_gmval_i[row]), imag(m_gmval_i[row]));
    } else if (m_gmindex[row] == 1) {
      vout.general(m_vl, "             (%2d,%2d)\n",
                   real(m_gmval_i[row]), imag(m_gmval_i[row]));
    } else if (m_gmindex[row] == 2) {
      vout.general(m_vl, "                      (%2d,%2d)\n",
                   real(m_gmval_i[row]), imag(m_gmval_i[row]));
    } else if (m_gmindex[row] == 3) {
      vout.general(m_vl, "                               (%2d,%2d)\n",
                   real(m_gmval_i[row]), imag(m_gmval_i[row]));
    } else {
      vout.crucial(m_vl, "Error at %s: gamma matrix not defined.\n", class_name.c_str());
      exit(EXIT_FAILURE);
    }
  }
}


//====================================================================
void GammaMatrix::set_values(int row)
{
  m_gmvalue[row] = cmplx((double)real(m_gmval_i[row]),
                         (double)imag(m_gmval_i[row]));
  // m_gmvalue[a] = \gamma_{ab} in complex<double> (SA)

  if (real(m_gmval_i[row]) != 0) { // if \gamma_{ab} is real (SA)
    m_gmindex_c[row] = 0;
    m_gmvalue_r[row] = real(m_gmvalue[row]);
    m_gmvalue_i[row] = real(m_gmvalue[row]);
  } else if (imag(m_gmval_i[row]) != 0) { // if \gamma_{ab} is complex (SA)
    m_gmindex_c[row] = 1;
    m_gmvalue_r[row] = -imag(m_gmvalue[row]);
    m_gmvalue_i[row] = imag(m_gmvalue[row]);
  } else {
    vout.crucial(m_vl, "Error at %s: wrong m_gmval_i.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }
}


//====================================================================
GammaMatrix GammaMatrix::mult(int n) const
{
  // n * \gamma_{ab} (SA)

  GammaMatrix gm2;

  for (int row = 0; row < m_Nd; ++row) {
    gm2.m_gmindex[row] = m_gmindex[row];
    gm2.m_gmval_i[row] = m_gmval_i[row] * icomplex(n, 0);
    gm2.set_values(row);
  }

  return gm2;
}


//====================================================================
//============================================================END=====
