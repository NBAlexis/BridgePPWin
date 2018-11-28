#include "BridgeLib_Private.h"

/*!
        @file    $Id:: mat_SU_N.cpp #$

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 00:49:46 #$

        @version $LastChangedRevision: 1561 $
*/

#include "mat_SU_N.h"
using namespace SU_N;

//==========================================================
Mat_SU_N& Mat_SU_N::reunit_SU3()
{
  // This implementation is applicable only to SU(3).
  //                             H.Matsufuru [5 Feb 2012]

  assert(m_Nc == 3);

  Mat_SU_N rhs(m_Nc);
  for (int i = 0; i < 2 * m_Nc * m_Nc; ++i) {
    rhs.va[i] = va[i];
  }

  double sn1 = 1.0 / (rhs.va[0] * rhs.va[0] + rhs.va[1] * rhs.va[1]
                      + rhs.va[2] * rhs.va[2] + rhs.va[3] * rhs.va[3]
                      + rhs.va[4] * rhs.va[4] + rhs.va[5] * rhs.va[5]);
  sn1 = sqrt(sn1);

  va[0] = rhs.va[0] * sn1;
  va[1] = rhs.va[1] * sn1;
  va[2] = rhs.va[2] * sn1;
  va[3] = rhs.va[3] * sn1;
  va[4] = rhs.va[4] * sn1;
  va[5] = rhs.va[5] * sn1;

  double sp1r = va[0] * rhs.va[6] + va[1] * rhs.va[7]
                + va[2] * rhs.va[8] + va[3] * rhs.va[9]
                + va[4] * rhs.va[10] + va[5] * rhs.va[11];
  double sp1i = va[0] * rhs.va[7] - va[1] * rhs.va[6]
                + va[2] * rhs.va[9] - va[3] * rhs.va[8]
                + va[4] * rhs.va[11] - va[5] * rhs.va[10];

  va[6]  = rhs.va[6] - sp1r * va[0] + sp1i * va[1];
  va[7]  = rhs.va[7] - sp1r * va[1] - sp1i * va[0];
  va[8]  = rhs.va[8] - sp1r * va[2] + sp1i * va[3];
  va[9]  = rhs.va[9] - sp1r * va[3] - sp1i * va[2];
  va[10] = rhs.va[10] - sp1r * va[4] + sp1i * va[5];
  va[11] = rhs.va[11] - sp1r * va[5] - sp1i * va[4];

  double sn2 = 1.0 / (va[6] * va[6] + va[7] * va[7]
                      + va[8] * va[8] + va[9] * va[9]
                      + va[10] * va[10] + va[11] * va[11]);
  sn2 = sqrt(sn2);

  va[6]  = va[6] * sn2;
  va[7]  = va[7] * sn2;
  va[8]  = va[8] * sn2;
  va[9]  = va[9] * sn2;
  va[10] = va[10] * sn2;
  va[11] = va[11] * sn2;

  va[12] = va[2] * va[10] - va[3] * va[11]
           - va[4] * va[8] + va[5] * va[9];

  va[13] = -va[2] * va[11] - va[3] * va[10]
           + va[4] * va[9] + va[5] * va[8];

  va[14] = va[4] * va[6] - va[5] * va[7]
           - va[0] * va[10] + va[1] * va[11];

  va[15] = -va[4] * va[7] - va[5] * va[6]
           + va[0] * va[11] + va[1] * va[10];

  va[16] = va[0] * va[8] - va[1] * va[9]
           - va[2] * va[6] + va[3] * va[7];

  va[17] = -va[0] * va[9] - va[1] * va[8]
           + va[2] * va[7] + va[3] * va[6];

  return *this;
}


//==========================================================
Mat_SU_N& Mat_SU_N::reunit_SU2()
{
  // This implementation is applicable only to SU(2).
  //                             H.Matsufuru [8 Mar 2012]

  assert(m_Nc == 2);

  double sn1 = 1.0 / (va[0] * va[0] + va[1] * va[1]
                      + va[2] * va[2] + va[3] * va[3]);
  sn1 = sqrt(sn1);

  va[0] *= sn1;
  va[1] *= sn1;
  va[2] *= sn1;
  va[3] *= sn1;

  va[4] = -va[2];
  va[5] = va[3];
  va[6] = va[0];
  va[7] = -va[1];

  return *this;
}


//==========================================================
Mat_SU_N& Mat_SU_N::reunit_general()
{
  assert(m_Nc > 3);

  Decompose_QR_Cmplx qr(m_Nc);
  qr.set_matrix(&va[0]);
  qr.get_Qu(&va[0]);

  Eigen_QR_Cmplx        eigen_qr(m_Nc);
  std::valarray<double> v       = eigen_qr.solve(&va[0]);
  double                det_arg = 0;
  for (int i = 0; i < m_Nc; ++i) {
    det_arg += std::atan2(v[2 * i + 1], v[2 * i]);
  }
  det_arg /= m_Nc;
  double inv_re = std::cos(det_arg);
  double inv_im = -std::sin(det_arg);

  for (int i = 0; i < m_Nc * m_Nc; ++i) {
    double re = va[2 * i];
    double im = va[2 * i + 1];

    va[2 * i]     = re * inv_re - im * inv_im;
    va[2 * i + 1] = re * inv_im + im * inv_re;
  }

  return *this;
}


//==========================================================
Mat_SU_N& Mat_SU_N::set_random_SU3(RandomNumbers *rand)
{
  // temporary implementation: only applicable to SU(3).
  //                            H.Matsufuru [5 Feb 2012]

  assert(m_Nc == 3);

  static const double PI  = 4.0 * atan(1.0);
  double              PI2 = 2.0 * PI;

  for (int j = 0; j < m_Nc; ++j) {
    int    j2    = j * 2 * m_Nc;
    double rand1 = rand->get();
    double rand2 = rand->get();
    double rand3 = rand->get();
    double rand4 = rand->get();
    double rand5 = rand->get();

    double c1   = 1.0 - 2.0 * rand1;
    double s1   = sqrt(1.0 - c1 * c1);
    double v1_2 = s1 * cos(PI2 * rand2);
    double v1_3 = s1 * sin(PI2 * rand2);

    va[j2 + 0] = c1 * cos(PI2 * rand3);
    va[j2 + 1] = c1 * sin(PI2 * rand3);
    va[j2 + 2] = v1_2 * cos(PI2 * rand4);
    va[j2 + 3] = v1_2 * sin(PI2 * rand4);
    va[j2 + 4] = v1_3 * cos(PI2 * rand5);
    va[j2 + 5] = v1_3 * sin(PI2 * rand5);
  }

  reunit();

  return *this;
}


//==========================================================
Mat_SU_N& Mat_SU_N::set_random_SU2(RandomNumbers *rand)
{
  assert(m_Nc == 2);

  printf("Error at Mat_SU_N: set_random for SU(2) is not supported.\n");
  exit(EXIT_FAILURE);
}


//==========================================================
Mat_SU_N& Mat_SU_N::set_random_general(RandomNumbers *rand)
{
  // implemented by S.Ueda.

  assert(m_Nc > 3);

  for (int j = 0; j < m_Nc * m_Nc; ++j) {
    rand->gauss(va[2 * j], va[2 * j + 1]);
  }

  reunit();

  return *this;
}


//==========================================================
//==================================================END=====
