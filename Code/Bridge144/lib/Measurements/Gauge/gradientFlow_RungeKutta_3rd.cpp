#include "BridgeLib_Private.h"

/*!
        @file    $Id:: gradientFlow_RungeKutta_3rd.cpp #$

        @brief

        @author  Yusuke Namekawa (namekawa)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 00:49:46 #$

        @version $LastChangedRevision: 1561 $
*/

#include "gradientFlow_RungeKutta_3rd.h"

const std::string GradientFlow_RungeKutta_3rd::class_name = "GradientFlow_RungeKutta_3rd";

//====================================================================
void GradientFlow_RungeKutta_3rd::flow(double& t, double& Estep, Field_G& U)
{
  const double a0 = 0.25;

  const double b0 = -17.0 / 36.0;
  const double b1 = 8.0 / 9.0;

  const double c0 = 17.0 / 36.0;
  const double c1 = -8.0 / 9.0;
  const double c2 = 0.75;

  //- aliases
  Field_G& w0 = U;
  Field_G& w3 = U;

  Field_G& w1 = m_w1;
  Field_G& w2 = m_w2;

  Field_G& z0 = m_z0;
  Field_G& z1 = m_z1;
  Field_G& z2 = m_z2;

  Field_G& zt = m_zt;

  //- step 0
  // calculate gradient of m_action Z_0 (SA)
  m_action->force(z0, w0);

  // zt = a0 * z0 = 1/4 z0
  copy(zt, z0);
  scal(zt, a0);

  // W_1=e^{Z_0/4}*U (SA)
  mult_exp_Field_G(w1, Estep, zt, w0, m_Nprec);

  //- step 1
  // Z_1 (SA)
  m_action->force(z1, w1);

  // zt = b1 * z1 + b0 * z0 = 8/9 z1 - 17/36 z0
  copy(zt, z1);
  scal(zt, b1);
  axpy(zt, b0, z0);

  // W_2=e^{8*Z_1/9-17*Z_0/36}*W_1 (SA)
  mult_exp_Field_G(w2, Estep, zt, w1, m_Nprec);

  //- step 2
  // Z_2 (SA)
  m_action->force(z2, w2);

  // zt = c2 * z2 + c1 * z1 * c0 * z0
  //    = 3/4 z2 - 8/9 z1 + 17/36 z0
  //    = 3/4 z2 - zt_prev
  copy(zt, z2);
  scal(zt, c2);
  axpy(zt, c1, z1);
  axpy(zt, c0, z0);

  // V_out=e^{3*Z_2/4-8*Z_1/9+17*Z_0/36}*W_2 (SA)
  mult_exp_Field_G(w3, Estep, zt, w2, m_Nprec);

  t += Estep;
}


//====================================================================
//============================================================END=====
