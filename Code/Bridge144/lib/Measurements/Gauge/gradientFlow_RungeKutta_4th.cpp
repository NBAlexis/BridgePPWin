#include "BridgeLib_Private.h"

/*!
        @file    $Id:: gradientFlow_RungeKutta_4th.cpp #$

        @brief

        @author  Yusuke Namekawa (namekawa)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 00:49:46 #$

        @version $LastChangedRevision: 1561 $
*/

#include "gradientFlow_RungeKutta_4th.h"

const std::string GradientFlow_RungeKutta_4th::class_name = "GradientFlow_RungeKutta_4th";

//====================================================================
void GradientFlow_RungeKutta_4th::flow(double& t, double& Estep, Field_G& U)
{
  //- aliases
  Field_G& w0 = U;
  Field_G& w5 = U;

  Field_G& w1 = m_w1;
  Field_G& w2 = m_w2;
  Field_G& w3 = m_w3;
  Field_G& w4 = m_w4;

  Field_G& z0 = m_z0;
  Field_G& z1 = m_z1;
  Field_G& z2 = m_z2;
  Field_G& z3 = m_z3;

  Field_G& zt = m_zt;

  //- step 0
  // calculate gradient of m_action Z_0 (SA)
  m_action->force(z0, w0);

  copy(zt, z0);
  scal(zt, 0.5);

  // W_1=e^{Z_0/2}*U (SA)
  mult_exp_Field_G(w1, Estep, zt, w0, m_Nprec);

  //- step 1
  // Z_1
  m_action->force(z1, w1);

  copy(zt, z1);
  scal(zt, 0.5);

  // W_2=e^{(1/2)*Z_1}*W_0
  mult_exp_Field_G(w2, Estep, zt, w0, m_Nprec);

  //- step 2
  // Z_2
  m_action->force(z2, w2);

  // Z_2 - (1/2)*Z_0
  copy(zt, z2);
  axpy(zt, -0.5, z0);

  // W_3=e^{Z_2 - (1/2)*Z_0}*W_1  NB. not W_0
  mult_exp_Field_G(w3, Estep, zt, w1, m_Nprec);

  //- step 3
  // Z_3
  m_action->force(z3, w3);

  // 3*Z_0+2*Z_1+2*Z_2-Z_3
  copy(zt, z0);
  scal(zt, 0.25);
  axpy(zt, 1.0 / 6.0, z1);
  axpy(zt, 1.0 / 6.0, z2);
  axpy(zt, -1.0 / 12.0, z3);

  // W_4=e^{(1/12)*(3*Z_0+2*Z_1+2*Z_2-Z_3)}*W_0
  mult_exp_Field_G(w4, Estep, zt, w0, m_Nprec);

  //- step 4 (extra step)
  // 1/12*(-Z_0+2*Z_1+2*Z_2+3*Z_3)
  copy(zt, z0);
  scal(zt, -1.0 / 12.0);
  axpy(zt, 1.0 / 6.0, z1);
  axpy(zt, 1.0 / 6.0, z2);
  axpy(zt, 0.25, z3);

  // V_t=e^{(1/12)*(-Z_0+2*Z_1+2*Z_2+3*Z_3)}*W_4  NB. not W_0
  mult_exp_Field_G(w5, Estep, zt, w4, m_Nprec);

  t += Estep;
}


//====================================================================
//============================================================END=====
