/*!
        @file    gradientFlow_RungeKutta_2nd.cpp

        @brief

        @author  Yusuke Namekawa (namekawa)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/
#include "BridgeLib_Private.h"
#include "gradientFlow_RungeKutta_2nd.h"

const std::string GradientFlow_RungeKutta_2nd::class_name = "GradientFlow_RungeKutta_2nd";

//====================================================================
void GradientFlow_RungeKutta_2nd::flow(double& t, double& Estep, Field_G& U)
{
  //- aliases
  Field_G& w0 = U;
  Field_G& w2 = U;

  Field_G& w1 = m_w1;

  Field_G& z0 = m_z0;
  Field_G& z1 = m_z1;

  //- step 0
  // calculate gradient of m_action Z_0 (SA)
  m_action->force(z0, w0);

  // W_1=e^{c2*Z_0}*U
  mult_exp_Field_G(w1, 0.5 * Estep, z0, w0, m_Nprec);

  //- step 1
  // Z_1
  m_action->force(z1, w1);

  // V_out=e^{Z_1}*W_0
  mult_exp_Field_G(w2, Estep, z1, w0, m_Nprec);

  t += Estep;
}


//====================================================================
//============================================================END=====
