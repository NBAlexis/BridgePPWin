#include "BridgeLib_Private.h"

/*!
        @file    $Id:: gradientFlow_RungeKutta_1st.cpp #$

        @brief

        @author  Yusuke Namekawa (namekawa)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 00:49:46 #$

        @version $LastChangedRevision: 1561 $
*/

#include "gradientFlow_RungeKutta_1st.h"

const std::string GradientFlow_RungeKutta_1st::class_name = "GradientFlow_RungeKutta_1st";

//====================================================================
void GradientFlow_RungeKutta_1st::flow(double& t, double& Estep, Field_G& U)
{
  //- aliases
  Field_G& w0 = U;
  Field_G& w1 = U;

  Field_G& z0 = m_z0;

  //- step 0
  // calculate gradient of m_action Z_0 (SA)
  m_action->force(z0, w0);

  // W_1=e^{Z_0}*U
  mult_exp_Field_G(w1, Estep, z0, w0, m_Nprec);  // assume it is ok for w = u case.

  t += Estep;
}


//====================================================================
//============================================================END=====
