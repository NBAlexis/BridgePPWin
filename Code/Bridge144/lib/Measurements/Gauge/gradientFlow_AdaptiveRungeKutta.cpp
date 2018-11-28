#include "BridgeLib_Private.h"

/*!
        @file    $Id:: gradientFlow_AdaptiveRungeKutta.cpp #$

        @brief

        @author  Yusuke Namekawa (namekawa)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 00:49:46 #$

        @version $LastChangedRevision: 1561 $
*/

#include "gradientFlow_AdaptiveRungeKutta.h"

const std::string GradientFlow_AdaptiveRungeKutta::class_name = "GradientFlow_AdaptiveRungeKutta";

//====================================================================
void GradientFlow_AdaptiveRungeKutta::flow(double& t, double& Estep, Field_G& U)
{
  while (true)
  {
    copy(U_rough, U);
    copy(U_fine, U);

    //- rough step
    double t2     = t;
    double estep2 = 2.0 * Estep;

    m_RK->flow(t2, estep2, U_rough);

    //- two original steps
    double t1     = t;
    double estep1 = Estep;

    m_RK->flow(t1, estep1, U_fine);
    m_RK->flow(t1, estep1, U_fine);

    double diff = max_diff_U(U_rough, U_fine);

    vout.general(m_vl, "  Estep,diff,m_tolerance = %e %e %e\n", Estep, diff, m_tolerance);

    if (diff < m_tolerance) {
      copy(U, U_fine);
      t += estep2;

      // extend stepsize
      Estep *= m_safety * pow(m_tolerance / diff, 1.0 / (Norder_RK() + 1));

      break;  // exit loop
    } else {
      // shrink stepsize
      Estep *= m_safety * pow(m_tolerance / diff, 1.0 / (Norder_RK() + 1));

      if (Estep < m_tolerance) {
        vout.crucial(m_vl, "Error at %s: too small Estep = %e\n", class_name.c_str(), Estep);
        exit(EXIT_FAILURE);
      }

      // proceed to next trial
    }
  }
}


//====================================================================
double GradientFlow_AdaptiveRungeKutta::max_diff_U(const Field_G& U1, const Field_G& U0) const
{
  int Nvol = U1.nvol();
  int Nex  = U0.nex();
  int Nc   = CommonParameters::Nc();

  Mat_SU_N u0(Nc), u1(Nc), u2(Nc);

  double max_norm = 0.0;

  for (int ex = 0; ex < Nex; ++ex) {
    for (int site = 0; site < Nvol; ++site) {
      U1.mat(u0, site, ex);
      U0.mat(u1, site, ex);

      u2 = u1 - u0;

      double norm = sqrt(u2.norm2()) / Nc;

      if (norm > max_norm) max_norm = norm;
    }
  }

  max_norm = Communicator::reduce_max(max_norm);

  return max_norm;
}


//====================================================================
//============================================================END=====
