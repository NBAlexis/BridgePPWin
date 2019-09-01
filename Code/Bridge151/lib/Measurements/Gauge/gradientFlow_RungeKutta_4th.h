/*!
        @file    gradientFlow_RungeKutta_4th.h

        @brief

        @author  Yusuke Namekawa (namekawa)

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/

#ifndef GRADIENTFLOW_RUNGEKUTTA_4TH_INCLUDED
#define GRADIENTFLOW_RUNGEKUTTA_4TH_INCLUDED

#include "gradientFlow_RungeKutta.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//! GradientFlow_RungeKutta_4th construction.

/*!
    This class BAPI implements 4th order Runge-Kutta
    in commutator-free method for GradientFlow.
    See E.Celledoni et al. FGCS 19, 341 (2003).
                               [10 Oct 2014 Y.Namekawa]
 */

class BAPI GradientFlow_RungeKutta_4th : public GradientFlow_RungeKutta
{
 public:
  static const std::string class_name;

 protected:
  Bridge::VerboseLevel m_vl;

 private:
  Action *m_action;
  int m_Nprec;
  int m_Ndim, m_Nvol;

  //- working area
  Field_G m_w1, m_w2, m_w3, m_w4;
  Field_G m_z0, m_z1, m_z2, m_z3, m_z4;
  Field_G m_zt;

 public:
  GradientFlow_RungeKutta_4th(Action *action, const int Nprec, const Bridge::VerboseLevel vl)
    : GradientFlow_RungeKutta(action, Nprec, vl)
  {
    m_action = action;
    m_Nprec  = Nprec;
    m_Ndim   = CommonParameters::Ndim();
    m_Nvol   = CommonParameters::Nvol();

    m_w1.reset(m_Nvol, m_Ndim);
    m_w2.reset(m_Nvol, m_Ndim);
    m_w3.reset(m_Nvol, m_Ndim);
    m_w4.reset(m_Nvol, m_Ndim);

    m_z0.reset(m_Nvol, m_Ndim);
    m_z1.reset(m_Nvol, m_Ndim);
    m_z2.reset(m_Nvol, m_Ndim);
    m_z3.reset(m_Nvol, m_Ndim);
    m_z4.reset(m_Nvol, m_Ndim);

    m_zt.reset(m_Nvol, m_Ndim);
  }

  ~GradientFlow_RungeKutta_4th() {}

  void flow(double& t, double& Estep, Field_G& U);

  int Norder_RK() const { return 4; }
};
#endif
