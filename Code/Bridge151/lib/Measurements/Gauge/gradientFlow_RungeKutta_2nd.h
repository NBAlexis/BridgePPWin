/*!
        @file    gradientFlow_RungeKutta_2nd.h

        @brief

        @author  Yusuke Namekawa (namekawa)

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/

#ifndef GRADIENTFLOW_RUNGEKUTTA_2ND_INCLUDED
#define GRADIENTFLOW_RUNGEKUTTA_2ND_INCLUDED

#include "gradientFlow_RungeKutta.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//! GradientFlow_RungeKutta_2nd construction.

/*!
    This class BAPI implements 2nd order Runge-Kutta method
    for GradientFlow. No commutator appears at 2nd order RK.
                               [10 Oct 2014 Y.Namekawa]
 */

class BAPI GradientFlow_RungeKutta_2nd : public GradientFlow_RungeKutta
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
  Field_G m_w1;
  Field_G m_z0, m_z1;

 public:
  GradientFlow_RungeKutta_2nd(Action *action, const int Nprec, const Bridge::VerboseLevel vl)
    : GradientFlow_RungeKutta(action, Nprec, vl)
  {
    m_action = action;
    m_Nprec  = Nprec;
    m_Ndim   = CommonParameters::Ndim();
    m_Nvol   = CommonParameters::Nvol();

    m_w1.reset(m_Nvol, m_Ndim);
    m_z0.reset(m_Nvol, m_Ndim);
    m_z1.reset(m_Nvol, m_Ndim);
  }

  ~GradientFlow_RungeKutta_2nd() {}

  void flow(double& t, double& Estep, Field_G& U);

  int Norder_RK() const { return 2; }
};
#endif
