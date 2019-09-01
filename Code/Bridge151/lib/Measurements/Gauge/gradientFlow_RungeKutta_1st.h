/*!
        @file    gradientFlow_RungeKutta_1st.h

        @brief

        @author  Yusuke Namekawa (namekawa)

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/

#ifndef GRADIENTFLOW_RUNGEKUTTA_1ST_INCLUDED
#define GRADIENTFLOW_RUNGEKUTTA_1ST_INCLUDED

#include "gradientFlow_RungeKutta.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//! GradientFlow_RungeKutta_1st construction.

/*!
    This class BAPI implements 1st order Runge-Kutta method
    for GradientFlow. No commutator appears at 1st order RK.
                               [10 Oct 2014 Y.Namekawa]
 */

class BAPI GradientFlow_RungeKutta_1st : public GradientFlow_RungeKutta
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
  Field_G m_z0;

 public:
  GradientFlow_RungeKutta_1st(Action *action, const int Nprec, const Bridge::VerboseLevel vl)
    : GradientFlow_RungeKutta(action, Nprec, vl)
  {
    m_action = action;
    m_Nprec  = Nprec;
    m_Ndim   = CommonParameters::Ndim();
    m_Nvol   = CommonParameters::Nvol();

    m_z0.reset(m_Nvol, m_Ndim);
  }

  ~GradientFlow_RungeKutta_1st() {}

  void flow(double& t, double& Estep, Field_G& U);

  int Norder_RK() const { return 1; }
};
#endif
