/*!
        @file    gradientFlow_AdaptiveRungeKutta.h

        @brief

        @author  Yusuke Namekawa (namekawa)

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/

#ifndef GRADIENTFLOW_ADAPTIVERUNGEKUTTA_INCLUDED
#define GRADIENTFLOW_ADAPTIVERUNGEKUTTA_INCLUDED

#include "gradientFlow_RungeKutta.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//! GradientFlow_AdaptiveRungeKutta construction.

/*!
    This class BAPI implements Runge-Kutta method with adaptive
    step size control by step doubling for GradientFlow.
    See C.W.Gear, Englewood Cliffs, NJ: Prentice-Hall (1971).
    RK_adaptive tends to deviate from the expected_result by O(10^{-11}),
    because RK_adaptive reflects rounding errors severely.
                               [01 May 2015 Y.Namekawa]
    Adaptive stepsize algorithm is organised by composite pattern
    by Aoyama-san.
 */

class BAPI GradientFlow_AdaptiveRungeKutta : public GradientFlow_RungeKutta
{
 public:
  static const std::string class_name;

 protected:
  Bridge::VerboseLevel m_vl;

 private:
  Action *m_action;

  int m_Ndim;
  int m_Nvol;

  GradientFlow_RungeKutta *m_RK;

  double m_tolerance;
  double m_safety;

  Field_G m_U_rough, m_U_fine;

 public:
  GradientFlow_AdaptiveRungeKutta(GradientFlow_RungeKutta *RK,
                                  const double tolerance, const double safety, const Bridge::VerboseLevel vl)
    : GradientFlow_RungeKutta(0, 0, vl),
    m_vl(vl),
    m_RK(RK),
    m_tolerance(tolerance),
    m_safety(safety)
  {
    m_Ndim = CommonParameters::Ndim();
    m_Nvol = CommonParameters::Nvol();

    m_U_rough.reset(m_Nvol, m_Ndim);
    m_U_fine.reset(m_Nvol, m_Ndim);
  }

  ~GradientFlow_AdaptiveRungeKutta() { if (m_RK) delete m_RK; }

  void flow(double& t, double& Estep, Field_G& U);

  int Norder_RK() const { return m_RK ? m_RK->Norder_RK() : 0; }

 private:
  double max_diff_U(const Field_G& U1, const Field_G& U0) const;
};
#endif
