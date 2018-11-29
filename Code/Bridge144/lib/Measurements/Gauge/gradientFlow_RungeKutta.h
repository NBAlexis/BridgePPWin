/*!
@file    $Id:: gradientFlow_RungeKutta.h #$

@brief

@author  Yusuke Namekawa (namekawa)

@date    $LastChangedDate:: 2017-02-24 18:35:38 #$

@version $LastChangedRevision: 1571 $
*/

#ifndef GRADIENTFLOW_RUNGEKUTTA_INCLUDED
#define GRADIENTFLOW_RUNGEKUTTA_INCLUDED

#include "Action/action.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//! GradientFlow_RungeKutta construction.

/*!
This class implements Runge-Kutta method
for GradientFlow.
[10 Oct 2014 Y.Namekawa]
*/

class BAPI GradientFlow_RungeKutta
{
protected:
    Bridge::VerboseLevel m_vl;

private:
    Action * m_action;
    int    m_Nprec;

    int m_Ndim;
    int m_Nvol;

public:
    GradientFlow_RungeKutta(Action *action, int Nprec, Bridge::VerboseLevel vl)
        : m_vl(vl),
        m_action(action),
        m_Nprec(Nprec),
        m_Ndim(CommonParameters::Ndim()),
        m_Nvol(CommonParameters::Nvol())
    {}

    virtual ~GradientFlow_RungeKutta() {}

    virtual void flow(double& t, double& Estep, Field_G& U) = 0;
    virtual int Norder_RK() const = 0;

private:
};
#endif
