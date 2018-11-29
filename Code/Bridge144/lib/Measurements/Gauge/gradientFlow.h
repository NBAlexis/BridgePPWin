/*!
@file    $Id:: gradientFlow.h #$

@brief

@author  Sinya Aoki

@date    $LastChangedDate:: 2017-02-24 18:35:38 #$

@version $LastChangedRevision: 1571 $
*/

#ifndef GRADIENTFLOW_INCLUDED
#define GRADIENTFLOW_INCLUDED

#include "gradientFlow_RungeKutta_1st.h"
#include "gradientFlow_RungeKutta_2nd.h"
#include "gradientFlow_RungeKutta_3rd.h"
#include "gradientFlow_RungeKutta_4th.h"
#include "gradientFlow_AdaptiveRungeKutta.h"

#include "staple_lex.h"


#include "IO/bridgeIO.h"
using Bridge::vout;

//! GradientFlow construction.

/*!
This is written by S. Aoki, based on staple.h
[01 July 2012 S.Aoki]
(Coding history will be recovered from trac.)
YAML is implemented.       [14 Nov 2012 Y.Namekawa]
4th order Runge-Kutta in commutator-free method
formulated by E.Celledoni et al. FGCS 19, 341 (2003),
as well as 1st and 2nd order Runge-Kutta are implemented.
[10 Oct 2014 Y.Namekawa]
unique_ptr is introduced to avoid memory leaks
[21 Mar 2015 Y.Namekawa]
Adaptive stepsize control is implemented.
[01 May 2015 Y.Namekawa]
*/


class BAPI GradientFlow
{
public:
    static const std::string class_name;

protected:
    Bridge::VerboseLevel m_vl;

private:
    int    m_Norder_RK;
    double m_Estep;
    int    m_Nprec;
    double m_tolerance;
    double m_safety;
    bool   m_is_adaptive;

    Action     *m_action;
    Staple_lex m_staple;

    GradientFlow_RungeKutta *m_impl;

public:
    GradientFlow(Action *action);
    GradientFlow(unique_ptr<Action>& action);

    ~GradientFlow();

private:
    // non-copyable
    GradientFlow(const GradientFlow&);
    GradientFlow& operator=(const GradientFlow&);

public:
    void set_parameters(const Parameters& params);
    void set_parameters(const int Norder_RK,
        const double Estep, const int Nprec,
        const int adaptive,
        const double tolerance, const double safety);
    void set_parameter_Norder_RK(const int order, const bool is_adaptive);

    void set_parameter_verboselevel(const Bridge::VerboseLevel vl) { m_vl = vl; }

    void initialize();

    double evolve(double& t, Field_G& U);
};
#endif
