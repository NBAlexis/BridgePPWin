/*!
@file    $Id:: fopr_Rational.h #$

@brief

@author  Hideo Matsufuru (matsufuru)
$LastChangedBy: aoym $

@date    $LastChangedDate:: 2017-02-24 18:35:38 #$

@version $LastChangedRevision: 1571 $
*/


#ifndef FOPR_RATIONAL_INCLUDED
#define FOPR_RATIONAL_INCLUDED

#include "fopr.h"
#include "Field/field_G.h"
#include "Solver/shiftsolver_CG.h"
#include "Tools/math_Rational.h"

#include "IO/bridgeIO.h"
using Bridge::vout;


//! Fermion operator for rational approximation.

/*!
This class generates fermion operator with rational approximation
for a given fermion operator (given to the constructer).
Shift-solver is used which is at present set to the CG solver
explicitly.
[05 Dec 2011 H.Matsufuru]
(Coding history will be recovered from trac.)
YAML is implemented.             [14 Nov 2012 Y.Namekawa]
unique_ptr is introduced to avoid memory leaks
[21 Mar 2015 Y.Namekawa]
*/

class BAPI Fopr_Rational : public Fopr
{
public:
    static const std::string class_name;

private:
    int    m_Np;              // number of poles in rational approx.
    int    m_n_exp, m_d_exp;  // numerator and denominator of the exponent
    double m_x_min, m_x_max;  // valid range of approximate sign function
    int    m_Niter;           // max iteration of shiftsolver
    double m_Stop_cond;       // stopping condition of shift solver

    Fopr           *m_fopr;
    Shiftsolver_CG *m_solver;

    double              m_a0;
    std::vector<double> m_cl;
    std::vector<double> m_bl;
    std::vector<Field>  m_xq;

public:
    Fopr_Rational(Fopr *fopr)
        : Fopr(), m_fopr(fopr) {}

    Fopr_Rational(unique_ptr<Fopr>& fopr)
        : Fopr(), m_fopr(fopr.get()) {}

    ~Fopr_Rational()
    {
        delete m_solver;
    }

    void set_parameters(const Parameters& params);
    void set_parameters(int Np, int n_exp, int d_exp, double x_min, double x_max,
        int Niter, double Stop_cond);

    void set_config(Field *U)
    {
        m_fopr->set_config(U);
    }

    void set_config(unique_ptr<Field_G>& U)
    {
        m_fopr->set_config(U.get());
    }

    void mult(Field& v, const Field& f);

    void mult_dag(Field& v, const Field& f)
    {
        mult(v, f);
    }

    double func(double x);

    int field_nvol() { return m_fopr->field_nvol(); }
    int field_nin() { return m_fopr->field_nin(); }
    int field_nex() { return m_fopr->field_nex(); }

private:
    void init_parameters();
};
#endif
