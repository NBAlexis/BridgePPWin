/*!
@file    $Id:: fprop_Standard_eo.h #$

@brief

@author  Satoru Ueda  (sueda)
$LastChangedBy: aoym $

@date    $LastChangedDate:: 2017-02-24 18:35:38 #$

@version $LastChangedRevision: 1571 $
*/

#ifndef FPROP_STANDARD_EO_INCLUDED
#define FPROP_STANDARD_EO_INCLUDED

#include "fprop.h"
#include "Field/field_G.h"

#include "Fopr/fopr_eo.h"
#include "Field/index_eo.h"

#include "Solver/solver.h"

//! Get quark propagator for Fopr with even-odd site index.

/*!
This is temporary implementation.
[28 Dec 2011 H.Matsufuru]
(Coding history will be recovered from trac.)
Modify this code to work.           [03 Mar 2013 Y.Namekawa]
Introduce unique_ptr to avoid memory leaks.
[21 Mar 2015 Y.Namekawa]
Add flop_count.                     [ 8 Aug 2016 Y.Namekawa]
*/

class BAPI Fprop_Standard_eo : public Fprop
{
public:
    static const std::string class_name;

private:
    Solver * m_solver;
    Field_G *m_Ueo;

    Index_eo *m_index;

public:
    Fprop_Standard_eo(Solver *solver)
        : Fprop(), m_solver(solver)
    {
        m_index = new Index_eo;
        m_Ueo = new Field_G(CommonParameters::Nvol(), CommonParameters::Ndim());
    }

    Fprop_Standard_eo(unique_ptr<Solver>& solver)
        : Fprop(), m_solver(solver.get())
    {
        m_index = new Index_eo;
        m_Ueo = new Field_G(CommonParameters::Nvol(), CommonParameters::Ndim());
    }

    ~Fprop_Standard_eo()
    {
        delete m_index;
        delete m_Ueo;
    }

private:
    // non-copyable
    Fprop_Standard_eo(const Fprop_Standard_eo&);
    Fprop_Standard_eo& operator=(const Fprop_Standard_eo&);

public:
    void set_config(Field *);

    void invert_D(Field&, const Field&, int&, double&);
    void invert_DdagD(Field&, const Field&, int&, double&);

    double flop_count();
};
#endif
