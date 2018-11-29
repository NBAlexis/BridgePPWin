/*!
@file    $Id:: shiftField_lex.h #$

@brief

@author  Hideo Matsufuru (matsufuru)
$LastChangedBy: namekawa $

@date    $LastChangedDate:: 2017-04-05 18:26:17 #$

@version $LastChangedRevision: 1608 $
*/

#ifndef SHIFTFIELD_LEX_INCLUDED
#define SHIFTFIELD_LEX_INCLUDED

#include "field.h"
#include "index_lex.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//! Methods to shift a field in the lexical site index.

/*!
This class BAPI defines the methods which shift a given Field
instance in the specified direction.
The forward shift means, e.g. in x-direction,
v(site) = w(site-\hat{mu}), where v is the shifted field
(output, first argument) and w the original field (input,
second argument).
[25 Dec 2011 H.Matsufuru]
*/

class BAPI ShiftField_lex {
public:
    static const std::string class_name;

private:
    int       m_Nx, m_Ny, m_Nz, m_Nt;
    Index_lex m_index_lex;

    Bridge::VerboseLevel m_vl;

public:
    ShiftField_lex() :
        m_Nx(CommonParameters::Nx()),
        m_Ny(CommonParameters::Ny()),
        m_Nz(CommonParameters::Nz()),
        m_Nt(CommonParameters::Nt()),
        m_vl(CommonParameters::Vlevel())
    {
    }

private:
    // non-copyable
    ShiftField_lex(const ShiftField_lex&);
    ShiftField_lex& operator=(const ShiftField_lex&);

public:
    void forward(Field&, const Field&, const int mu);
    void backward(Field&, const Field&, const int mu);

    void forward(Field&, const Field&, const int boundary_condition, const int mu);
    void backward(Field&, const Field&, const int boundary_condition, const int mu);

private:
    void up_x(Field *, const Field *, const int boundary_condition);
    void up_y(Field *, const Field *, const int boundary_condition);
    void up_z(Field *, const Field *, const int boundary_condition);
    void up_t(Field *, const Field *, const int boundary_condition);

    void dn_x(Field *, const Field *, const int boundary_condition);
    void dn_y(Field *, const Field *, const int boundary_condition);
    void dn_z(Field *, const Field *, const int boundary_condition);
    void dn_t(Field *, const Field *, const int boundary_condition);
};
#endif
