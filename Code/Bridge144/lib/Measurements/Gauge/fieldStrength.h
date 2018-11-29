/*!
@file    $Id:: fieldStrength.h #$

@brief

@author  Yusuke Namekawa  (namekawa)
$LastChangedBy: aoym $

@date    $LastChangedDate:: 2017-02-24 18:35:38 #$

@version $LastChangedRevision: 1571 $
*/

#ifndef FIELDSTRENGTH_INCLUDED
#define FIELDSTRENGTH_INCLUDED

#include "staple_lex.h"
#include "Field/shiftField_lex.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//! field strength construction.

/*!
This class constructs a field strength, Fmunu,
defined by a clover leaf on the lattice.
[03 Mar 2016 Y.Namekawa]
*/

class BAPI FieldStrength
{
public:
    static const std::string class_name;

protected:
    Bridge::VerboseLevel m_vl;

private:
    ShiftField_lex m_shift;
    Staple_lex     m_staple;

public:
    FieldStrength()
        : m_vl(CommonParameters::Vlevel()) {}

    virtual ~FieldStrength() {}

private:
    // non-copyable
    FieldStrength(const FieldStrength&);
    FieldStrength& operator=(const FieldStrength&);

public:
    void construct_Fmunu_1x1(Field_G& Fmunu, const int mu, const int nu, const Field_G& U);
    void construct_Fmunu_1x2(Field_G& Fmunu, const int mu, const int nu, const Field_G& U);
};
#endif
