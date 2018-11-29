/*!
@file    $Id:: energyDensity.h #$

@brief

@author  Yusuke Namekawa  (namekawa)
$LastChangedBy: aoym $

@date    $LastChangedDate:: 2017-02-24 18:35:38 #$

@version $LastChangedRevision: 1571 $
*/

#ifndef ENERGYDENSITY_INCLUDED
#define ENERGYDENSITY_INCLUDED

#include "fieldStrength.h"


#include "IO/bridgeIO.h"
using Bridge::vout;

//! energy density measurement.

/*!
This class measures an energy density.
[03 Mar 2016 Y.Namekawa]
*/



class BAPI EnergyDensity
{
public:
    static const std::string class_name;

protected:
    Bridge::VerboseLevel m_vl;

private:
    std::string m_filename_output;

    double m_c_plaq;
    double m_c_rect;

    FieldStrength m_field_strength;
    Staple_lex    m_staple;


public:
    EnergyDensity()
        : m_vl(CommonParameters::Vlevel())
    {
        m_filename_output = "stdout";
    }

    virtual ~EnergyDensity() {}

private:
    // non-copyable
    EnergyDensity(const EnergyDensity&);
    EnergyDensity& operator=(const EnergyDensity&);


public:
    //! setting parameters.
    virtual void set_parameters(const Parameters& params);
    void set_parameters(double c_plaq, double c_rect);

    void set_parameter_verboselevel(const Bridge::VerboseLevel vl) { m_vl = vl; }

    double E_plaq(const Field_G& U);
    double E_clover(const Field_G& U);
};
#endif
