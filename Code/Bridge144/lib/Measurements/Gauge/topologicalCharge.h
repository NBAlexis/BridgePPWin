/*!
@file    $Id:: topologicalCharge.h #$

@brief

@author  Yusuke Namekawa  (namekawa)
$LastChangedBy: aoym $

@date    $LastChangedDate:: 2017-02-24 00:49:46 #$

@version $LastChangedRevision: 1561 $
*/

#ifndef TOPOLOGICALCHARGE_INCLUDED
#define TOPOLOGICALCHARGE_INCLUDED

#include "fieldStrength.h"


#include "IO/bridgeIO.h"
using Bridge::vout;

//! Topological Charge measurement.

/*!
This class measures a topological charge
defined by a clover leaf on the lattice.
[01 Jan 2014 Y.Namekawa]
*/


class BAPI TopologicalCharge
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


public:
    TopologicalCharge()
        : m_vl(CommonParameters::Vlevel())
    {
        m_filename_output = "stdout";
    }

    virtual ~TopologicalCharge() {}

private:
    // non-copyable
    TopologicalCharge(const TopologicalCharge&);
    TopologicalCharge& operator=(const TopologicalCharge&);

public:
    //! setting parameters.
    virtual void set_parameters(const Parameters& params);
    void set_parameters(double c_plaq, double c_rect);

    void set_parameter_verboselevel(const Bridge::VerboseLevel vl) { m_vl = vl; }

    //! main function to measure Topological Charge.
    double measure(Field_G& U);


private:
    double contract_epsilon_tensor(Field_G& Fmunu_1, Field_G& Fmunu_2);
};
#endif
