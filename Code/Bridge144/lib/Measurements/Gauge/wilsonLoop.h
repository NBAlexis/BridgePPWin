/*!
@file    $Id:: wilsonLoop.h #$

@brief

@author  Hideo Matsufuru (matsufuru)
$LastChangedBy: aoym $

@date    $LastChangedDate:: 2017-02-24 18:35:38 #$

@version $LastChangedRevision: 1571 $
*/

#ifndef WILSONLOOP_INCLUDED
#define WILSONLOOP_INCLUDED

#include "Field/field_G.h"
#include "Parameters/parameters.h"
#include "Field/shiftField_eo.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//! Wilson loop measurement.

/*!
This class measures Wilson loops in spatial off-axis directions
as well as on-axis for a given gauge configuration.
The construction of off-axis loops originates from the Fortran code
written by Takashi Umeda (1997), while details were modified.
The calculation is performed in the following steps.
(1) on each node, extended gauge configuration necessary to
determine the whole Wilson loop is set.
(2) temporal gauge fixing is performed.
(3) spatial products of link variables are obtained, from which
Wilson loops are constructed.
[28 May 2012 H.Matsufuru]
(Coding history will be recovered from trac.)
YAML is implemented.                  [14 Nov 2012 Y.Namekawa]
*/


class BAPI WilsonLoop
{
public:
    static const std::string class_name;

protected:
    Bridge::VerboseLevel m_vl;

private:
    std::string m_filename_output;

    //! parameters set by user
    int m_Nspc_size;  //!< spatial size of loop
    int m_Ntmp_size;  //!< spatial size of loop
    int m_Ntype;      //!< number of measured loop-type

                      //! internal data members
    int m_Ntype_max;  //!< maximum size of loop-type
    int m_Nx_ext;     //!< size of extended gauge config.
    int m_Ny_ext;     //!< size of extended gauge config.
    int m_Nz_ext;     //!< size of extended gauge config.
    int m_Nt_ext;     //!< size of extended gauge config.
    int m_Nvol_ext;   //!< volume of extended gauge config.

    typedef std::vector<int>   unit_vec;
    std::vector<unit_vec> m_Nunit;
    std::vector<int>      m_Nmax;

public:
    WilsonLoop()
        : m_vl(CommonParameters::Vlevel())
    {
        init();
    }

    virtual ~WilsonLoop() {}

private:
    // non-copyable
    WilsonLoop(const WilsonLoop&);
    WilsonLoop& operator=(const WilsonLoop&);

public:
    //! setting parameters.
    virtual void set_parameters(const Parameters& params);
    void set_parameters(int Nspc_size, int Ntmp_size, int Ntype);

    void set_parameter_verboselevel(const Bridge::VerboseLevel vl) { m_vl = vl; }

    //! main function to measure Wilson loops.
    double measure(Field_G& U);

    //! index for Wilson loop variable.
    int index_wloop(int i_spc, int i_tmp, int i_type)
    {
        return i_spc + m_Nspc_size * (i_tmp + m_Ntmp_size * i_type);
    }

private:
    //! initial setup independent of parameters.
    void init();

    //! temporal gauge fixing of extended gauge field.
    double calc_wloop(Field_G& Uspc, int t_ext);

    //! redefinition of product of spatial link variables.
    void redef_Uspc(Field_G& Uspc, Field_G& Uext,
        int j, int nu, std::vector<int>& unit_v);

    //! setup of extended gauge field.
    void set_extfield(Field_G& Uext, Field_G& Uorg);

    //! temporal gauge fixing of extended gauge field.
    void gfix_temporal(Field_G& Uext);
};
#endif
