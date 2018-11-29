/*!
@file    $Id:: polyakovLoop.h #$

@brief

@author  Hideo Matsufuru (matsufuru)
$LastChangedBy: aoym $

@date    $LastChangedDate:: 2017-02-24 18:35:38 #$

@version $LastChangedRevision: 1571 $
*/

#ifndef POLYAKOVLOOP_INCLUDED
#define POLYAKOVLOOP_INCLUDED

#include <cassert>
#include "Parameters/parameters.h"
#include "Field/field_G.h"
#include "Field/shiftField_lex.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//! Polyakov loop measurement.

/*!
This class determines the Polyakov loop for a given gauge
configuration.
It is planed that the Polyakov loop correlators are also
measured, but still not implemented.
For the latter case, set_parameters() is prepared (the definition
of correlator constellation is the same as WilsonLoop).
[22 Aug 2012 H.Matsufuru]
(Coding history will be recovered from trac.)
YAML is implemented.         [14 Nov 2012 Y.Namekawa]
*/


class BAPI PolyakovLoop
{
public:
    static const std::string class_name;

protected:
    Bridge::VerboseLevel m_vl;

private:
    std::string m_filename_output;

    //! parameters set by user
    int m_Nspc_size;  //!< spatial size of loop
    int m_Ntype;      //!< number of measured loop-type

                      //! internal data members
    int m_Ntype_max;  //!< maximum size of loop-type
    int m_Nx_ext;     //!< size of extended gauge config.
    int m_Ny_ext;     //!< size of extended gauge config.
    int m_Nz_ext;     //!< size of extended gauge config.
    int m_Nt_ext;     //!< size of extended gauge config.
    int m_Nvol_ext;   //!< volume of extended gauge config.

    typedef std::vector<int>   unitvec;
    std::vector<unitvec> m_Nunit;
    std::vector<int>     m_Nmax;

public:

    PolyakovLoop()
        : m_vl(CommonParameters::Vlevel()), m_Nspc_size(0), m_Ntype(0)
    {
        init();
    }

    virtual ~PolyakovLoop() {}

private:
    // non-copyable
    PolyakovLoop(const PolyakovLoop&);
    PolyakovLoop& operator=(const PolyakovLoop&);

public:
    //! setting parameters: only for Polyakov loop correlators.
    virtual void set_parameters(const Parameters& params);
    void set_parameters(int Nspc_size, int Ntype);

    void set_parameter_verboselevel(const Bridge::VerboseLevel vl) { m_vl = vl; }

    //! Polyakov loop measurement
    dcomplex measure_ploop(Field_G& U);

    //! Polyakov loop correlator measurement (not implemented).
    double measure_ploop_corr(Field_G& U);

private:

    //! Polyakov loop measurement
    void calc_ploop(Field_G& P, Field_G& U);

    //! initial setup independent of parameters.
    void init();
};
#endif
