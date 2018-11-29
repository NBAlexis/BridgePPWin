/*!
@file    $Id:: staple.h #$

@brief

@author  Yusuke Namekawa (namekawa)
$LastChangedBy: aoym $

@date    $LastChangedDate:: 2017-02-24 00:49:46 #$

@version $LastChangedRevision: 1561 $
*/

#ifndef STAPLE_INCLUDED
#define STAPLE_INCLUDED

#include "Parameters/parameters.h"
#include "Field/field_G.h"

#include "IO/bridgeIO.h"

#ifdef USE_FACTORY
#include "Tools/factory.h"
#endif

//! Base class for Staple construction.

/*!
This class defines interface of Staple-type classes.
[24 Jan 2017 Y.Namekawa]
*/

class BAPI Staple
{
protected:
    Bridge::VerboseLevel m_vl;

public:
    Staple()
        : m_vl(CommonParameters::Vlevel()) {}

    virtual ~Staple() {}

private:
    // non-copyable
    Staple(const Staple&);
    Staple& operator=(const Staple&);

public:
    //! setting parameters.
    virtual void set_parameters(const Parameters& params) = 0;

    void set_parameter_verboselevel(const Bridge::VerboseLevel vl) { m_vl = vl; }

    //! constructs upper staple in mu-nu plane.
    virtual void upper(Field_G&, const Field_G&, const int mu, const int nu) = 0;

    //! constructs lower staple in mu-nu plane.
    virtual void lower(Field_G&, const Field_G&, const int mu, const int nu) = 0;

    //! constructs staple in mu-direction (summing up nu-direction).
    virtual void staple(Field_G&, const Field_G&, const int mu) = 0;

    //! calculates plaquette value.
    virtual double plaquette(const Field_G&) = 0;

    //! calculates spatial plaquette value.
    virtual double plaq_s(const Field_G&) = 0;

    //! calculates temporal plaquette value.
    virtual double plaq_t(const Field_G&) = 0;

#ifdef USE_FACTORY
public:
    typedef Staple *(*ProductCreator)();
    typedef FactoryTemplate<Staple, ProductCreator>   Factory;

    static Staple *New(const IdentifierType& subtype)
    {
        ProductCreator p = Factory::Find(subtype);

        return p ? (*p)() : 0;
    }
#endif
};
#endif
