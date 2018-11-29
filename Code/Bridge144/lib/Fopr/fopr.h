/*!
@file    $Id:: fopr.h #$

@brief

@author  Hideo Matsufuru (matsufuru)
$LastChangedBy: matufuru $

@date    $LastChangedDate:: 2017-11-30 20:28:21 #$

@version $LastChangedRevision: 1682 $
*/


#ifndef FOPR_INCLUDED
#define FOPR_INCLUDED

#include "Field/field_F.h"
#include "Parameters/parameters.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

#ifdef USE_FACTORY
#include "Tools/factory.h"
#include "Tools/director.h"
#endif


//! Base class of fermion operator family.

/*!
In Bridge-OF, the fermion operator implies an operator which
transforms a field to other field, irrespective of physical
formulation.
This class defines the interface of the fermion operators.
At present, void functions mult(v,w) and mult_dag(v,w) is
not purely virtual, because some of subclass have not
implemented them yet.
[20 Dec 2011 H.Matsufuru]
unique_ptr is introduced to avoid memory leaks.
[21 Mar 2015 Y.Namekawa]
parameters factory is introduced. [21 Apr 2015 Y.Namekawa]
mult_gm5 is added.                [30 Nov 2017 H.Matsufuru]
*/

class BAPI Fopr
{
public:

    Fopr()
        : m_vl(CommonParameters::Vlevel()) {}

    virtual ~Fopr() {}

private:
    // non-copyable
    Fopr(const Fopr&);
    Fopr& operator=(const Fopr&);

public:
    virtual void set_parameters(const Parameters&) = 0;

    void set_parameter_verboselevel(const Bridge::VerboseLevel vl) { m_vl = vl; }

    //! setting pointer to the gauge configuration.
    virtual void set_config(Field *) = 0;
    virtual void set_config(unique_ptr<Field_G>&) = 0;

    //! \brief multiplies fermion operator to a given field (2nd argument)
    //   and set the resultant field to the 1st argument.
    virtual void mult(Field&, const Field&) = 0;

    //! hermitian conjugate of mult(Field&, const Field&).
    virtual void mult_dag(Field&, const Field&) {}

    //! execute mult with specified mode (unchanging internal mode). [23 May 2016 H.Matsufuru].
    virtual void mult(Field&, const Field&, const std::string mode) {}

    //! execute mult_dag with specified mode (unchanging internal mode). [23 May 2016 H.Matsufuru].
    virtual void mult_dag(Field&, const Field&, const std::string mode) {}

    //! gamma_5 multiplication. [31 Mar 2017 H.Matsufuru]
    virtual void mult_gm5(Field&, const Field&)
    {
        vout.crucial(m_vl, "Fopr: mult_gm5 not implemented.\n");
        exit(EXIT_FAILURE);
    }

    //! nearest neighbor hopping term: temporary entry [H.Matsufuru]
    virtual void mult_up(int, Field&, const Field&) {}
    virtual void mult_dn(int, Field&, const Field&) {}

    //! \brief setting the mode of multiplication if necessary.
    //!  Default implementation here is just to avoid irrelevant call.
    virtual void set_mode(std::string mode)
    {
        vout.general(m_vl, "Fopr: set_mode not implemented.\n");
    }

    //! only for Fopr_Overlap
    // virtual void set_lowmodes(int Nsbt, std::vector<double> *ev,
    //                           std::vector<Field> *vk);

    virtual std::string get_mode() const
    {
        vout.general(m_vl, "Fopr: get_mode not implemented.\n");
        return std::string();
    }

    //! returns the volume for which the fermion operator is defined.
    virtual int field_nvol() = 0;

    //! returns the on-site d.o.f. for which the fermion operator is defined.
    virtual int field_nin() = 0;

    //! returns the external d.o.f. for which the fermion operator is defined.
    virtual int field_nex() = 0;

    //! returns the flops per site.
    virtual double flop_count() { return 0.0; }

    //! returns the flops per site for specified mode. [23 May 2016 H.Matsufuru]
    virtual double flop_count(const std::string mode)
    {
        return 0.0;
    }

protected:
    Bridge::VerboseLevel m_vl;


#ifdef USE_FACTORY
public:
    typedef Fopr *(*ProductCreator_noarg)();
    typedef Fopr *(*ProductCreator_fopr)(Fopr *fopr);
    typedef Fopr *(*ProductCreator_fopr_director)(Fopr *fopr, Director *director);
    typedef Fopr *(*ProductCreator_string)(const std::string& arg);

    typedef FactoryTemplate<Fopr, ProductCreator_noarg>           Factory_noarg;
    typedef FactoryTemplate<Fopr, ProductCreator_fopr>            Factory_fopr;
    typedef FactoryTemplate<Fopr, ProductCreator_fopr_director>   Factory_fopr_director;
    typedef FactoryTemplate<Fopr, ProductCreator_string>          Factory_string;

    static Fopr *New(const IdentifierType& subtype)
    {
        ProductCreator_noarg p = Factory_noarg::Find(subtype);

        return p ? (*p)() : 0;
    }

    static Fopr *New(const IdentifierType& subtype, Fopr *fopr)
    {
        ProductCreator_fopr p = Factory_fopr::Find(subtype);

        return p ? (*p)(fopr) : 0;
    }

    static Fopr *New(const IdentifierType& subtype, unique_ptr<Fopr>& fopr)
    {
        ProductCreator_fopr p = Factory_fopr::Find(subtype);

        return p ? (*p)(fopr.get()) : 0;
    }

    static Fopr *New(const IdentifierType& subtype, Fopr *fopr, Director *director)
    {
        ProductCreator_fopr_director p = Factory_fopr_director::Find(subtype);

        return p ? (*p)(fopr, director) : 0;
    }

    static Fopr *New(const IdentifierType& subtype, unique_ptr<Fopr>& fopr, unique_ptr<Director>& director)
    {
        ProductCreator_fopr_director p = Factory_fopr_director::Find(subtype);

        return p ? (*p)(fopr.get(), director.get()) : 0;
    }

    static Fopr *New(const IdentifierType& subtype, const std::string& arg)
    {
        ProductCreator_string p = Factory_string::Find(subtype);

        return p ? (*p)(arg) : 0;
    }
#endif
};
#endif
