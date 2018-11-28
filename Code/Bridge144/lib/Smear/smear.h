/*!
        @file    $Id:: smear.h #$

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 00:49:46 #$

        @version $LastChangedRevision: 1561 $
*/


#ifndef SMEAR_INCLUDED
#define SMEAR_INCLUDED

#include "Parameters/parameters.h"
#include "projection.h"

#include "IO/bridgeIO.h"

#ifdef USE_FACTORY
#include "Tools/factory.h"
#endif


//! base class for smearing of link variables.

/*!
                            [28 Dec 2011 H.Matsufuru]
    unique_ptr is introduced to avoid memory leaks
                            [21 Mar 2015 Y.Namekawa]
 */

class Smear
{
 protected:
  Bridge::VerboseLevel m_vl;

 public:
  Smear()
    : m_vl(CommonParameters::Vlevel()) {}

  virtual ~Smear() {}

 private:
  Smear(const Smear&);
  Smear& operator=(const Smear&);

 public:
  virtual void smear(Field_G&, const Field_G&) = 0;

  virtual void set_parameters(const Parameters&) = 0;

  void set_parameter_verboselevel(const Bridge::VerboseLevel vl) { m_vl = vl; }

#ifdef USE_FACTORY
 public:
  typedef Smear *(*ProductCreator)(Projection *);
  typedef FactoryTemplate<Smear, ProductCreator>   Factory;

  static Smear *New(const IdentifierType& subtype, Projection *proj)
  {
    ProductCreator p = Factory::Find(subtype);

    return p ? (*p)(proj) : 0;
  }

  static Smear *New(const IdentifierType& subtype, unique_ptr<Projection>& proj)
  {
    ProductCreator p = Factory::Find(subtype);

    return p ? (*p)(proj.get()) : 0;
  }
#endif
};
#endif
