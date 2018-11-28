/*!
        @file    $Id:: forceSmear.h #$

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 00:49:46 #$

        @version $LastChangedRevision: 1561 $
*/

#ifndef FORCESMEAR_INCLUDED
#define FORCESMEAR_INCLUDED

#include "Smear/projection.h"

#include "IO/bridgeIO.h"

#ifdef USE_FACTORY
#include "Tools/factory.h"
#endif


//! Base class for force calculation of smeared operators.

/*!
                                     [28 Dec 2011 H.Matsufuru]
    unique_ptr is introduced to avoid memory leaks
                                     [21 Mar 2015 Y.Namekawa]
 */

class ForceSmear
{
 public:

  ForceSmear()
    : m_vl(CommonParameters::Vlevel()) {}

  virtual ~ForceSmear() {}

 private:
  // non-copyable
  ForceSmear(const ForceSmear&);
  ForceSmear& operator=(const ForceSmear&);

 public:
  virtual void set_parameters(const Parameters&) = 0;

  void set_parameter_verboselevel(const Bridge::VerboseLevel vl) { m_vl = vl; }

  virtual void force_udiv(Field_G&, const Field_G&, const Field_G&) {}

 protected:
  Bridge::VerboseLevel m_vl;


#ifdef USE_FACTORY
 public:
  typedef ForceSmear *(*ProductCreator)(Projection *);
  typedef FactoryTemplate<ForceSmear, ProductCreator>   Factory;

  static ForceSmear *New(const IdentifierType& subtype, Projection *proj)
  {
    ProductCreator p = Factory::Find(subtype);

    return p ? (*p)(proj) : 0;
  }

  static ForceSmear *New(const IdentifierType& subtype, unique_ptr<Projection>& proj)
  {
    ProductCreator p = Factory::Find(subtype);

    return p ? (*p)(proj.get()) : 0;
  }
#endif
};
#endif
