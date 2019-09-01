/*!
        @file    projection.h

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/

#ifndef PROJECTION_INCLUDED
#define PROJECTION_INCLUDED

#include "Field/field_G.h"
#include "Parameters/parameters.h"

#ifdef USE_FACTORY
#include "Tools/factory.h"
#endif

//! base class BAPI for projection operator into gauge group.

/*!
                               [07 Apr 2012 H.Matsufuru]
 */

class BAPI Projection
{
 protected:
  Bridge::VerboseLevel m_vl;

 public:
  Projection()
    : m_vl(CommonParameters::Vlevel()) {}

  virtual ~Projection() {}

 private:
  // non-copyable
  Projection(const Projection&);
  Projection& operator=(const Projection&);

 public:
  //! projection V = P[alpha, C, U]
  virtual void project(Field_G& v,
                       const double alpha,
                       const Field_G& C, const Field_G& U) = 0;

  //! determination of fields for force calculation
  virtual void force_recursive(Field_G& Xi, Field_G& iTheta,
                               const double alpha, const Field_G& Sigmap,
                               const Field_G& C, const Field_G& U) = 0;

  virtual void set_parameters(const Parameters& param) = 0;

  void set_parameter_verboselevel(const Bridge::VerboseLevel vl) { m_vl = vl; }


#ifdef USE_FACTORY
 public:
  typedef Projection *(*ProductCreator)();
  typedef FactoryTemplate<Projection, ProductCreator> Factory;

  static Projection *New(const IdentifierType& subtype)
  {
    ProductCreator p = Factory::Find(subtype);

    return p ? (*p)() : 0;
  }

#ifdef USE_FACTORY_AUTOREGISTER
#else
  static bool init_factory();
#endif
#endif
};
#endif
