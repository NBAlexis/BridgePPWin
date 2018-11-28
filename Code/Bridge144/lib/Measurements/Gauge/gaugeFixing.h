/*!
        @file    $Id:: gaugeFixing.h #$

        @brief

        @author  Yusuke Namekawa (namekawa)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 00:49:46 #$

        @version $LastChangedRevision: 1561 $
*/

#ifndef GAUGEFIXING_INCLUDED
#define GAUGEFIXING_INCLUDED

#include "Field/field_G.h"
#include "Parameters/parameters.h"
#include "IO/bridgeIO.h"

#ifdef USE_FACTORY
#include "Tools/factory.h"
#endif

//! gauge fixing.

/*
    This class fixes a gauge of the configuration
                                        [10 Oct 2012 Y.Namekawa]
    unique_ptr is introduced to avoid memory leaks
                                        [21 Mar 2015 Y.Namekawa]
    Staple and RandomNumbers are moved into gaugeFixing
                                        [30 Mar 2016 Y.Namekawa]
 */

class GaugeFixing
{
 public:
  GaugeFixing()
    : m_vl(CommonParameters::Vlevel()) {}

  virtual ~GaugeFixing() {}

 private:
  // non-copyable
  GaugeFixing(const GaugeFixing&);
  GaugeFixing& operator=(const GaugeFixing&);

 public:
  virtual void set_parameters(const Parameters& params) = 0;

  void set_parameter_verboselevel(const Bridge::VerboseLevel vl) { m_vl = vl; }

  virtual void fix(Field_G& Ufix, const Field_G& Uorg) = 0;

 protected:
  Bridge::VerboseLevel m_vl;


#ifdef USE_FACTORY
 public:
  typedef GaugeFixing *(*ProductCreator)();
  typedef FactoryTemplate<GaugeFixing, ProductCreator>   Factory;

  static GaugeFixing *New(const IdentifierType& subtype)
  {
    ProductCreator p = Factory::Find(subtype);

    return p ? (*p)() : 0;
  }
#endif
};
#endif
