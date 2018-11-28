/*!
        @file    $Id:: force_G.h #$

        @brief

        @author  Yusuke Namekawa (namekawa)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 18:35:38 #$

        @version $LastChangedRevision: 1571 $
*/

#ifndef FORCE_G_INCLUDED
#define FORCE_G_INCLUDED

#include "Field/field_G.h"
#include "Parameters/parameters.h"

#ifdef USE_FACTORY
#include "Tools/factory.h"
#endif

//! Base class of gauge force calculation.

/*!
    This class defines the interface of gauge force calculation.
                                       [3 Mar 2016 Y.Namekawa]
*/

class Force_G
{
 protected:
  Field_G              *m_U;
  Bridge::VerboseLevel m_vl;

 public:
  Force_G()
    : m_U(0),
      m_vl(CommonParameters::Vlevel()) {}

  virtual ~Force_G() {}

 private:
  // non-copyable
  Force_G(const Force_G&);
  Force_G& operator=(const Force_G&);

 public:
  virtual void set_parameters(const Parameters&) = 0;

  void set_parameter_verboselevel(const Bridge::VerboseLevel vl)
  {
    m_vl = vl;
  }

  void set_config(Field *U)
  {
    m_U = (Field_G *)U;
  }

  void set_config(Field_G *U)
  {
    m_U = U;
  }

  virtual void force_core(Field&) = 0;

  virtual void force_core(Field& v, Field *U)
  {
    set_config(U);
    force_core(v);
  }

  virtual void force_core(Field& v, Field_G *U)
  {
    set_config(U);
    force_core(v);
  }

#ifdef USE_FACTORY
 public:
  typedef Force_G *(*ProductCreator)();
  typedef FactoryTemplate<Force_G, ProductCreator>   Factory;

  static Force_G *New(const IdentifierType& subtype)
  {
    ProductCreator p = Factory::Find(subtype);

    return p ? (*p)() : 0;
  }
#endif
};
#endif
