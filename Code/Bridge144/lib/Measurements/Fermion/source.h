/*!
        @file    $Id:: source.h #$

        @brief

        @author  Tatsumi Aoyama  (aoym)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 00:49:46 #$

        @version $LastChangedRevision: 1561 $
*/

#ifndef SOURCE_INCLUDED
#define SOURCE_INCLUDED

#include "Parameters/parameters.h"
#include "Field/field.h"
#include "IO/bridgeIO.h"

#ifdef USE_FACTORY
#include "Tools/factory.h"
#endif

//! Base class of source for a linear solver.

/*!
    (Coding history will be recovered from trac.)
    Parameters_Source_All is implemented for source_selector.
                                      [2 Feb 2013 Y.Namekawa]
    Bug fix: initialization of verbose level was added.
                                     [23 May 2016 H.Matsufuru]
 */

class Source
{
 protected:
  Bridge::VerboseLevel m_vl;

 public:
  Source()
    : m_vl(CommonParameters::Vlevel()) {}

  virtual ~Source() {}

 private:
  // non-copyable
  Source(const Source&);
  Source& operator=(const Source&);

 public:

  virtual void set_parameters(const Parameters&) = 0;

  void set_parameter_verboselevel(const Bridge::VerboseLevel vl) { m_vl = vl; }

  virtual void set(Field&, int) = 0;

#ifdef USE_FACTORY
 public:
  typedef Source *(*ProductCreator)();
  typedef FactoryTemplate<Source, ProductCreator>   Factory;

  static Source *New(const IdentifierType& subtype)
  {
    ProductCreator p = Factory::Find(subtype);

    return p ? (*p)() : 0;
  }
#endif
};
#endif /* SOURCE_INCLUDED */
