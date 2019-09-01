/*!
        @file    action.h

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/


#ifndef ACTION_INCLUDED
#define ACTION_INCLUDED

#include "Field/field_G.h"
#include "Tools/randomNumbers.h"
#include "Parameters/parameters.h"
#include "IO/bridgeIO.h"

#ifdef USE_FACTORY
#include "Tools/factory.h"
#endif


//! Base class BAPI of HMC action class BAPI family.

/*!
   This class BAPI defines interface of Action-element_type classes.
                                   [28 Dec 2011 H.Matsufuru]
   Factory is introduced.          [21 Mar 2015 Y.Namekawa]
 */

class BAPI Action
{
 public:

  Action()
    : m_vl(CommonParameters::Vlevel()) {}

  virtual ~Action() {}

 private:
  // non-copyable
  Action(const Action&);
  Action& operator=(const Action&);

 public:
  virtual void set_parameters(const Parameters& param) = 0;

  void set_parameter_verboselevel(const Bridge::VerboseLevel vl) { m_vl = vl; }

  //! setting pointer to the gauge configuration.
  virtual void set_config(Field *U) = 0;

  //! Langevis step.
  virtual double langevin(RandomNumbers *) = 0;

  //! calculate Hamiltonian of this action term.
  virtual double calcH() = 0;

  //! returns force for molcular dynamical update of conjugate momenta.
  //virtual const Field force() = 0;
  virtual void force(Field&) = 0;

  virtual void force(Field& v, Field& U)
  {
    set_config(&U);
    force(v);
  }

 protected:
  Bridge::VerboseLevel m_vl;

#ifdef USE_FACTORY
 public:
  typedef Action *(*ProductCreator)();
  typedef FactoryTemplate<Action, ProductCreator> Factory;

  static Action *New(const IdentifierType& subtype)
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
