/*!
        @file    $Id:: action.h #$

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 00:49:46 #$

        @version $LastChangedRevision: 1561 $
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


//! Base class of HMC action class family.

/*!
   This class defines interface of Action-type classes.
                                   [28 Dec 2011 H.Matsufuru]
   Factory is introduced.          [21 Mar 2015 Y.Namekawa]
 */

class Action
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
  typedef FactoryTemplate<Action, ProductCreator>   Factory;

  static Action *New(const IdentifierType& subtype)
  {
    ProductCreator p = Factory::Find(subtype);

    return p ? (*p)() : 0;
  }
#endif
};
#endif
