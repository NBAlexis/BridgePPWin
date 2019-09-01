/*!
        @file    integrator.h

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/

#ifndef INTEGRATOR_INCLUDED
#define INTEGRATOR_INCLUDED

#include "Action/action.h"
#include "Tools/director.h"

#include "IO/bridgeIO.h"

//! Base class BAPI of Integrator class BAPI family.

/*!
   This class BAPI defines the interface of Integrator-element_type classes.
                                        [25 Dec 2011 H.Matsufuru]
 */

class BAPI Integrator
{
 public:
  Integrator()
    : m_vl(CommonParameters::Vlevel()) {}

  virtual ~Integrator() {}

 private:
  // non-copyable
  Integrator(const Integrator&);
  Integrator& operator=(const Integrator&);

 public:
  virtual void set_parameters(const Parameters&) = 0;

  virtual void set_parameter_Nstep(const int Nstep) {}
  virtual void set_parameter_Nsteps(const std::vector<int>& Nsteps) {}

  void set_parameter_verboselevel(const Bridge::VerboseLevel vl) { m_vl = vl; }

  virtual void evolve(const double step_size, Field_G& iP, Field_G& U) = 0;


  // cache management
  virtual void invalidate_cache() = 0;


 protected:
  Bridge::VerboseLevel m_vl;
};
#endif
