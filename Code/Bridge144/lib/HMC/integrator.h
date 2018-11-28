/*!
        @file    $Id:: integrator.h #$

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 00:49:46 #$

        @version $LastChangedRevision: 1561 $
*/

#ifndef INTEGRATOR_INCLUDED
#define INTEGRATOR_INCLUDED

#include "Action/action.h"
#include "Tools/director.h"

#include "IO/bridgeIO.h"

//! Base class of Integrator class family.

/*!
   This class defines the interface of Integrator-type classes.
                                        [25 Dec 2011 H.Matsufuru]
 */

class Integrator
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
