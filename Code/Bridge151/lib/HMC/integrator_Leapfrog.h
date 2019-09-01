/*!
        @file    integrator_Leapfrog.h

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/

#ifndef INTEGRATOR_LEAPFROG_INCLUDED
#define INTEGRATOR_LEAPFROG_INCLUDED

#include "integrator.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//! Standard leapfrog integrator to compose MD integrator.

/*!
                                     [28 Dec 2011 H.Matsufuru]
    (Coding history will be recovered from trac.)
    YAML is implemented.             [03 Mar 2013 Y.Namekawa]
 */


class BAPI Integrator_Leapfrog : public Integrator
{
 public:
  static const std::string class_name;

 private:
  int m_level;             // level number
  int m_Nstep;             // number of steps

  Integrator *m_update_p;  // momentum updator
  Integrator *m_update_U;  // link variable updator or next level integrator

 public:

  //! constructor
  Integrator_Leapfrog(Integrator *update_p, Integrator *update_U)
    : m_level(0), m_Nstep(0),
    m_update_p(update_p), m_update_U(update_U)
  {
  }

  //! destructor
  ~Integrator_Leapfrog()
  {
  }

  void set_parameters(const Parameters& params);
  void set_parameters(const int level, const int Nstep);

  void set_parameter_level(const int level);

  void set_parameter_Nstep(const int Nstep);
  void set_parameter_Nsteps(const std::vector<int>& Nsteps);


  void evolve(const double step_size, Field_G& iP, Field_G& U);

  // cache management
  void invalidate_cache()
  {
    m_update_p->invalidate_cache();
    m_update_U->invalidate_cache();
  }
};
#endif
