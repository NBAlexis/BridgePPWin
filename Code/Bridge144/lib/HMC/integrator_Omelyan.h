/*!
        @file    $Id:: integrator_Omelyan.h #$

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 18:35:38 #$

        @version $LastChangedRevision: 1571 $
*/

#ifndef INTEGRATOR_OMELYAN_INCLUDED
#define INTEGRATOR_OMELYAN_INCLUDED

#include "integrator.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//! Omelyan integrator to compose MD integrator.

/*!
                             [28 Dec 2011 H.Matsufuru]
    (Coding history will be recovered from trac.)
    YAML is implemented.     [03 Mar 2013 Y.Namekawa]
 */


class Integrator_Omelyan : public Integrator
{
 public:
  static const std::string class_name;

 private:
  int    m_level;          // level number
  int    m_Nstep;          // number of steps
  double m_lambda;         // Omelyan's optimized parameter

  Integrator *m_update_p;  // momentum updator
  Integrator *m_update_U;  // link variable updator or next level integrator

 public:

  //! constructor
  Integrator_Omelyan(Integrator *update_p, Integrator *update_U)
    : m_level(0), m_Nstep(0),
      m_lambda(0.0),
      m_update_p(update_p), m_update_U(update_U)
  {
  }

  //! destructor
  ~Integrator_Omelyan() {}

  void set_parameters(const Parameters& params);
  void set_parameters(int level, int Nstep, double lambda_Omelyan);

  void set_parameter_level(const int level);

  void set_parameter_Nstep(const int Nstep);
  void set_parameter_Nsteps(const std::vector<int>& Nsteps);

  void set_parameter_lambda(const double lambda_omelyan);

  void evolve(const double step_size, Field_G& iP, Field_G& U);

  // cache management
  void invalidate_cache()
  {
    m_update_p->invalidate_cache();
    m_update_U->invalidate_cache();
  }
};
#endif
