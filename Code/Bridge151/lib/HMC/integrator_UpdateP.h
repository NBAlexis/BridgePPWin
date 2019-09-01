/*!
        @file    integrator_UpdateP.h

        @brief

        @author  Tatsumi Aoyama (aoym)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/

#ifndef INTEGRATOR_UPDATEP_INCLUDED
#define INTEGRATOR_UPDATEP_INCLUDED

#include "integrator.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//! Integrator of conjugate momenta for given link variables.

/*!
    This class BAPI updates conjugate momenta by force terms provided
    by action classes. it contains a list of actions that belong
    to a level of nested integrators.
    The value of force is cached in order to re-use it for the
    next half-step in multi-level integrator.
    cache should be invalidated when the link variable is updated
    through the call of invalidate_cache() method.

    [2015/04/17 T.Aoyama]
 */


class BAPI Integrator_UpdateP : public Integrator
{
 public:
  static const std::string class_name;

 private:
  std::vector<Action *> m_action;     // actions

  // work area
  Field_G m_workfield;

  // cache
  Field_G m_force;

  bool m_cache_valid;


 public:
  //! constructor when no director is necessary
  Integrator_UpdateP(const std::vector<Action *>& action)
    : m_action(action), m_cache_valid(false)
  {
    m_force.reset(CommonParameters::Nvol(), CommonParameters::Ndim());
    m_workfield.reset(CommonParameters::Nvol(), CommonParameters::Ndim());
  }

  //! destructor
  ~Integrator_UpdateP() {}

  void set_parameters(const Parameters& params);
  void set_parameters();

  void evolve(const double step_size, Field_G& iP, Field_G& U);

 public:
  void invalidate_cache() { m_cache_valid = false; }

  bool is_cache_valid() const { return m_cache_valid; }
  void cache_validated() { m_cache_valid = true; }
};
#endif
