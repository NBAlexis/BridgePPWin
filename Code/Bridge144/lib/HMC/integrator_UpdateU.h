/*!
        @file    $Id:: integrator_UpdateU.h #$

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 18:35:38 #$

        @version $LastChangedRevision: 1571 $
*/

#ifndef INTEGRATOR_UPDATEU_INCLUDED
#define INTEGRATOR_UPDATEU_INCLUDED

#include "integrator.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//! Integrator of link variable for a given conjugate momenta.

/*!
    In the present implementation, exponential of matrix is determined
    by Taylor series, and whose degree (m_Nprec) is explicitly
    specified in the class definition.
                                          [25 Dec 2011 H.Matsufuru]
    (Coding history will be recovered from trac.)
    YAML is implemented.                  [03 Mar 2013 Y.Namekawa]
 */


class BAPI Integrator_UpdateU : public Integrator
{
 public:
  static const std::string class_name;

  static const int default_Nprec = 8;

 private:
  int m_Nprec;                          // precision of exponential series

  std::vector<Director *>   m_director; // directors
  std::vector<Integrator *> m_integs;   // list of update_p to send notification

 public:
  Integrator_UpdateU(const std::vector<Director *>& director = std::vector<Director *>())
    : m_Nprec(default_Nprec), m_director(director), m_integs()
  {}

  //! destructor
  ~Integrator_UpdateU() {}

  void set_parameters(const Parameters& params);
  void set_parameters(int Nprec);

  void set_parameter_Nprec(const int Nprec);


  void evolve(const double step_size, Field_G& iP, Field_G& U);


  // cache management
  void invalidate_cache()
  {
    // lowest level. do nothing.
  }

  void append_notify(Integrator *const integ)
  {
    m_integs.push_back(integ);
  }

  void notify_update();
};
#endif
