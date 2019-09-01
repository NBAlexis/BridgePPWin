/*!
        @file    shiftsolver.h

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/

#ifndef SHIFTSOLVER_INCLUDED
#define SHIFTSOLVER_INCLUDED

#include "Fopr/fopr.h"

#include "IO/bridgeIO.h"


//! Shiftsolver class BAPI as an abstract base class BAPI for multi-shift solvers.

class BAPI Shiftsolver
{
 public:
  Shiftsolver()
    : m_vl(CommonParameters::Vlevel()) {}

  virtual ~Shiftsolver() {}

 private:
  // non-copyable
  Shiftsolver(const Shiftsolver&);
  Shiftsolver& operator=(const Shiftsolver&);

 public:
  virtual void set_parameters(const Parameters&) = 0;

  void set_parameter_verboselevel(const Bridge::VerboseLevel vl) { m_vl = vl; }

  virtual void solve(std::vector<Field>& solution,
                     const std::vector<double>& shift,
                     const Field& source,
                     int& Nconv, double& diff) = 0;

 protected:
  Bridge::VerboseLevel m_vl;
};
#endif
