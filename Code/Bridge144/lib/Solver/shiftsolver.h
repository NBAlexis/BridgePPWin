/*!
        @file    $Id:: shiftsolver.h #$

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 00:49:46 #$

        @version $LastChangedRevision: 1561 $
*/

#ifndef SHIFTSOLVER_INCLUDED
#define SHIFTSOLVER_INCLUDED

#include "defs.h"
#include "Parameters/parameters.h"
#include "Parameters/commonParameters.h"
#include "IO/bridgeIO.h"
#include "Field/field.h"
#include "Fopr/fopr.h"

//! Shiftsolver class as an abstract base class for multi-shift solvers.

class Shiftsolver
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

  virtual void solve(
    std::vector<Field>& solution,
    const std::vector<double>& shift,
    const Field& source,
    int& Nconv, double& diff) = 0;

 protected:
  Bridge::VerboseLevel m_vl;
};
#endif
