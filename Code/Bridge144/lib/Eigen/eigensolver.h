/*!
        @file    $Id:: eigensolver.h #$

        @brief

        @author  Tatsumi Aoyama (aoym)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 00:49:46 #$

        @version $LastChangedRevision: 1561 $
*/

#ifndef EIGENSOLVER_INCLUDED
#define EIGENSOLVER_INCLUDED

#include "defs.h"
#include "Parameters/parameters.h"
#include "Parameters/commonParameters.h"
#include "Field/field.h"

//! Eigensolver class for abstract base class of eigen solvers.

/**
   Eigensolver class provides an abstract base class for solvers
   of eigenvalues and eigenvectors.
 */

class Eigensolver
{
 public:

  Eigensolver()
    : m_vl(CommonParameters::Vlevel()) {}

  virtual ~Eigensolver() {}

 private:
  // non-copyable
  Eigensolver(const Eigensolver&);
  Eigensolver& operator=(const Eigensolver&);

 public:
  virtual void set_parameters(const Parameters&) = 0;

  void set_parameter_verboselevel(Bridge::VerboseLevel vl) { m_vl = vl; }

  virtual void solve(std::vector<double>& TDa,
                     std::vector<Field>& vk,
                     int& Nsbt, int& Nconv,
                     const Field& b) = 0;

 protected:
  Bridge::VerboseLevel m_vl;
};
#endif
