/*!
        @file    eigensolver.h

        @brief

        @author  Tatsumi Aoyama (aoym)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/

#ifndef EIGENSOLVER_INCLUDED
#define EIGENSOLVER_INCLUDED

#include "Field/field.h"
#include "Parameters/parameters.h"

//! Eigensolver class BAPI for abstract base class BAPI of eigen solvers.

/**
   Eigensolver class BAPI provides an abstract base class BAPI for solvers
   of eigenvalues and eigenvectors.
 */

class BAPI Eigensolver
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

  void set_parameter_verboselevel(const Bridge::VerboseLevel vl) { m_vl = vl; }

  virtual void solve(std::vector<double>& TDa,
                     std::vector<Field>& vk,
                     int& Nsbt, int& Nconv,
                     const Field& b) = 0;

 protected:
  Bridge::VerboseLevel m_vl;
};
#endif
