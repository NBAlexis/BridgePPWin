/*!
        @file    $Id:: solver.h #$

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 00:49:46 #$

        @version $LastChangedRevision: 1561 $
*/

#ifndef SOLVER_INCLUDED
#define SOLVER_INCLUDED

#include "Fopr/fopr.h"

#include "ResourceManager/threadManager_OpenMP.h"

#ifdef USE_FACTORY
#include "Tools/factory.h"
#endif



//! Base class for linear solver class family.

/*!
                               [28 Dec 2011 H.Matsufuru]
    Introduce unique_ptr to avoid memory leaks.
                               [21 Mar 2015 Y.Namekawa]
    Add restart.               [22 Feb 2016 Y.Namekawa]
    Add flop_count.            [ 8 Aug 2016 Y.Namekawa]
 */

class Solver
{
 public:
  Solver()
    : m_vl(CommonParameters::Vlevel()) {}

  virtual ~Solver() {}

 private:
  Solver(const Solver&);
  Solver& operator=(const Solver&);

 public:
  virtual void set_parameters(const Parameters& params) = 0;
  virtual void set_parameters(const int Niter, const int Nrestart, const double Stop_cond) = 0;

  void set_parameter_verboselevel(const Bridge::VerboseLevel vl) { m_vl = vl; }

  virtual void solve(Field& solution, const Field& source,
                     int& Nconv, double& diff) = 0;

  virtual Fopr *get_fopr() = 0;

  virtual double flop_count() = 0;

 protected:
  Bridge::VerboseLevel m_vl;

#ifdef USE_FACTORY
 public:
  typedef Solver *(*ProductCreator)(Fopr *);
  typedef FactoryTemplate<Solver, ProductCreator>   Factory;

  static Solver *New(const IdentifierType& subtype, Fopr *fopr)
  {
    ProductCreator p = Factory::Find(subtype);

    return p ? (*p)(fopr) : 0;
  }

  static Solver *New(const IdentifierType& subtype, unique_ptr<Fopr>& fopr)
  {
    ProductCreator p = Factory::Find(subtype);

    return p ? (*p)(fopr.get()) : 0;
  }
#endif
};
#endif
