/*!
        @file    solver.h

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/

#ifndef SOLVER_INCLUDED
#define SOLVER_INCLUDED

#include "Fopr/fopr.h"

#include "ResourceManager/threadManager_OpenMP.h"

#ifdef USE_FACTORY
#include "Tools/factory.h"
#endif

//! Base class BAPI for linear solver class BAPI family.

/*!
                               [28 Dec 2011 H.Matsufuru]
    Introduce unique_ptr to avoid memory leaks.
                               [21 Mar 2015 Y.Namekawa]
    Add restart.               [22 Feb 2016 Y.Namekawa]
    Add flop_count.            [ 8 Aug 2016 Y.Namekawa]
    Add use_init_guess.       [ 7 Jul 2017 Y.Namekawa]
 */

class BAPI Solver
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
  virtual void set_parameters(const int Niter, const int Nrestart, const double Stop_cond, const bool use_init_guess) = 0;

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
  typedef FactoryTemplate<Solver, ProductCreator> Factory;

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

#ifdef USE_FACTORY_AUTOREGISTER
#else
  static bool init_factory();
#endif
#endif
};
#endif
