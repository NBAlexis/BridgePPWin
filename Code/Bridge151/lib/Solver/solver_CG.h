/*!
        @file    solver_CG.h

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/

#ifndef SOLVER_CG_INCLUDED
#define SOLVER_CG_INCLUDED

#include "solver.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//! Standard Conjugate Gradient solver algorithm.

/*!
    This solver class BAPI implements the standard Conjugate Gradient
    solver algorithm.
                                   [22 Dec H.Matsufuru]
    (Coding history will be recovered from trac.)
    YAML is implemented.           [14 Nov 2012 Y.Namekawa]
    Multi-threaded.                [17 Jul 2014 Y.Namekawa]
    Introduce unique_ptr to avoid memory leaks.
                                   [21 Mar 2015 Y.Namekawa]
    Add restart.                   [22 Feb 2016 Y.Namekawa]
    Add flop_count.                [ 8 Aug 2016 Y.Namekawa]
    Add use_init_guess.           [ 7 Jul 2017 Y.Namekawa]
 */

class BAPI Solver_CG : public Solver
{
 public:
  static const std::string class_name;

 private:
  Fopr *m_fopr;

  int m_Niter;
  int m_Nrestart;
  double m_Stop_cond;
  bool m_use_init_guess;

  //- working area
  Field m_s, m_r, m_x, m_p;

  int m_Nrestart_count;
  int m_Nconv_count;

 public:
  Solver_CG(Fopr *fopr)
    : Solver(), m_fopr(fopr)
  {
    m_use_init_guess = false;
    m_Nrestart_count = 0;
    m_Nconv_count    = 0;
  }

  Solver_CG(unique_ptr<Fopr>& fopr)
    : Solver(), m_fopr(fopr.get())
  {
    m_use_init_guess = false;
    m_Nrestart_count = 0;
    m_Nconv_count    = 0;
  }

  ~Solver_CG() {}

  void set_parameters(const Parameters& params);
  void set_parameters(const int Niter, const int Nrestart, const double Stop_cond);
  void set_parameters(const int Niter, const int Nrestart, const double Stop_cond, const bool use_init_guess);

  void solve(Field& solution, const Field& source,
             int& Nconv, double& diff);

  Fopr *get_fopr() { return m_fopr; }

  double flop_count();

 private:
  void reset_field(const Field&);

  void solve_init(const Field&, double&);
  void solve_step(double&);

#ifdef USE_FACTORY
 private:
  static Solver *create_object(Fopr *fopr)
  {
    return new Solver_CG(fopr);
  }

 public:
  static bool register_factory()
  {
    return Solver::Factory::Register("CG", create_object);
  }
#endif
};
#endif
