/*!
        @file    $Id:: solver_CG.h #$

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 18:35:38 #$

        @version $LastChangedRevision: 1571 $
*/

#ifndef SOLVER_CG_INCLUDED
#define SOLVER_CG_INCLUDED

#include "solver.h"

#include "IO/bridgeIO.h"
using Bridge::vout;


//! Standard Conjugate Gradient solver algorithm.

/*!
    This solver class implements the standard Conjugate Gradient
    solver algorithm.
                                   [22 Dec H.Matsufuru]
    (Coding history will be recovered from trac.)
    YAML is implemented.           [14 Nov 2012 Y.Namekawa]
    Multi-threaded.                [17 Jul 2014 Y.Namekawa]
    Introduce unique_ptr to avoid memory leaks.
                                   [21 Mar 2015 Y.Namekawa]
    Add restart.                   [22 Feb 2016 Y.Namekawa]
    Add flop_count.                [ 8 Aug 2016 Y.Namekawa]
 */

class Solver_CG : public Solver
{
 public:
  static const std::string class_name;

 private:
  Fopr *m_fopr;

  int    m_Niter;
  int    m_Nrestart;
  double m_Stop_cond;

  //- working area
  Field m_s, m_r, m_x, m_p;

  int m_Nrestart_count;
  int m_Nconv_count;

 public:
  Solver_CG(Fopr *fopr)
    : Solver(), m_fopr(fopr)
  {
    m_Nrestart_count = 0;
    m_Nconv_count    = 0;
  }

  Solver_CG(unique_ptr<Fopr>& fopr)
    : Solver(), m_fopr(fopr.get())
  {
    m_Nrestart_count = 0;
    m_Nconv_count    = 0;
  }

  ~Solver_CG() {}

  void set_parameters(const Parameters& params);
  void set_parameters(const int Niter, const int Nrestart, const double Stop_cond);

  void solve(Field& solution, const Field& source,
             int& Nconv, double& diff);

  Fopr *get_fopr() { return m_fopr; }

  double flop_count();

 private:
  void reset_field(const Field&);

  void solve_init(const Field&, double&);
  void solve_step(double&);
};
#endif
