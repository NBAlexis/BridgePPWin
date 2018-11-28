/*!
        @file    $Id:: solver_BiCGStab_Cmplx.h #$

        @brief

        @author  Yusuke Namekawa (namekawa)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 18:35:38 #$

        @version $LastChangedRevision: 1571 $
*/

#ifndef SOLVER_BICGSTAB_CMPLX_INCLUDED
#define SOLVER_BICGSTAB_CMPLX_INCLUDED

#include "solver.h"

#include "IO/bridgeIO.h"
using Bridge::vout;



//! BiCGStab algorithm with complex variables.

/*!
   This class implements BiCGStab algorithm for nonhermitian matrix.
   The product of vectors is treated in complex.
                                   [12 Feb 2012 Y.Namekawa]
    (Coding history will be recovered from trac.)
    YAML is implemented.           [14 Nov 2012 Y.Namekawa]
    Multi-threaded.                [10 Jul 2014 H.Matsufuru]
    Introduce unique_ptr to avoid memory leaks.
                                   [21 Mar 2015 Y.Namekawa]
    Add restart.                   [22 Feb 2016 Y.Namekawa]
    Add flop_count.                [ 8 Aug 2016 Y.Namekawa]
 */

class Solver_BiCGStab_Cmplx : public Solver
{
 public:
  static const std::string class_name;

 private:
  Fopr *m_fopr;

  int    m_Niter, m_Nrestart;
  double m_Stop_cond;

  //- working area
  dcomplex m_rho_prev, m_alpha_prev, m_omega_prev;
  Field    m_s, m_r, m_x, m_rh, m_p, m_v, m_t;

  int m_Nrestart_count;
  int m_Nconv_count;

 public:
  Solver_BiCGStab_Cmplx(Fopr *fopr)
    : Solver(), m_fopr(fopr)
  {
    m_Nrestart_count = 0;
    m_Nconv_count    = 0;
  }

  Solver_BiCGStab_Cmplx(unique_ptr<Fopr>& fopr)
    : Solver(), m_fopr(fopr.get())
  {
    m_Nrestart_count = 0;
    m_Nconv_count    = 0;
  }

  ~Solver_BiCGStab_Cmplx() {}

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
