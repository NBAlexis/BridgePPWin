/*!
        @file    $Id:: solver_BiCGStab_L_Cmplx.h #$

        @brief

        @author  Yusuke Namekawa (namekawa)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 18:35:38 #$

        @version $LastChangedRevision: 1571 $
*/

#ifndef SOLVER_BICGSTAB_L_CMPLX_INCLUDED
#define SOLVER_BICGSTAB_L_CMPLX_INCLUDED

#include "solver.h"

#include "IO/bridgeIO.h"
using Bridge::vout;


//! BiCGStab(L) algorithm.

/*!
    This class implements BiCGStab(L) algorithm for nonhermitian matrix.
    The product of vectors is treated in complex.
    See G.L.G.Sleijpen and D.R.Fokkema, Elec.Trans.Numer.Anal.
    1 (1993) 11.
                                   [22 Jan 2012 Y.Namekawa]
    (Coding history will be recovered from trac.)
    YAML is implemented.           [14 Nov 2012 Y.Namekawa]
    Multi-threaded.                [17 Jul 2014 Y.Namekawa]
    Introduce unique_ptr to avoid memory leaks.
                                   [21 Mar 2015 Y.Namekawa]
    Add restart.                   [22 Feb 2016 Y.Namekawa]
    Add flop_count.                [ 8 Aug 2016 Y.Namekawa]
 */

class Solver_BiCGStab_L_Cmplx : public Solver
{
 public:
  static const std::string class_name;

 private:
  Fopr *m_fopr;

  int    m_Niter, m_Nrestart;
  double m_Stop_cond;

  int    m_N_L;
  double m_Tol_L;

  //- working area
  dcomplex m_rho_prev, m_alpha_prev;

  std::vector<Field> m_u, m_r;
  Field              m_s, m_x, m_r_init, m_v;

  int m_Nrestart_count;
  int m_Nconv_count;

 public:
  Solver_BiCGStab_L_Cmplx(Fopr *fopr)
    : Solver(), m_fopr(fopr)
  {
    m_Nrestart_count = 0;
    m_Nconv_count    = 0;
  }

  Solver_BiCGStab_L_Cmplx(unique_ptr<Fopr>& fopr)
    : Solver(), m_fopr(fopr.get())
  {
    m_Nrestart_count = 0;
    m_Nconv_count    = 0;
  }

  ~Solver_BiCGStab_L_Cmplx() {}

  void set_parameters(const Parameters& params);
  void set_parameters(const int Niter, const int Nrestart, const double Stop_cond);
  void set_parameters_L(const int N_L);

  void solve(Field& solution, const Field& source,
             int& Nconv, double& diff);

  Fopr *get_fopr() { return m_fopr; }

  double flop_count();

 private:
  void reset_field(const Field&);

  void solve_init(const Field&, double&);
  void solve_step(double&);

  int index_ij(int i, int j)
  {
    return i + m_N_L * j;
  }
};
#endif
