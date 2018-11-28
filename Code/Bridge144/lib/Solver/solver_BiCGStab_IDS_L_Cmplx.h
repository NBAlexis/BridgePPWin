/*!
        @file    $Id:: solver_BiCGStab_IDS_L_Cmplx.h #$

        @brief

        @author  Yusuke Namekawa (namekawa)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 18:35:38 #$

        @version $LastChangedRevision: 1571 $
*/

#ifndef SOLVER_BICGSTAB_IDS_L_CMPLX_INCLUDED
#define SOLVER_BICGSTAB_IDS_L_CMPLX_INCLUDED

#include "solver.h"

#include "IO/bridgeIO.h"
using Bridge::vout;


//! BiCGStab(IDS_L) algorithm.

/*!
    This class implements BiCGStab(IDS_L) algorithm for nonhermitian matrix.
    The product of vectors is treated in complex.
    See S.Itoh and Y.Namekawa, J.Comp.Appl.Math. 159 (2003) 65.
                                   [22 Jan 2012 Y.Namekawa]
    (Coding history will be recovered from trac.)
    YAML is implemented.           [14 Nov 2012 Y.Namekawa]
    Multi-threaded.                [17 Jul 2014 Y.Namekawa]
    Introduce unique_ptr to avoid memory leaks.
                                   [21 Mar 2015 Y.Namekawa]
    Add restart.                   [22 Feb 2016 Y.Namekawa]
    Add flop_count.                [ 8 Aug 2016 Y.Namekawa]
 */

class Solver_BiCGStab_IDS_L_Cmplx : public Solver
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
  double   m_rr_prev;
  int      m_N_L_prev;

  std::vector<Field> m_u, m_r;
  Field              m_s, m_x, m_r_init, m_v;

  int m_Nrestart_count;
  int m_Niter_count;
  int m_Nconv_count;
  int m_N_L_part_count;

 public:
  Solver_BiCGStab_IDS_L_Cmplx(Fopr *fopr)
    : Solver(), m_fopr(fopr)
  {
    m_Nrestart_count = 0;
    m_Niter_count    = 0;
    m_Nconv_count    = 0;
    m_N_L_part_count = 0;
  }

  Solver_BiCGStab_IDS_L_Cmplx(unique_ptr<Fopr>& fopr)
    : Solver(), m_fopr(fopr.get())
  {
    m_Nrestart_count = 0;
    m_Niter_count    = 0;
    m_Nconv_count    = 0;
    m_N_L_part_count = 0;
  }

  ~Solver_BiCGStab_IDS_L_Cmplx() {}

  void set_parameters(const Parameters& params);
  void set_parameters(const int Niter, const int Nrestart, const double Stop_cond);
  void set_parameters_DS_L(const int N_L, const double Tol_L);

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
