/*!
        @file    $Id:: solver_GMRES_m_Cmplx.h #$

        @brief

        @author  Yusuke Namekawa (namekawa)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 18:35:38 #$

        @version $LastChangedRevision: 1571 $
*/

#ifndef SOLVER_GMRES_m_CMPLX_INCLUDED
#define SOLVER_GMRES_m_CMPLX_INCLUDED

#include "solver.h"

#include <valarray>

#include "IO/bridgeIO.h"
using Bridge::vout;


//! GMRES(m) algorithm with complex variables.

/*!
   This class implements GMRES(m) algorithm for nonhermitian matrix.
   The product of vectors is treated in complex.
   See Y.Saad and M.H.Schultz, SIAM J.Sci.Stat.Comput. 7 (1986) 856.
                                   [8 Aug 2012  Y.Namekawa]
    (Coding history will be recovered from trac.)
    YAML is implemented.           [14 Nov 2012 Y.Namekawa]
    Multi-threaded.                [17 Jul 2014 Y.Namekawa]
    Introduce unique_ptr to avoid memory leaks.
                                   [21 Mar 2015 Y.Namekawa]
    Add restart.                   [22 Feb 2016 Y.Namekawa]
    Add flop_count.                [ 8 Aug 2016 Y.Namekawa]
 */

class Solver_GMRES_m_Cmplx : public Solver
{
 public:
  static const std::string class_name;

 private:
  Fopr *m_fopr;

  int    m_Niter;
  int    m_Nrestart;
  double m_Stop_cond;

  int m_N_M;

  //- working area
  double m_beta_prev;

  std::vector<Field> m_v;
  Field              m_s, m_r, m_x, m_v_tmp;

  int m_Nrestart_count;
  int m_Nconv_count;

 public:
  Solver_GMRES_m_Cmplx(Fopr *fopr)
    : Solver(), m_fopr(fopr)
  {
    m_Nrestart_count = 0;
    m_Nconv_count    = 0;
  }

  Solver_GMRES_m_Cmplx(unique_ptr<Fopr>& fopr)
    : Solver(), m_fopr(fopr.get())
  {
    m_Nrestart_count = 0;
    m_Nconv_count    = 0;
  }

  ~Solver_GMRES_m_Cmplx() {}

  void set_parameters(const Parameters& params);
  void set_parameters(const int Niter, const int Nrestart, const double Stop_cond);
  void set_parameters_GMRES_m(const int N_M);

  void solve(Field& solution, const Field& source,
             int& Nconv, double& diff);

  Fopr *get_fopr() { return m_fopr; }

  double flop_count();

 private:
  void reset_field(const Field&);

  void solve_init(const Field&, double&);
  void solve_step(const Field&, double&);

  void innerprod_c(double& prod_r, double& prod_i,
                   const Field& v, const Field& w);

  void min_J(std::valarray<dcomplex>& y,
             std::valarray<dcomplex>& h);

  int index_ij(int i, int j)
  {
    return i + (m_N_M + 1) * j;
  }
};
#endif
