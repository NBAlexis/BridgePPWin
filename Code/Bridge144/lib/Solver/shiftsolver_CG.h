/*!
        @file    $Id:: shiftsolver_CG.h #$

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 18:35:38 #$

        @version $LastChangedRevision: 1571 $
*/

#ifndef SHIFTSOLVER_CG_INCLUDED
#define SHIFTSOLVER_CG_INCLUDED

#include "shiftsolver.h"

#include "IO/bridgeIO.h"
using Bridge::vout;


//! Multishift Conjugate Gradient solver.

/*!
                                    [23 Dec 2011  H.Matsufuru]
    unique_ptr is introduced to avoid memory leaks
                                    [21 Mar 2015 Y.Namekawa]
 */

class Shiftsolver_CG : public Shiftsolver
{
 public:
  static const std::string class_name;

 private:
  Fopr *m_fopr;

  int    m_Niter;
  double m_Stop_cond;

  std::vector<Field>  m_x, m_p;
  Field               m_r, m_s;
  std::vector<double> m_zeta1, m_zeta2, m_csh2, m_pp;

  double m_snorm, m_alpha_p, m_beta_p;
  int    m_Nshift2;

  double m_sigma0;

 public:
  Shiftsolver_CG(Fopr *fopr)
    : Shiftsolver(), m_fopr(fopr) {}

  Shiftsolver_CG(Fopr *fopr, int niter, double stop_cond)
    : Shiftsolver(), m_fopr(fopr)
  {
    set_parameters(niter, stop_cond);
  }

  Shiftsolver_CG(unique_ptr<Fopr>& fopr)
    : Shiftsolver(), m_fopr(fopr.get()) {}

  Shiftsolver_CG(unique_ptr<Fopr>& fopr, int niter, double stop_cond)
    : Shiftsolver(), m_fopr(fopr.get())
  {
    set_parameters(niter, stop_cond);
  }

  ~Shiftsolver_CG() {}

  void set_parameters(const Parameters& params);
  void set_parameters(const int niter, const double stop_cond);

  void solve(
    std::vector<Field>& solution,
    const std::vector<double>& shift,
    const Field& source,
    int& Nconv, double& diff);

 private:

  void solve_init(double&);
  void solve_step(double&);

  void reset_field(const Field& b, const std::vector<double>& sigma, const int Nshift);
};
#endif
