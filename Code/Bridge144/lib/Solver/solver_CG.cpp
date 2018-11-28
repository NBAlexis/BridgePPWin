#include "BridgeLib_Private.h"

/*!
        @file    $Id:: solver_CG.cpp #$

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 18:35:38 #$

        @version $LastChangedRevision: 1571 $
*/

#include "solver_CG.h"


#ifdef USE_FACTORY
namespace {
  Solver *create_object(Fopr *fopr)
  {
    return new Solver_CG(fopr);
  }


  bool init = Solver::Factory::Register("CG", create_object);
}
#endif


const std::string Solver_CG::class_name = "Solver_CG";

//====================================================================
void Solver_CG::set_parameters(const Parameters& params)
{
  const string str_vlevel = params.get_string("verbose_level");

  m_vl = vout.set_verbose_level(str_vlevel);

  //- fetch and check input parameters
  int    Niter, Nrestart;
  double Stop_cond;

  int err = 0;
  err += params.fetch_int("maximum_number_of_iteration", Niter);
  err += params.fetch_int("maximum_number_of_restart", Nrestart);
  err += params.fetch_double("convergence_criterion_squared", Stop_cond);

  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameter not found.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }

  set_parameters(Niter, Nrestart, Stop_cond);
}


//====================================================================
void Solver_CG::set_parameters(const int Niter, const int Nrestart, const double Stop_cond)
{
  ThreadManager_OpenMP::assert_single_thread(class_name);

  //- print input parameters
  vout.general(m_vl, "%s:\n", class_name.c_str());
  vout.general(m_vl, "  Niter     = %d\n", Niter);
  vout.general(m_vl, "  Nrestart  = %d\n", Nrestart);
  vout.general(m_vl, "  Stop_cond = %8.2e\n", Stop_cond);

  //- range check
  int err = 0;
  err += ParameterCheck::non_negative(Niter);
  err += ParameterCheck::non_negative(Nrestart);
  err += ParameterCheck::square_non_zero(Stop_cond);

  if (err) {
    vout.crucial(m_vl, "Error at %s: parameter range check failed.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }

  //- store values
  m_Niter     = Niter;
  m_Nrestart  = Nrestart;
  m_Stop_cond = Stop_cond;
}


//====================================================================
void Solver_CG::solve(Field& xq, const Field& b,
                      int& Nconv, double& diff)
{
  double bnorm2 = b.norm2();
  int    bsize  = b.size();

  vout.paranoiac(m_vl, "%s: solver starts\n", class_name.c_str());
  vout.paranoiac(m_vl, "  norm of b = %16.8e\n", bnorm2);
  vout.paranoiac(m_vl, "  size of b = %d\n", bsize);

  bool   is_converged = false;
  int    Nconv2       = 0;
  double diff2        = 1.0;
  double rr;

  reset_field(b);
  copy(m_s, b);  // s = b;
  solve_init(b, rr);
  Nconv2 += 1;

  vout.detailed(m_vl, "    iter: %8d  %22.15e\n", Nconv2, rr / bnorm2);


  for (int i_restart = 0; i_restart < m_Nrestart; i_restart++) {
    for (int iter = 0; iter < m_Niter; iter++) {
      if (rr / bnorm2 < m_Stop_cond) break;

      solve_step(rr);
      Nconv2 += 1;

      vout.detailed(m_vl, "    iter: %8d  %22.15e\n", Nconv2, rr / bnorm2);
    }

    //- calculate true residual
    m_fopr->mult(m_s, m_x);  // s  = m_fopr->mult(x);
    axpy(m_s, -1.0, b);      // s -= b;
    diff2 = m_s.norm2();

    if (diff2 / bnorm2 < m_Stop_cond) {
      vout.detailed(m_vl, "%s: converged.\n", class_name.c_str());
      vout.detailed(m_vl, "  iter(final): %8d  %22.15e\n", Nconv2, diff2 / bnorm2);

      is_converged = true;
      break;
    } else {
      //- restart with new approximate solution
      copy(m_s, m_x); // s = x;
      solve_init(b, rr);

      vout.detailed(m_vl, "%s: restarted.\n", class_name.c_str());
    }
  }


  if (!is_converged) {
    vout.crucial(m_vl, "Error at %s: not converged.\n", class_name.c_str());
    vout.crucial(m_vl, "  iter(final): %8d  %22.15e\n", Nconv2, diff2 / bnorm2);
    exit(EXIT_FAILURE);
  }


  copy(xq, m_x); // xq = x;

#pragma omp barrier
#pragma omp master
  {
    diff  = sqrt(diff2 / bnorm2);
    Nconv = Nconv2;
  }
#pragma omp barrier
}


//====================================================================
void Solver_CG::reset_field(const Field& b)
{
#pragma omp barrier
#pragma omp master
  {
    int Nin  = b.nin();
    int Nvol = b.nvol();
    int Nex  = b.nex();

    if ((m_s.nin() != Nin) || (m_s.nvol() != Nvol) || (m_s.nex() != Nex)) {
      m_s.reset(Nin, Nvol, Nex);
      m_r.reset(Nin, Nvol, Nex);
      m_x.reset(Nin, Nvol, Nex);
      m_p.reset(Nin, Nvol, Nex);
    }
  }
#pragma omp barrier

  vout.detailed(m_vl, "    %s: field size reset.\n", class_name.c_str());
}


//====================================================================
void Solver_CG::solve_init(const Field& b, double& rr)
{
  copy(m_x, m_s);  // x = s;

  // r = b - A x
  copy(m_r, b);            // r = b;
  m_fopr->mult(m_s, m_x);  // s  = m_fopr->mult(x);
  axpy(m_r, -1.0, m_s);    // r -= s;

  copy(m_p, m_r);          // p  = r;
  rr = m_r.norm2();        // rr = r * r;
}


//====================================================================
void Solver_CG::solve_step(double& rr)
{
  double rr_prev = rr;

  m_fopr->mult(m_s, m_p);     // s = m_fopr->mult(p);

  double pap = dot(m_p, m_s); // pap  = p * s;
  double cr  = rr_prev / pap;

  axpy(m_x, cr, m_p);   // x += cr * p;
  axpy(m_r, -cr, m_s);  // r -= cr * s;

  rr = m_r.norm2();     // rr = r * r;

  double rr_ratio = rr / rr_prev;
  aypx(rr_ratio, m_p, m_r);  // p  = (rr / rr_prev) * p + r
}


//====================================================================
double Solver_CG::flop_count()
{
  int    NPE = CommonParameters::NPE();
  double eps = CommonParameters::epsilon_criterion();

  //- NB1 Nin = 2 * Nc * Nd, Nex = 1  for field_F
  //- NB2 Nvol = CommonParameters::Nvol()/2 for eo
  int Nin  = m_x.nin();
  int Nvol = m_x.nvol();
  int Nex  = m_x.nex();

  double flop_fopr = m_fopr->flop_count() / (Nvol * NPE);

  if (flop_fopr < eps) {
    vout.crucial(m_vl, "Warning at %s: no fopr->flop_count() is available, setting flop = 0.0.\n", class_name.c_str());
    return 0.0;
  }

  double flop_axpy = static_cast<double>(Nin * Nex * 2);
  double flop_dot  = static_cast<double>(Nin * Nex * 2);  // (Nin * Nex * 4) for Cmplx
  double flop_norm = static_cast<double>(Nin * Nex * 2);

  double flop_init          = flop_fopr + flop_axpy + flop_norm;
  double flop_step          = flop_fopr + flop_dot + 3 * flop_axpy + flop_norm;
  double flop_true_residual = flop_fopr + flop_axpy + flop_norm;

  double flop = (flop_norm + flop_init + flop_step * m_Nconv_count + flop_true_residual
                 + flop_init * m_Nrestart_count) * (Nvol * NPE);


  return flop;
}


//====================================================================
//============================================================END=====
