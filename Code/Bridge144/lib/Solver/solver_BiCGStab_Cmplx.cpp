#include "BridgeLib_Private.h"

/*!
        @file    $Id:: solver_BiCGStab_Cmplx.cpp #$

        @brief

        @author  Yusuke Namekawa (namekawa)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 18:35:38 #$

        @version $LastChangedRevision: 1571 $
*/

#include "solver_BiCGStab_Cmplx.h"


#ifdef USE_FACTORY
namespace {
  Solver *create_object(Fopr *fopr)
  {
    return new Solver_BiCGStab_Cmplx(fopr);
  }


  bool init = Solver::Factory::Register("BiCGStab_Cmplx", create_object);
}
#endif


const std::string Solver_BiCGStab_Cmplx::class_name = "Solver_BiCGStab_Cmplx";

//====================================================================
void Solver_BiCGStab_Cmplx::set_parameters(const Parameters& params)
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
    vout.crucial(m_vl, "Error at %s: input parameter not found.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

  set_parameters(Niter, Nrestart, Stop_cond);
}


//====================================================================
void Solver_BiCGStab_Cmplx::set_parameters(const int Niter, const int Nrestart,
                                           const double Stop_cond)
{
  ThreadManager_OpenMP::assert_single_thread(class_name);

  //- print input parameters
  vout.general(m_vl, "%s: input parameters\n", class_name.c_str());
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
void Solver_BiCGStab_Cmplx::solve(Field& xq, const Field& b,
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
      Nconv2 += 2;

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
      copy(m_s, m_x);  // s = x;
      solve_init(b, rr);

      vout.detailed(m_vl, "%s: restarted.\n", class_name.c_str());
    }
  }


  if (!is_converged) {
    vout.crucial(m_vl, "Error at %s: not converged.\n", class_name.c_str());
    vout.crucial(m_vl, "  iter(final): %8d  %22.15e\n", Nconv2, diff2 / bnorm2);
    exit(EXIT_FAILURE);
  }


  copy(xq, m_x);  // xq = x;

#pragma omp barrier
#pragma omp master
  {
    diff  = sqrt(diff2 / bnorm2);
    Nconv = Nconv2;
  }
#pragma omp barrier
}


//====================================================================
void Solver_BiCGStab_Cmplx::reset_field(const Field& b)
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
      m_v.reset(Nin, Nvol, Nex);
      m_t.reset(Nin, Nvol, Nex);
      m_rh.reset(Nin, Nvol, Nex);
    }
  }
#pragma omp barrier

  vout.detailed(m_vl, "    %s: field size reset.\n", class_name.c_str());
}


//====================================================================
void Solver_BiCGStab_Cmplx::solve_init(const Field& b, double& rr)
{
  copy(m_x, m_s);  // x = s;

  // r = b - A x_0
  m_fopr->mult(m_v, m_s); // v = m_fopr->mult(s);
  copy(m_r, b);           // r  = b;
  axpy(m_r, -1.0, m_v);   // r -= v;
  copy(m_rh, m_r);        // rh = r;

  rr = m_r.norm2();       // rr = r * r;

  m_p.set(0.0);           // p = 0.0
  m_v.set(0.0);           // v = 0.0

#pragma omp barrier
#pragma omp master
  {
    m_rho_prev   = cmplx(1.0, 0.0);
    m_alpha_prev = cmplx(1.0, 0.0);
    m_omega_prev = cmplx(1.0, 0.0);
  }
#pragma omp barrier
}


//====================================================================
void Solver_BiCGStab_Cmplx::solve_step(double& rr)
{
  dcomplex rho  = dotc(m_rh, m_r);  // rho = rh * r;
  dcomplex beta = rho * m_alpha_prev / (m_rho_prev * m_omega_prev);

  // p = r + beta * (p - m_omega_prev * v);
  axpy(m_p, -m_omega_prev, m_v);     // p += - m_omega_prev * v;
  aypx(beta, m_p, m_r);              // p  = beta * p + r;

  m_fopr->mult(m_v, m_p);            // v  = m_fopr->mult(p);

  dcomplex aden  = dotc(m_rh, m_v);  // aden = rh * v;
  dcomplex alpha = rho / aden;

  copy(m_s, m_r);                    // s  = r
  axpy(m_s, -alpha, m_v);            // s += - alpha * v;

  m_fopr->mult(m_t, m_s);            // t  = m_fopr->mult(s);

  double   omega_d = dot(m_t, m_t);  // omega_d = t * t;
  dcomplex omega_n = dotc(m_t, m_s); // omega_n = t * s;
  dcomplex omega   = omega_n / omega_d;

  axpy(m_x, omega, m_s);   // x += omega * s;
  axpy(m_x, alpha, m_p);   // x += alpha * p;

  copy(m_r, m_s);          // r  = s
  axpy(m_r, -omega, m_t);  // r += - omega * t;

  rr = m_r.norm2();        // rr = r * r;

#pragma omp barrier
#pragma omp master
  {
    m_rho_prev   = rho;
    m_alpha_prev = alpha;
    m_omega_prev = omega;
  }
#pragma omp barrier
}


//====================================================================
double Solver_BiCGStab_Cmplx::flop_count()
{
  int    NPE = CommonParameters::NPE();
  double eps = CommonParameters::epsilon_criterion();

  //- NB1 Nin = 2 * Nc * Nd, Nex = 1  for field_F
  //- NB2 Nvol = CommonParameters::Nvol()/2 for eo
  int Nin  = m_x.nin();
  int Nvol = m_x.nvol();
  int Nex  = m_x.nex();

  double flop_fopr = m_fopr->flop_count();

  if (flop_fopr < eps) {
    vout.crucial(m_vl, "Warning at %s: no fopr->flop_count() is available, setting flop = 0.0.\n", class_name.c_str());
    return 0.0;
  }

  double flop_axpy = static_cast<double>(Nin * Nex * 2) * (Nvol * NPE);
  double flop_dotc = static_cast<double>(Nin * Nex * 4) * (Nvol * NPE);
  double flop_norm = static_cast<double>(Nin * Nex * 2) * (Nvol * NPE);

  int N_iter = (m_Nconv_count - 1) / 2;

  double flop_init          = flop_fopr + flop_axpy + flop_norm;
  double flop_step          = 2 * flop_fopr + 3 * flop_dotc + 6 * flop_axpy + 2 * flop_norm;
  double flop_true_residual = flop_fopr + flop_axpy + flop_norm;

  double flop = flop_norm + flop_init + flop_step * N_iter + flop_true_residual
                + flop_init * m_Nrestart_count;


  return flop;
}


//====================================================================
//============================================================END=====
