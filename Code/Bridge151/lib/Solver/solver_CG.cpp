/*!
        @file    solver_CG.cpp

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/
#include "BridgeLib_Private.h"
#include "solver_CG.h"

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init = Solver_CG::register_factory();
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
  bool   use_init_guess;

  int err = 0;
  err += params.fetch_int("maximum_number_of_iteration", Niter);
  err += params.fetch_int("maximum_number_of_restart", Nrestart);
  err += params.fetch_double("convergence_criterion_squared", Stop_cond);
  err += params.fetch_bool("use_initial_guess", use_init_guess);

  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameter not found.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }

  set_parameters(Niter, Nrestart, Stop_cond, use_init_guess);
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
void Solver_CG::set_parameters(const int Niter, const int Nrestart, const double Stop_cond, const bool use_init_guess)
{
  ThreadManager_OpenMP::assert_single_thread(class_name);

  //- print input parameters
  vout.general(m_vl, "%s:\n", class_name.c_str());
  vout.general(m_vl, "  Niter          = %d\n", Niter);
  vout.general(m_vl, "  Nrestart       = %d\n", Nrestart);
  vout.general(m_vl, "  Stop_cond      = %8.2e\n", Stop_cond);
  vout.general(m_vl, "  use_init_guess = %s\n", use_init_guess ? "true" : "false");

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
  m_Niter          = Niter;
  m_Nrestart       = Nrestart;
  m_Stop_cond      = Stop_cond;
  m_use_init_guess = use_init_guess;
}


//====================================================================
void Solver_CG::solve(Field& xq, const Field& b,
                      int& Nconv, double& diff)
{
  const double bnorm2 = b.norm2();
  const int    bsize  = b.size();

  vout.paranoiac(m_vl, "%s: solver starts\n", class_name.c_str());
  vout.paranoiac(m_vl, "  norm of b = %16.8e\n", bnorm2);
  vout.paranoiac(m_vl, "  size of b = %d\n", bsize);

  bool   is_converged = false;
  int    Nconv2       = 0;
  double diff2        = 1.0;  // superficial initialization
  double rr;

  int Nconv_unit = 1;
  // if (m_fopr->get_mode() == "DdagD" || m_fopr->get_mode() == "DDdag") {
  //   Nconv_unit = 2;
  // }

  reset_field(b);

  if (m_use_init_guess) {
    copy(m_s, xq);  // s = xq;
  } else {
    copy(m_s, b);   // s = b;
  }
  solve_init(b, rr);
  Nconv2 += Nconv_unit;

  vout.detailed(m_vl, "    iter: %8d  %22.15e\n", Nconv2, rr / bnorm2);


  for (int i_restart = 0; i_restart < m_Nrestart; i_restart++) {
    for (int iter = 0; iter < m_Niter; iter++) {
      if (rr / bnorm2 < m_Stop_cond) break;

      solve_step(rr);
      Nconv2 += Nconv_unit;

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

      m_Nrestart_count = i_restart;
      m_Nconv_count    = Nconv2;

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
    const int Nin  = b.nin();
    const int Nvol = b.nvol();
    const int Nex  = b.nex();

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
  const double rr_prev = rr;

  m_fopr->mult(m_s, m_p);           // s = m_fopr->mult(p);

  const double pap = dot(m_p, m_s); // pap  = p * s;
  const double cr  = rr_prev / pap;

  axpy(m_x, cr, m_p);   // x += cr * p;
  axpy(m_r, -cr, m_s);  // r -= cr * s;

  rr = m_r.norm2();     // rr = r * r;

  const double rr_ratio = rr / rr_prev;
  aypx(rr_ratio, m_p, m_r);  // p  = (rr / rr_prev) * p + r
}


//====================================================================
double Solver_CG::flop_count()
{
  const int NPE = CommonParameters::NPE();

  //- NB1 Nin = 2 * Nc * Nd, Nex = 1  for field_F
  //- NB2 Nvol = CommonParameters::Nvol()/2 for eo
  const int Nin    = m_x.nin();
  const int Nvol   = m_x.nvol();
  const int Nex    = m_x.nex();
  const int e_type = m_x.field_element_type();


  const double gflop_fopr = m_fopr->flop_count();

  if (gflop_fopr < CommonParameters::epsilon_criterion()) {
    vout.crucial(m_vl, "Warning at %s: no fopr->flop_count() is available, setting flop = 0\n", class_name.c_str());
    return 0.0;
  }

  const double gflop_axpy = Nin * Nex * 2 * ((Nvol * NPE) / 1.0e+9);
  const double gflop_norm = Nin * Nex * 2 * ((Nvol * NPE) / 1.0e+9);

  double gflop_dot = 0.0;                       // superficial initialization
  if (e_type == Element_type::REAL) {           // element_type REAL = 1
    gflop_dot = Nin * Nex * 2 * ((Nvol * NPE) / 1.0e+9);
  } else if (e_type == Element_type::COMPLEX) { // element_type COMPLEX = 2
    gflop_dot = Nin * Nex * 4 * ((Nvol * NPE) / 1.0e+9);
  } // NB. The other e_type is not allowed by construction in field.h

  const double gflop_init          = gflop_fopr + gflop_axpy + gflop_norm;
  const double gflop_step          = gflop_fopr + gflop_dot + 3 * gflop_axpy + gflop_norm;
  const double gflop_true_residual = gflop_fopr + gflop_axpy + gflop_norm;

  const double gflop = gflop_norm + gflop_init
                       + gflop_step * (m_Nconv_count - 1)
                       + gflop_true_residual * (m_Nrestart_count + 1)
                       + gflop_init * m_Nrestart_count;

  return gflop;
}


//====================================================================
//============================================================END=====
