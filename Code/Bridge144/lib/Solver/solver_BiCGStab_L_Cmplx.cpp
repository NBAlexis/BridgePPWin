#include "BridgeLib_Private.h"

/*!
        @file    $Id:: solver_BiCGStab_L_Cmplx.cpp #$

        @brief

        @author  Yusuke Namekawa (namekawa)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 18:35:38 #$

        @version $LastChangedRevision: 1571 $
*/

#include "solver_BiCGStab_L_Cmplx.h"


#ifdef USE_FACTORY
namespace {
  Solver *create_object(Fopr *fopr)
  {
    return new Solver_BiCGStab_L_Cmplx(fopr);
  }


  bool init = Solver::Factory::Register("BiCGStab_L_Cmplx", create_object);
}
#endif


const std::string Solver_BiCGStab_L_Cmplx::class_name = "Solver_BiCGStab_L_Cmplx";

//====================================================================
void Solver_BiCGStab_L_Cmplx::set_parameters(const Parameters& params)
{
  const string str_vlevel = params.get_string("verbose_level");

  m_vl = vout.set_verbose_level(str_vlevel);

  //- fetch and check input parameters
  int    Niter, Nrestart;
  double Stop_cond;
  int    N_L;

  int err = 0;
  err += params.fetch_int("maximum_number_of_iteration", Niter);
  err += params.fetch_int("maximum_number_of_restart", Nrestart);
  err += params.fetch_double("convergence_criterion_squared", Stop_cond);
  err += params.fetch_int("number_of_orthonormal_vectors", N_L);

  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameter not found.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }

  set_parameters(Niter, Nrestart, Stop_cond);
  set_parameters_L(N_L);
}


//====================================================================
void Solver_BiCGStab_L_Cmplx::set_parameters(const int Niter, const int Nrestart, const double Stop_cond)
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
void Solver_BiCGStab_L_Cmplx::set_parameters_L(const int N_L)
{
  //- print input parameters
  vout.general(m_vl, "  N_L   = %d\n", N_L);

  //- range check
  int err = 0;
  err += ParameterCheck::non_negative(N_L);

  if (err) {
    vout.crucial(m_vl, "Error at %s: parameter range check failed.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }

  //- store values
  m_N_L = N_L;
}


//====================================================================
void Solver_BiCGStab_L_Cmplx::solve(Field& xq, const Field& b,
                                    int& Nconv, double& diff)
{
  double bnorm2 = b.norm2();
  int    bsize  = b.size();

  vout.paranoiac(m_vl, "%s: starts\n", class_name.c_str());
  vout.paranoiac(m_vl, "  norm of b = %16.8e\n", bnorm2);
  vout.paranoiac(m_vl, "  size of b = %d\n", bsize);

  bool   is_converged = false;
  int    Nconv2       = 0;
  double diff2        = 1.0;
  double rr;

  reset_field(b);
  copy(m_s, b); // s = b;
  solve_init(b, rr);
  Nconv2 += 1;

  vout.detailed(m_vl, "    iter: %8d  %22.15e\n", Nconv2, rr / bnorm2);


  for (int i_restart = 0; i_restart < m_Nrestart; i_restart++) {
    for (int iter = 0; iter < m_Niter; iter++) {
      if (rr / bnorm2 < m_Stop_cond) break;

      solve_step(rr);
      Nconv2 += 2 * m_N_L;

      vout.detailed(m_vl, "    iter: %8d  %22.15e\n", Nconv2, rr / bnorm2);
    }

    //- calculate true residual
    m_fopr->mult(m_s, m_x); // s  = m_fopr->mult(x);
    axpy(m_s, -1.0, b);     // s -= b;
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
void Solver_BiCGStab_L_Cmplx::reset_field(const Field& b)
{
#pragma omp barrier
#pragma omp master
  {
    int Nin  = b.nin();
    int Nvol = b.nvol();
    int Nex  = b.nex();

    if ((m_s.nin() != Nin) || (m_s.nvol() != Nvol) || (m_s.nex() != Nex)) {
      m_s.reset(Nin, Nvol, Nex);
      m_x.reset(Nin, Nvol, Nex);
      m_r_init.reset(Nin, Nvol, Nex);
      m_v.reset(Nin, Nvol, Nex);
    }

    m_u.resize(m_N_L + 1);
    m_r.resize(m_N_L + 1);

    for (int i = 0; i < m_N_L + 1; ++i) {
      m_u[i].reset(Nin, Nvol, Nex);
      m_r[i].reset(Nin, Nvol, Nex);
    }
  }
#pragma omp barrier

  vout.paranoiac(m_vl, "    %s: field size reset.\n", class_name.c_str());
}


//====================================================================
void Solver_BiCGStab_L_Cmplx::solve_init(const Field& b, double& rr)
{
  copy(m_x, m_s);  // x = s;

  for (int i = 0; i < m_N_L + 1; ++i) {
    m_r[i].set(0.0);        // r[i] = 0.0;
    m_u[i].set(0.0);        // u[i] = 0.0;
  }

  // r[0] = b - A x_0
  m_fopr->mult(m_v, m_s);   // m_v = m_fopr->mult(s);
  copy(m_r[0], b);          // r[0]  = b;
  axpy(m_r[0], -1.0, m_v);  // r[0] -= m_v;

  copy(m_r_init, m_r[0]);   // r_init = r[0];
  rr = m_r[0].norm2();      // rr     = r[0] * r[0];

#pragma omp barrier
#pragma omp master
  {
    m_rho_prev = cmplx(-1.0, 0.0);

    // NB. m_alpha_prev = 0.0 \neq 1.0
    m_alpha_prev = cmplx(0.0, 0.0);
  }
#pragma omp barrier
}


//====================================================================
void Solver_BiCGStab_L_Cmplx::solve_step(double& rr)
{
  dcomplex rho_prev2   = m_rho_prev;
  dcomplex alpha_prev2 = m_alpha_prev;

  for (int j = 0; j < m_N_L; ++j) {
    dcomplex rho = dotc(m_r[j], m_r_init);  // dcomplex rho  = r[j] * r_init;
    rho = conj(rho);

    dcomplex beta = alpha_prev2 * (rho / rho_prev2);

    rho_prev2 = rho;

    for (int i = 0; i < j + 1; ++i) {
      aypx(-beta, m_u[i], m_r[i]);  // u[i] = - beta * u[i] + r[i];
    }

    m_fopr->mult(m_u[j + 1], m_u[j]);  // u[j+1] = m_fopr->mult(u[j]);

    dcomplex gamma = dotc(m_u[j + 1], m_r_init);
    alpha_prev2 = rho_prev2 / conj(gamma);

    for (int i = 0; i < j + 1; ++i) {
      axpy(m_r[i], -alpha_prev2, m_u[i + 1]);  // r[i] -= alpha_prev * u[i+1];
    }

    m_fopr->mult(m_r[j + 1], m_r[j]);  // r[j+1] = m_fopr->mult(r[j]);

    axpy(m_x, alpha_prev2, m_u[0]);    // x += alpha_prev * u[0];
  }


  std::vector<double>   sigma(m_N_L + 1);
  std::vector<dcomplex> gamma_prime(m_N_L + 1);

  // NB. tau(m_N_L,m_N_L+1), not (m_N_L+1,m_N_L+1)
  std::vector<dcomplex> tau(m_N_L * (m_N_L + 1));
  int ij, ji;

  for (int j = 1; j < m_N_L + 1; ++j) {
    for (int i = 1; i < j; ++i) {
      ij = index_ij(i, j);

      dcomplex r_ji = dotc(m_r[j], m_r[i]);
      tau[ij] = conj(r_ji) / sigma[i];  // tau[ij]  = (r[j] * r[i]) / sigma[i];
      axpy(m_r[j], -tau[ij], m_r[i]);   // r[j]    -= tau[ij] * r[i];
    }

    sigma[j] = m_r[j].norm2();  // sigma[j] = r[j] * r[j];

    dcomplex r_0j = dotc(m_r[0], m_r[j]);
    gamma_prime[j] = conj(r_0j) / sigma[j]; // gamma_prime[j] = (r[0] * r[j]) / sigma[j];
  }


  std::vector<dcomplex> gamma(m_N_L + 1);
  dcomplex              c_tmp;

  gamma[m_N_L] = gamma_prime[m_N_L];

  for (int j = m_N_L - 1; j > 0; --j) {
    c_tmp = cmplx(0.0, 0.0);

    for (int i = j + 1; i < m_N_L + 1; ++i) {
      ji     = index_ij(j, i);
      c_tmp += tau[ji] * gamma[i];
    }

    gamma[j] = gamma_prime[j] - c_tmp;
  }


  // NB. gamma_double_prime(m_N_L), not (m_N_L+1)
  std::vector<dcomplex> gamma_double_prime(m_N_L);

  for (int j = 1; j < m_N_L; ++j) {
    c_tmp = cmplx(0.0, 0.0);

    for (int i = j + 1; i < m_N_L; ++i) {
      ji     = index_ij(j, i);
      c_tmp += tau[ji] * gamma[i + 1];
    }

    gamma_double_prime[j] = gamma[j + 1] + c_tmp;
  }


  axpy(m_x, gamma[1], m_r[0]);                    // x    += gamma[          1] * r[    0];
  axpy(m_r[0], -gamma_prime[m_N_L], m_r[m_N_L]);  // r[0] -= gamma_prime[m_N_L] * r[m_N_L];
  axpy(m_u[0], -gamma[m_N_L], m_u[m_N_L]);        // u[0] -= gamma[      m_N_L] * u[m_N_L];

  for (int j = 1; j < m_N_L; ++j) {
    axpy(m_x, gamma_double_prime[j], m_r[j]);     // x    += gamma_double_prime[j] * r[j];
    axpy(m_r[0], -gamma_prime[j], m_r[j]);        // r[0] -= gamma_prime[       j] * r[j];
    axpy(m_u[0], -gamma[j], m_u[j]);              // u[0] -= gamma[             j] * u[j];
  }

  rr = m_r[0].norm2();  // rr = r[0] * r[0];

#pragma omp barrier
#pragma omp master
  {
    m_rho_prev   = rho_prev2;
    m_alpha_prev = alpha_prev2;

    m_rho_prev *= -gamma_prime[m_N_L];
  }
#pragma omp barrier
}


//====================================================================
double Solver_BiCGStab_L_Cmplx::flop_count()
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

  int N_iter = (m_Nconv_count - 1) / (2 * m_N_L);

  int N_L_part = 0;
  for (int j = 0; j < m_N_L; ++j) {
    for (int i = 0; i < j + 1; ++i) {
      N_L_part += 1;
    }
  }

  double flop_init           = flop_fopr + flop_axpy + flop_norm;
  double flop_step_BiCG_part = 2 * m_N_L * flop_fopr
                               + 2 * m_N_L * flop_dotc
                               + (m_N_L + 2 * N_L_part) * flop_axpy;
  double flop_step_L_part = (N_L_part + m_N_L) * flop_dotc
                            + (N_L_part + 3 * m_N_L) * flop_axpy
                            + (m_N_L + 1) * flop_norm;
  double flop_step          = flop_step_BiCG_part + flop_step_L_part;
  double flop_true_residual = flop_fopr + flop_axpy + flop_norm;

  double flop = flop_norm + flop_init + flop_step * N_iter + flop_true_residual
                + flop_init * m_Nrestart_count;


  return flop;
}


//====================================================================
//============================================================END=====
