/*!
        @file    solver_GMRES_m_Cmplx.cpp

        @brief

        @author  Yusuke Namekawa (namekawa)
                 $LastChangedBy: namekawa $

        @date    $LastChangedDate:: 2019-03-07 03:01:48 #$

        @version $LastChangedRevision: 1953 $
*/
#include "BridgeLib_Private.h"
#include "solver_GMRES_m_Cmplx.h"

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init = Solver_GMRES_m_Cmplx::register_factory();
}
#endif

const std::string Solver_GMRES_m_Cmplx::class_name = "Solver_GMRES_m_Cmplx";

//====================================================================
void Solver_GMRES_m_Cmplx::set_parameters(const Parameters& params)
{
  const std::string str_vlevel = params.get_string("verbose_level");

  m_vl = vout.set_verbose_level(str_vlevel);

  //- fetch and check input parameters
  int    Niter, Nrestart;
  double Stop_cond;
  bool   use_init_guess;
  int    N_M;

  int err = 0;
  err += params.fetch_int("maximum_number_of_iteration", Niter);
  err += params.fetch_int("maximum_number_of_restart", Nrestart);
  err += params.fetch_double("convergence_criterion_squared", Stop_cond);
  err += params.fetch_bool("use_initial_guess", use_init_guess);
  err += params.fetch_int("number_of_orthonormal_vectors", N_M);

  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameter not found.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }

  set_parameters(Niter, Nrestart, Stop_cond, use_init_guess);
  set_parameters_GMRES_m(N_M);
}


//====================================================================
void Solver_GMRES_m_Cmplx::set_parameters(const int Niter, const int Nrestart, const double Stop_cond)
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
void Solver_GMRES_m_Cmplx::set_parameters(const int Niter, const int Nrestart, const double Stop_cond, const bool use_init_guess)
{
  ThreadManager_OpenMP::assert_single_thread(class_name);

  //- print input parameters
  vout.general(m_vl, "%s: input parameters\n", class_name.c_str());
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
void Solver_GMRES_m_Cmplx::set_parameters_GMRES_m(const int N_M)
{
  //- print input parameters
  vout.general(m_vl, "  N_M   = %d\n", N_M);

  //- range check
  int err = 0;
  err += ParameterCheck::non_negative(N_M);

  if (err) {
    vout.crucial(m_vl, "Error at %s: parameter range check failed.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }

  //- store values
  m_N_M = N_M;
}


//====================================================================
void Solver_GMRES_m_Cmplx::solve(Field& xq, const Field& b,
                                 int& Nconv, double& diff)
{
  const double bnorm2 = b.norm2();
  const int    bsize  = b.size();

  vout.paranoiac(m_vl, "%s: starts\n", class_name.c_str());
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

      solve_step(b, rr);
      Nconv2 += Nconv_unit * m_N_M;

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
void Solver_GMRES_m_Cmplx::reset_field(const Field& b)
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

      m_v_tmp.reset(Nin, Nvol, Nex);

      m_v.resize(m_N_M + 1);

      for (int i = 0; i < m_N_M + 1; ++i) {
        m_v[i].reset(Nin, Nvol, Nex);
      }
    }
  }
#pragma omp barrier

  vout.paranoiac(m_vl, "  %s: field size reset.\n", class_name.c_str());
}


//====================================================================
void Solver_GMRES_m_Cmplx::solve_init(const Field& b, double& rr)
{
  copy(m_x, m_s);  // x  = s;

  for (int i = 0; i < m_N_M + 1; ++i) {
    m_v[i].set(0.0);           // m_v[i] = 0.0;
  }

  // r = b - A x_0
  m_fopr->mult(m_v_tmp, m_s);  // v_tmp = m_fopr->mult(s);
  copy(m_r, b);                // r  = b;
  axpy(m_r, -1.0, m_v_tmp);    // r -= v_tmp;

  rr = m_r.norm2();            // rr = r * r;

#pragma omp barrier
#pragma omp master
  m_beta_prev = sqrt(rr);
#pragma omp barrier

  //- v[0] = (1.0 / m_beta_prev) * r;
  copy(m_v[0], m_r);                  // v[0] = r;
  scal(m_v[0], (1.0 / m_beta_prev));  // v[0] = (1.0 / beta_p) * v[0];
}


//====================================================================
void Solver_GMRES_m_Cmplx::solve_step(const Field& b, double& rr)
{
  std::valarray<dcomplex> h((m_N_M + 1) * m_N_M), y(m_N_M);

  h = cmplx(0.0, 0.0);
  y = cmplx(0.0, 0.0);


  for (int j = 0; j < m_N_M; ++j) {
    m_fopr->mult(m_v_tmp, m_v[j]);  // v_tmp = m_fopr->mult(v[j]);

    for (int i = 0; i < j + 1; ++i) {
      int ij = index_ij(i, j);
      h[ij] = dotc(m_v[i], m_v_tmp);  // h[ij] = (v[i], A v[j]);
    }

    //- v[j+1] = A v[j] - \Sum_{i=0}^{j-1} h[i,j] * v[i]
    m_v[j + 1] = m_v_tmp;

    for (int i = 0; i < j + 1; ++i) {
      int ij = index_ij(i, j);
      axpy(m_v[j + 1], -h[ij], m_v[i]);  // v[j+1] -= h[ij] * v[i];
    }

    double v_norm2 = m_v[j + 1].norm2();

    int j1j = index_ij(j + 1, j);
    h[j1j] = cmplx(sqrt(v_norm2), 0.0);

    scal(m_v[j + 1], 1.0 / sqrt(v_norm2));  // v[j+1] /= sqrt(v_norm2);
  }


  // Compute y, which minimizes J := |r_new| = |beta_p - h * y|
  min_J(y, h);


  // x += Sum_{i=0}^{N_M-1} y[i] * v[i];
  for (int i = 0; i < m_N_M; ++i) {
    axpy(m_x, y[i], m_v[i]);  // x += y[i] * v[i];
  }


  // r = b - m_fopr->mult(x);
  copy(m_s, m_x);  // s = x;
  solve_init(b, rr);
}


//====================================================================
void Solver_GMRES_m_Cmplx::min_J(std::valarray<dcomplex>& y,
                                 std::valarray<dcomplex>& h)
{
  // Compute y, which minimizes J := |r_new| = |beta_p - h * y|

  std::valarray<dcomplex> g(m_N_M + 1);

  g    = cmplx(0.0, 0.0);
  g[0] = cmplx(m_beta_prev, 0.0);

  for (int i = 0; i < m_N_M; ++i) {
    int    ii    = index_ij(i, i);
    double h_1_r = abs(h[ii]);

    int    i1i   = index_ij(i + 1, i);
    double h_2_r = abs(h[i1i]);

    double denomi = sqrt(h_1_r * h_1_r + h_2_r * h_2_r);

    dcomplex cs = h[ii] / denomi;
    dcomplex sn = h[i1i] / denomi;

    for (int j = i; j < m_N_M; ++j) {
      int ij  = index_ij(i, j);
      int i1j = index_ij(i + 1, j);

      dcomplex const_1_c = conj(cs) * h[ij] + sn * h[i1j];
      dcomplex const_2_c = -sn * h[ij] + cs * h[i1j];

      h[ij]  = const_1_c;
      h[i1j] = const_2_c;
    }

    dcomplex const_1_c = conj(cs) * g[i] + sn * g[i + 1];
    dcomplex const_2_c = -sn * g[i] + cs * g[i + 1];

    g[i]     = const_1_c;
    g[i + 1] = const_2_c;
  }


  for (int i = m_N_M - 1; i > -1; --i) {
    for (int j = i + 1; j < m_N_M; ++j) {
      int ij = index_ij(i, j);
      g[i] -= h[ij] * y[j];
    }

    int ii = index_ij(i, i);
    y[i] = g[i] / h[ii];
  }
}


//====================================================================
double Solver_GMRES_m_Cmplx::flop_count()
{
  const int NPE = CommonParameters::NPE();

  //- NB1 Nin = 2 * Nc * Nd, Nex = 1  for field_F
  //- NB2 Nvol = CommonParameters::Nvol()/2 for eo
  const int Nin  = m_x.nin();
  const int Nvol = m_x.nvol();
  const int Nex  = m_x.nex();

  const double gflop_fopr = m_fopr->flop_count();

  if (gflop_fopr < CommonParameters::epsilon_criterion()) {
    vout.crucial(m_vl, "Warning at %s: no fopr->flop_count() is available, setting flop = 0\n", class_name.c_str());
    return 0.0;
  }

  const double gflop_axpy = (Nin * Nex * 2) * ((Nvol * NPE) / 1.0e+9);
  const double gflop_dotc = (Nin * Nex * 4) * ((Nvol * NPE) / 1.0e+9);
  const double gflop_norm = (Nin * Nex * 2) * ((Nvol * NPE) / 1.0e+9);
  const double gflop_scal = (Nin * Nex * 2) * ((Nvol * NPE) / 1.0e+9);

  int N_M_part = 0;
  for (int j = 0; j < m_N_M; ++j) {
    for (int i = 0; i < j + 1; ++i) {
      N_M_part += 1;
    }
  }

  const double gflop_init = gflop_fopr + gflop_axpy + gflop_norm;
  const double gflop_step = m_N_M * gflop_fopr + N_M_part * gflop_dotc
                            + (N_M_part + m_N_M) * gflop_axpy
                            + gflop_init;
  const double gflop_true_residual = gflop_fopr + gflop_axpy + gflop_norm;

  const int    N_iter = (m_Nconv_count - 1) / m_N_M;
  const double gflop  = gflop_norm + gflop_init
                        + gflop_step * N_iter + gflop_true_residual * (m_Nrestart_count + 1)
                        + gflop_init * m_Nrestart_count;

  return gflop;
}


//====================================================================
//============================================================END=====
