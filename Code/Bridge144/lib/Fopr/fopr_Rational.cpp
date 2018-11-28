#include "BridgeLib_Private.h"

/*!
        @file    $Id:: fopr_Rational.cpp #$

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 00:49:46 #$

        @version $LastChangedRevision: 1561 $
*/

#include "fopr_Rational.h"


#ifdef USE_FACTORY
namespace {
  Fopr *create_object(Fopr *fopr)
  {
    return new Fopr_Rational(fopr);
  }


  bool init = Fopr::Factory_fopr::Register("Rational", create_object);
}
#endif



const std::string Fopr_Rational::class_name = "Fopr_Rational";

//====================================================================
void Fopr_Rational::set_parameters(const Parameters& params)
{
  const string str_vlevel = params.get_string("verbose_level");

  m_vl = vout.set_verbose_level(str_vlevel);

  //- fetch and check input parameters
  int    Np, n_exp, d_exp;
  double x_min, x_max;
  int    Niter;
  double Stop_cond;

  int err = 0;
  err += params.fetch_int("number_of_poles", Np);
  err += params.fetch_int("exponent_numerator", n_exp);
  err += params.fetch_int("exponent_denominator", d_exp);
  err += params.fetch_double("lower_bound", x_min);
  err += params.fetch_double("upper_bound", x_max);
  err += params.fetch_int("maximum_number_of_iteration", Niter);
  err += params.fetch_double("convergence_criterion_squared", Stop_cond);

  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameter not found.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }


  set_parameters(Np, n_exp, d_exp, x_min, x_max, Niter, Stop_cond);
}


//====================================================================
void Fopr_Rational::set_parameters(int Np, int n_exp, int d_exp,
                                   double x_min, double x_max,
                                   int Niter, double Stop_cond)
{
  //- print input parameters
  vout.general(m_vl, "%s:\n", class_name.c_str());
  vout.general(m_vl, "  Np        = %d\n", Np);
  vout.general(m_vl, "  n_exp     = %d\n", n_exp);
  vout.general(m_vl, "  d_exp     = %d\n", d_exp);
  vout.general(m_vl, "  x_min     = %12.8f\n", x_min);
  vout.general(m_vl, "  x_max     = %12.8f\n", x_max);
  vout.general(m_vl, "  Niter     = %d\n", Niter);
  vout.general(m_vl, "  Stop_cond = %8.2e\n", Stop_cond);

  //- range check
  int err = 0;
  err += ParameterCheck::non_zero(Np);
  err += ParameterCheck::non_zero(n_exp);
  err += ParameterCheck::non_zero(d_exp);
  // NB. x_min,x_max=0 is allowed.
  err += ParameterCheck::non_zero(Niter);
  err += ParameterCheck::square_non_zero(Stop_cond);

  if (err) {
    vout.crucial(m_vl, "Error at %s: parameter range check failed.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }

  //- store values
  m_Np        = Np;
  m_n_exp     = n_exp;
  m_d_exp     = d_exp;
  m_x_min     = x_min;
  m_x_max     = x_max;
  m_Niter     = Niter;
  m_Stop_cond = Stop_cond;

  //- post-process
  init_parameters();
}


//====================================================================
void Fopr_Rational::init_parameters()
{
  int Nin  = m_fopr->field_nin();
  int Nvol = m_fopr->field_nvol();
  int Nex  = m_fopr->field_nex();

  int Nshift = m_Np;

  double x_min2 = m_x_min * m_x_min;
  double x_max2 = m_x_max * m_x_max;

  m_cl.resize(m_Np);
  m_bl.resize(m_Np);

  m_xq.resize(m_Np);
  for (int i = 0; i < Nshift; ++i) {
    m_xq[i].reset(Nin, Nvol, Nex);
  }

  m_solver = new Shiftsolver_CG(m_fopr, m_Niter, m_Stop_cond);

  Math_Rational rational;
  rational.set_parameters(m_Np, m_n_exp, m_d_exp, x_min2, x_max2);
  rational.get_parameters(m_a0, m_bl, m_cl);

  vout.general(m_vl, " a0     = %18.14f\n", m_a0);
  for (int i = 0; i < m_Np; i++) {
    vout.general(m_vl, " bl[%d] = %18.14f  cl[%d] = %18.14f\n",
                 i, m_bl[i], i, m_cl[i]);
  }
}


//====================================================================
void Fopr_Rational::mult(Field& v, const Field& b)
{
  int    Nconv;
  double diff;

  assert(v.nin() == b.nin());
  assert(v.nvol() == b.nvol());
  assert(v.nex() == b.nex());

  vout.general(m_vl, "    Shift solver in rational function\n");
  vout.general(m_vl, "      Number of shift values = %d\n", m_Np);

  m_fopr->set_mode("DdagD");
  m_solver->solve(m_xq, m_cl, b, Nconv, diff);

  vout.general(m_vl, "      diff(max) = %22.15e  \n", diff);

  copy(v, b);
  scal(v, m_a0);               // v *= m_a0;
  for (int i = 0; i < m_Np; i++) {
    axpy(v, m_bl[i], m_xq[i]); // v += m_bl[i] * m_xq[i];
  }
}


//====================================================================
double Fopr_Rational::func(double x)
{
  double y = m_a0;

  for (int k = 0; k < m_Np; ++k) {
    y += m_bl[k] / (x + m_cl[k]);
  }

  return y;
}


//====================================================================
//============================================================END=====
