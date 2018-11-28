#include "BridgeLib_Private.h"

/*!
        @file    $Id:: gradientFlow.cpp #$

        @brief

        @author  Sinya Aoki
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 18:35:38 #$

        @version $LastChangedRevision: 1571 $
*/

#include "gradientFlow.h"

const std::string GradientFlow::class_name = "GradientFlow";

//====================================================================
GradientFlow::GradientFlow(Action *action)
  : m_vl(CommonParameters::Vlevel()),
    m_action(action),
    m_impl(0)
{
  initialize();
}


//====================================================================
GradientFlow::GradientFlow(unique_ptr<Action>& action)
  : m_vl(CommonParameters::Vlevel()),
    m_action(action.get()),
    m_impl(0)
{
  initialize();
}


//====================================================================
void GradientFlow::initialize()
{
  m_Nprec = 0;
  m_Estep = 0.0;
}


//====================================================================
GradientFlow::~GradientFlow()
{
  if (m_impl) delete m_impl;
}


//====================================================================
void GradientFlow::set_parameters(const Parameters& params)
{
  const string str_vlevel = params.get_string("verbose_level");

  m_vl = vout.set_verbose_level(str_vlevel);

  //- fetch and check input parameters
  int    Norder_RK;
  double Estep;
  int    Nprec;
  int    adaptive;
  double tolerance, safety;

  int err = 0;
  err += params.fetch_int("order_of_RungeKutta", Norder_RK);
  err += params.fetch_double("step_size", Estep);
  err += params.fetch_int("order_of_approx_for_exp_iP", Nprec);

  err += params.fetch_int("adaptive", adaptive);

  if (adaptive != 0) {  // mandatory only if adaptive is on.
    err += params.fetch_double("tolerance", tolerance);
    err += params.fetch_double("safety_factor", safety);
  } else {
    params.fetch_double("tolerance", tolerance);
    params.fetch_double("safety_factor", safety);
  }

  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameter not found\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }

  set_parameters(Norder_RK, Estep, Nprec, adaptive, tolerance, safety);
}


//====================================================================
void GradientFlow::set_parameters(
  const int Norder_RK,
  const double Estep, const int Nprec,
  const int adaptive,
  const double tolerance, const double safety)
{
  //- print input parameters
  vout.general(m_vl, "%s:\n", class_name.c_str());
  vout.general(m_vl, "  Norder_RK  = %d\n", Norder_RK);
  vout.general(m_vl, "  Estep      = %10.6f\n", Estep);
  vout.general(m_vl, "  Nprec      = %d\n", Nprec);
  vout.general(m_vl, "  adaptive   = %d\n", adaptive);
  vout.general(m_vl, "  tolerance  = %10.6e\n", tolerance);
  vout.general(m_vl, "  safety     = %10.6f\n", safety);


  //- range check
  int err = 0;
  err += ParameterCheck::non_negative(Norder_RK);
  err += ParameterCheck::square_non_zero(Estep);
  err += ParameterCheck::non_negative(Nprec);
  err += ParameterCheck::non_negative(tolerance);
  err += ParameterCheck::non_negative(safety);

  if (err) {
    vout.crucial(m_vl, "Error at %s: parameter range check failed\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }

  //- store values
  m_Norder_RK = Norder_RK; // order of Runge-Kutta
  m_Estep     = Estep;     // step size (SA)
  m_Nprec     = Nprec;     // order of approximation for e^{iP} (SA)

  m_is_adaptive = (adaptive == 0) ? false : true;
  m_tolerance   = tolerance;
  m_safety      = safety;

  set_parameter_Norder_RK(Norder_RK, m_is_adaptive);
}


//====================================================================
void GradientFlow::set_parameter_Norder_RK(const int Norder_RK, bool is_adaptive)
{
  m_Norder_RK = Norder_RK;

  switch (Norder_RK)
  {
  case 1:
    m_impl = new GradientFlow_RungeKutta_1st(m_action, m_Nprec, m_vl);
    break;

  case 2:
    m_impl = new GradientFlow_RungeKutta_2nd(m_action, m_Nprec, m_vl);
    break;

  case 3:
    m_impl = new GradientFlow_RungeKutta_3rd(m_action, m_Nprec, m_vl);
    break;

  case 4:
    m_impl = new GradientFlow_RungeKutta_4th(m_action, m_Nprec, m_vl);
    break;

  default:
    vout.crucial(m_vl, "Error at %s: order of Runge-Kutta is out of range\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }

  if (is_adaptive) {
    m_impl = new GradientFlow_AdaptiveRungeKutta(m_impl, m_tolerance, m_safety, m_vl);
  }
}


//====================================================================
double GradientFlow::evolve(double& t, Field_G& U)
{
  double plaq = m_staple.plaquette(U);

  vout.general(m_vl, "\n");
  vout.general(m_vl, "  plaq_org = %.16f\n", plaq); // write-out

  //- time evolution (SA)
  m_impl->flow(t, m_Estep, U);

  plaq = m_staple.plaquette(U);

  vout.general(m_vl, "  (t, plaq) = %.8f %.16f\n", t, plaq);

  return t * t * 36.0 * (1.0 - plaq);
}


//====================================================================
//============================================================END=====
