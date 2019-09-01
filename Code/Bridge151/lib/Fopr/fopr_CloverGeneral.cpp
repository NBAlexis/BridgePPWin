/*!
        @file    fopr_CloverGeneral.cpp

        @brief

        @author  Yusuke Namekawa (namekawa)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/
#include "BridgeLib_Private.h"
#include "fopr_CloverGeneral.h"

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init = Fopr_CloverGeneral::register_factory();
}
#endif

const std::string Fopr_CloverGeneral::class_name = "Fopr_CloverGeneral";

//====================================================================
void Fopr_CloverGeneral::init(const std::string repr)
{
  m_Nvol = CommonParameters::Nvol();
  m_Ndim = CommonParameters::Ndim();
  m_Nc   = CommonParameters::Nc();
  m_Nd   = CommonParameters::Nd();
  m_NinF = 2 * m_Nc * m_Nd;

  m_boundary.resize(m_Ndim);

  m_U = 0;

  m_repr = repr;

  m_fopr_w   = new Fopr_WilsonGeneral(repr);
  m_fopr_csw = new Fopr_CloverTerm_General(repr);

  m_v1.reset(m_NinF, m_Nvol, 1);
  m_v2.reset(m_NinF, m_Nvol, 1);
}


//====================================================================
void Fopr_CloverGeneral::tidyup()
{
  delete m_fopr_w;
  delete m_fopr_csw;
}


//====================================================================
void Fopr_CloverGeneral::set_parameters(const Parameters& params)
{
  const std::string str_vlevel = params.get_string("verbose_level");

  m_vl = vout.set_verbose_level(str_vlevel);

  //- fetch and check input parameters
  double           kappa_s, kappa_t;
  double           nu_s, r_s;
  double           cSW_s, cSW_t;
  std::vector<int> bc;

  int err = 0;
  err += params.fetch_double("hopping_parameter_spatial", kappa_s);
  err += params.fetch_double("hopping_parameter_temporal", kappa_t);
  err += params.fetch_double("dispersion_parameter_spatial", nu_s);
  err += params.fetch_double("Wilson_parameter_spatial", r_s);
  err += params.fetch_double("clover_coefficient_spatial", cSW_s);
  err += params.fetch_double("clover_coefficient_temporal", cSW_t);

  err += params.fetch_int_vector("boundary_condition", bc);

  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameter not found.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }

  set_parameters(kappa_s, kappa_t, nu_s, r_s, cSW_s, cSW_t, bc);
}


//====================================================================
void Fopr_CloverGeneral::set_parameters(const double kappa_s, const double kappa_t,
                                        const double nu_s, const double r_s,
                                        const double cSW_s, const double cSW_t,
                                        const std::vector<int> bc)
{
  //- print input parameters
  vout.general(m_vl, "%s:\n", class_name.c_str());
  vout.general(m_vl, "  kappa_s = %12.8f\n", kappa_s);
  vout.general(m_vl, "  kappa_t = %12.8f\n", kappa_t);
  vout.general(m_vl, "  nu_s    = %12.8f\n", nu_s);
  vout.general(m_vl, "  r_s     = %12.8f\n", r_s);
  vout.general(m_vl, "  cSW_s   = %12.8f\n", cSW_s);
  vout.general(m_vl, "  cSW_t   = %12.8f\n", cSW_t);
  for (int mu = 0; mu < m_Ndim; ++mu) {
    vout.general(m_vl, "  boundary[%d] = %2d\n", mu, bc[mu]);
  }

  //- range check
  // NB. kappa,cSW == 0 is allowed.
  assert(bc.size() == m_Ndim);

  //- store values
  m_kappa_s = kappa_s;
  m_kappa_t = kappa_t;
  m_nu_s    = nu_s;
  m_r_s     = r_s;
  m_cSW_s   = cSW_s;
  m_cSW_t   = cSW_t;

  // m_boundary.resize(m_Ndim);  // NB. already resized in init.
  m_boundary = bc;

  //- propagate parameters to components
  m_fopr_w->set_parameters(m_kappa_s, m_kappa_t, m_nu_s, m_r_s, m_boundary);
  m_fopr_csw->set_parameters(m_kappa_s, m_kappa_t, m_cSW_s, m_cSW_t, m_boundary);
}


//====================================================================
void Fopr_CloverGeneral::D(Field& w, const Field& f)
{
  assert(f.nex() == 1);

  m_fopr_w->D(w, f);
  m_fopr_csw->mult_sigmaF(m_v1, f);
  axpy(w, -1.0, m_v1);  //  w -= m_v1;

#pragma omp barrier
}


//====================================================================
void Fopr_CloverGeneral::Ddag(Field& w, const Field& f)
{
  mult_gm5(w, f);
  D(m_v2, w);
  mult_gm5(w, m_v2);
}


//====================================================================
void Fopr_CloverGeneral::DdagD(Field& w, const Field& f)
{
  D(m_v2, f);
  mult_gm5(w, m_v2);
  D(m_v2, w);
  mult_gm5(w, m_v2);
}


//====================================================================
void Fopr_CloverGeneral::DDdag(Field& w, const Field& f)
{
  mult_gm5(m_v2, f);
  D(w, m_v2);
  mult_gm5(m_v2, w);
  D(w, m_v2);
}


//====================================================================
void Fopr_CloverGeneral::H(Field& w, const Field& f)
{
  D(m_v2, f);
  mult_gm5(w, m_v2);
}


//====================================================================
void Fopr_CloverGeneral::mult_isigma(Field_F& v, const Field_F& w,
                                     const int mu, const int nu)
{
  m_fopr_csw->mult_isigma(v, w, mu, nu);
}


//====================================================================
double Fopr_CloverGeneral::flop_count()
{
  // Counting of floating point operations in giga unit.
  // defined only for D, Dag, H, DDdag, DdagD which can be called
  // from the solver algorithms.
  // Since the flop_count() of Fopr_Wilson_eo defines flop of
  // (1 - Meo*Moe), flop of clover term is twice added together with
  // contribution of addition.

  const int Nvol = CommonParameters::Nvol();
  const int NPE  = CommonParameters::NPE();

  const double gflop_w = 2 * m_fopr_w->flop_count();

  double gflop_csw = m_fopr_csw->flop_count();

  gflop_csw += 2 * m_Nc * m_Nd / Nvol / NPE / 1.0e+9;

  double gflop = gflop_w + gflop_csw;

  // Additional twice mult of clover term
  if ((m_mode == "DdagD") || (m_mode == "DDdag")) gflop += gflop_csw;

  return gflop;
}


//====================================================================
//============================================================END=====
