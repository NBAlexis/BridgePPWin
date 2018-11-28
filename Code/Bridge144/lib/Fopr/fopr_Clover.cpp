#include "BridgeLib_Private.h"

/*!
        @file    $Id:: fopr_Clover.cpp #$

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 00:49:46 #$

        @version $LastChangedRevision: 1561 $
*/

#include "fopr_Clover.h"


using std::string;

#ifdef USE_FACTORY
namespace {
  Fopr *create_object()
  {
    return new Fopr_Clover();
  }


  Fopr *create_object_with_arg(const std::string& repr)
  {
    return new Fopr_Clover(repr);
  }


  bool init1 = Fopr::Factory_noarg::Register("Clover", create_object);
  bool init2 = Fopr::Factory_string::Register("Clover", create_object_with_arg);
}
#endif



const std::string Fopr_Clover::class_name = "Fopr_Clover";

//====================================================================
void Fopr_Clover::init(string repr)
{
  m_Nvol = CommonParameters::Nvol();
  m_Ndim = CommonParameters::Ndim();
  m_Nc   = CommonParameters::Nc();
  m_Nd   = CommonParameters::Nd();
  m_NinF = 2 * m_Nc * m_Nd;

  m_boundary.resize(m_Ndim);

  m_U = 0;

  m_repr = repr;

  m_fopr_w   = new Fopr_Wilson(repr);
  m_fopr_csw = new Fopr_CloverTerm(repr);

  m_v1.reset(m_NinF, m_Nvol, 1);
  m_v2.reset(m_NinF, m_Nvol, 1);
}


//====================================================================
void Fopr_Clover::tidyup()
{
  delete m_fopr_w;
  delete m_fopr_csw;
}


//====================================================================
void Fopr_Clover::set_parameters(const Parameters& params)
{
  const string str_vlevel = params.get_string("verbose_level");

  m_vl = vout.set_verbose_level(str_vlevel);

  //- fetch and check input parameters
  double           kappa, cSW;
  std::vector<int> bc;

  int err = 0;
  err += params.fetch_double("hopping_parameter", kappa);
  err += params.fetch_double("clover_coefficient", cSW);
  err += params.fetch_int_vector("boundary_condition", bc);

  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameter not found.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }

  set_parameters(kappa, cSW, bc);
}


//====================================================================
void Fopr_Clover::set_parameters(double kappa, double cSW, std::vector<int> bc)
{
  //- print input parameters
  vout.general(m_vl, "%s:\n", class_name.c_str());
  vout.general(m_vl, "  kappa = %12.8f\n", kappa);
  vout.general(m_vl, "  cSW   = %12.8f\n", cSW);
  for (int mu = 0; mu < m_Ndim; ++mu) {
    vout.general(m_vl, "  boundary[%d] = %2d\n", mu, bc[mu]);
  }

  //- range check
  // NB. kappa,cSW == 0 is allowed.
  assert(bc.size() == m_Ndim);

  //- store values
  m_kappa = kappa;
  m_cSW   = cSW;

  for (int mu = 0; mu < m_Ndim; ++mu) {
    m_boundary[mu] = bc[mu];
  }

  //- propagate parameters to components
  m_fopr_w->set_parameters(m_kappa, m_boundary);
  m_fopr_csw->set_parameters(m_kappa, m_cSW, m_boundary);
}


//====================================================================
double Fopr_Clover::flop_count()
{
  // Counting of floating point operations.
  // defined only for D, Dag, H, DDdag, DdagD which can be called
  // from the solver algorithms.
  // Since the flop_count() of Fopr_Wilson_eo defines flop of
  // (1 - Meo*Moe), flop of clover term is twice added together with
  // contribution of addition.

  int Lvol = CommonParameters::Lvol();

  double flop_w   = m_fopr_w->flop_count();
  double flop_csw = m_fopr_csw->flop_count();

  flop_csw += static_cast<double>(2 * m_Nc * m_Nd * Lvol);

  double flop = flop_w + flop_csw;

  if ((m_mode == "DdagD") || (m_mode == "DDdag")) flop += flop_csw;
  // for additional twice mult of clover term.

  return flop;
}


//====================================================================
void Fopr_Clover::D(Field& w, const Field& f)
{
  assert(f.nex() == 1);

  m_fopr_w->D(w, f);
  m_fopr_csw->mult_sigmaF(m_v1, f);
  axpy(w, -1.0, m_v1);  //  w -= m_v1;

#pragma omp barrier
}


//====================================================================
void Fopr_Clover::Ddag(Field& w, const Field& f)
{
  mult_gm5(w, f);
  D(m_v2, w);
  mult_gm5(w, m_v2);
}


//====================================================================
void Fopr_Clover::DdagD(Field& w, const Field& f)
{
  D(m_v2, f);
  mult_gm5(w, m_v2);
  D(m_v2, w);
  mult_gm5(w, m_v2);
}


//====================================================================
void Fopr_Clover::DDdag(Field& w, const Field& f)
{
  mult_gm5(m_v2, f);
  D(w, m_v2);
  mult_gm5(m_v2, w);
  D(w, m_v2);
}


//====================================================================
void Fopr_Clover::H(Field& w, const Field& f)
{
  D(m_v2, f);
  mult_gm5(w, m_v2);
}


//====================================================================
void Fopr_Clover::mult_isigma(Field_F& v, const Field_F& w,
                              const int mu, const int nu)
{
  m_fopr_csw->mult_isigma(v, w, mu, nu);
}


//====================================================================
//============================================================END=====
