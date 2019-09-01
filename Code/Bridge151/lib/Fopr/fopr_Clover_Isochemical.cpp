/*!
        @file    fopr_Clover_Isochemical.cpp

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/
#include "BridgeLib_Private.h"
#include "fopr_Clover_Isochemical.h"

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init = Fopr_Clover_Isochemical::register_factory();
}
#endif

const std::string Fopr_Clover_Isochemical::class_name = "Fopr_Clover_Isochemical";

//====================================================================
void Fopr_Clover_Isochemical::set_parameters(const Parameters& params)
{
  const std::string str_vlevel = params.get_string("verbose_level");

  m_vl = vout.set_verbose_level(str_vlevel);

  //- fetch and check input parameters
  double           kappa, cSW, mu;
  std::vector<int> bc;

  int err = 0;
  err += params.fetch_double("hopping_parameter", kappa);
  err += params.fetch_double("clover_coefficient", cSW);
  err += params.fetch_double("isospin_chemical_potential", mu);
  err += params.fetch_int_vector("boundary_condition", bc);

  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameter not found.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }


  set_parameters(kappa, cSW, mu, bc);
}


//====================================================================
void Fopr_Clover_Isochemical::set_parameters(const double kappa, const double cSW, const double mu,
                                             const std::vector<int> bc)
{
  //- print input parameters
  vout.general(m_vl, "%s:\n", class_name.c_str());
  vout.general(m_vl, "  kappa = %12.8f\n", kappa);
  vout.general(m_vl, "  cSW   = %12.8f\n", cSW);
  vout.general(m_vl, "  mu    = %12.8f\n", mu);
  for (int dir = 0; dir < m_Ndim; ++dir) {
    vout.general(m_vl, "  boundary[%d] = %2d\n", dir, bc[dir]);
  }

  //- range check
  // NB. kappa,cSW,mu == 0 is allowed.
  assert(bc.size() == m_Ndim);

  //- store values
  m_kappa = kappa;
  m_cSW   = cSW;
  m_mu    = mu;

  // m_boundary.resize(m_Ndim);  // NB. already resized in init.
  m_boundary = bc;

  //- propagate parameters
  m_fopr_w->set_parameters(m_kappa, m_mu, m_boundary);
  m_fopr_csw->set_parameters(m_kappa, m_cSW, m_boundary);
}


//====================================================================
void Fopr_Clover_Isochemical::init(const std::string repr)
{
  m_Nvol = CommonParameters::Nvol();
  m_Ndim = CommonParameters::Ndim();
  m_Nc   = CommonParameters::Nc();
  m_Nd   = CommonParameters::Nd();
  m_NinF = 2 * m_Nc * m_Nd;

  m_boundary.resize(m_Ndim);

  m_U = 0;

  m_repr = repr;

  m_fopr_w   = new Fopr_Wilson_Isochemical(repr);
  m_fopr_csw = new Fopr_CloverTerm(repr);

  m_w1.reset(m_NinF, m_Nvol, 1);
  m_w2.reset(m_NinF, m_Nvol, 1);
}


//====================================================================
void Fopr_Clover_Isochemical::tidyup()
{
  delete m_fopr_w;
  delete m_fopr_csw;
}


//====================================================================
void Fopr_Clover_Isochemical::set_config(Field *U)
{
  m_U = (Field_G *)U;

  m_fopr_w->set_config(U);
  m_fopr_csw->set_config(U);
}


//====================================================================
void Fopr_Clover_Isochemical::D(Field& w, const Field& f)
{
  assert(f.nex() == 1);

  m_fopr_w->D(w, f);
  m_fopr_csw->mult_sigmaF(m_w1, f);
  axpy(w, -1.0, m_w1);  //  w -= m_w1;

#pragma omp barrier
}


//====================================================================
void Fopr_Clover_Isochemical::Ddag(Field& w, const Field& f)
{
  assert(f.nex() == 1);

  m_fopr_w->Ddag(w, f);
  m_fopr_csw->mult_sigmaF(m_w1, f);
  axpy(w, -1.0, m_w1);  //  w -= m_w1;

#pragma omp barrier
}


//====================================================================
void Fopr_Clover_Isochemical::DdagD(Field& w, const Field& f)
{
  assert(f.nex() == 1);

  D(m_w2, f);
  Ddag(w, m_w2);
}


//====================================================================
void Fopr_Clover_Isochemical::H(Field& w, const Field& f)
{
  assert(f.nex() == 1);

  D(m_w2, f);
  mult_gm5(w, m_w2);
}


//====================================================================
void Fopr_Clover_Isochemical::Hdag(Field& w, const Field& f)
{
  assert(f.nex() == 1);

  mult_gm5(m_w2, f);
  Ddag(w, m_w2);
}


//====================================================================
void Fopr_Clover_Isochemical::mult_isigma(Field_F& v, const Field_F& w,
                                          const int mu, const int nu)
{
  m_fopr_csw->mult_isigma(v, w, mu, nu);
}


//====================================================================
double Fopr_Clover_Isochemical::flop_count()
{
  //- Counting of floating point operations in giga unit.
  //  not implemented, yet.

  const double gflop = 0.0;

  return gflop;
}


//====================================================================
//============================================================END=====
