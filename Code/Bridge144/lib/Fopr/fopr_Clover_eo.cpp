#include "BridgeLib_Private.h"

/*!
        @file    $Id:: fopr_Clover_eo.cpp #$

        @brief

        @author  Satoru Ueda  (sueda)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 00:49:46 #$

        @version $LastChangedRevision: 1561 $
*/

#include "fopr_Clover_eo.h"

#include "ResourceManager/threadManager_OpenMP.h"


#ifdef USE_FACTORY
namespace {
  Fopr *create_object_with_repr(const std::string& repr)
  {
    return new Fopr_Clover_eo(repr);
  }


  bool init = Fopr::Factory_string::Register("Clover_eo", create_object_with_repr);
}
#endif

//====================================================================


const std::string Fopr_Clover_eo::class_name = "Fopr_Clover_eo";

//====================================================================
void Fopr_Clover_eo::init(const std::string repr)
{
  m_Nc    = CommonParameters::Nc();
  m_Nd    = CommonParameters::Nd();
  m_Ndim  = CommonParameters::Ndim();
  m_NinF  = 2 * m_Nc * m_Nd;
  m_Nvol  = CommonParameters::Nvol();
  m_Nvol2 = m_Nvol / 2;

  m_boundary.resize(m_Ndim);

  m_Ueo = new Field_G(m_Nvol, m_Ndim);

  m_fopr_w   = new Fopr_Wilson_eo(repr);
  m_fopr_csw = new Fopr_CloverTerm_eo(repr);

  // working field (Field)
  m_w1.reset(m_NinF, m_Nvol2, 1);

  // working field (Field_F)
  m_vF1.reset(m_Nvol2, 1);
  m_vF2.reset(m_Nvol2, 1);
  m_vF3.reset(m_Nvol2, 1);
}


//====================================================================
void Fopr_Clover_eo::tidyup()
{
  delete m_fopr_w;
  delete m_fopr_csw;

  delete m_Ueo;
}


//====================================================================
void Fopr_Clover_eo::set_parameters(const Parameters& params)
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
void Fopr_Clover_eo::set_parameters(const double kappa, const double cSW,
                                    const std::vector<int> bc)
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

  assert(bc.size() == m_Ndim);
  for (int mu = 0; mu < m_Ndim; ++mu) {
    m_boundary[mu] = bc[mu];
  }

  //- propagate parameters to components
  m_fopr_w->set_parameters(m_kappa, m_boundary);
  m_fopr_csw->set_parameters(m_kappa, m_cSW, m_boundary);
}


//====================================================================
void Fopr_Clover_eo::set_config(Field *U)
{
  int Ndim = CommonParameters::Ndim();
  int Nvol = CommonParameters::Nvol();

  m_idx.convertField(*m_Ueo, *U);

  m_fopr_w->set_config(U);
  m_fopr_csw->set_config(m_Ueo);
}


//====================================================================
void Fopr_Clover_eo::D(Field& v, const Field& f)
{
  Meo(m_vF2, f, 1);
  Meo(v, m_vF2, 0);
  aypx(-1.0, v, f);

#pragma omp barrier
}


//====================================================================
void Fopr_Clover_eo::Ddag(Field& v, const Field& f)
{
  Mdageo(m_vF2, f, 1);
  Mdageo(m_vF3, m_vF2, 0);
  copy(v, f);
  axpy(v, -1.0, m_vF3);

#pragma omp barrier
}


//====================================================================
void Fopr_Clover_eo::DdagD(Field& v, const Field& f)
{
  D(m_w1, f);
  Ddag(v, m_w1);
}


//====================================================================
void Fopr_Clover_eo::DDdag(Field& v, const Field& f)
{
  Ddag(m_w1, f);
  D(v, m_w1);
}


//====================================================================
void Fopr_Clover_eo::H(Field& v, const Field& f)
{
  D(m_w1, f);
  mult_gm5(v, m_w1);
}


//====================================================================
void Fopr_Clover_eo::mult_gm5(Field& v, const Field& f)
{
  m_fopr_w->mult_gm5(v, f);
}


//====================================================================
void Fopr_Clover_eo::MeoMoe(Field& v, const Field& f)
{
  Meo(m_vF2, f, 1);
  Meo(m_vF3, m_vF2, 0);
  axpy(v, -1.0, m_vF3);

#pragma omp barrier
}


//====================================================================
double Fopr_Clover_eo::flop_count()
{
  // Counting of floating point operations.
  // defined only for D, Dag, H, DDdag, DdagD which can be called
  // from the solver algorithms.
  // Since the flop_count() of Fopr_Wilson_eo defines flop of
  // (1 - Meo*Moe), flop of clover term is twice added together with
  // contribution of addition.

  int Lvol = CommonParameters::Lvol();

  double flop_w = m_fopr_w->flop_count();
  // this is for aypx + Meo * 2 with Wilson_eo.

  double flop_csw = m_fopr_csw->flop_count();
  // this is for inversion of (1 - clover term).

  double flop = flop_w + 2.0 * flop_csw;

  if ((m_mode == "DdagD") || (m_mode == "DDdag")) flop += 2.0 * flop_csw;
  // for additional twice mult of clover term.

  return flop;
}


//====================================================================
void Fopr_Clover_eo::Meo(Field& v, const Field& f, const int ieo)
{
  // ieo=0: even <-- odd
  // ieo=1: odd  <-- even

  m_fopr_w->Meo(m_vF1, f, ieo);
  m_fopr_csw->mult_csw_inv(v, m_vF1, ieo);
}


//====================================================================
void Fopr_Clover_eo::Mdageo(Field_F& v, const Field_F& f, const int ieo)
{
  m_fopr_w->mult_gm5(m_vF1, f);
  m_fopr_w->Meo(v, m_vF1, ieo);
  m_fopr_w->mult_gm5(m_vF1, v);

  m_fopr_csw->mult_csw_inv(v, m_vF1, ieo);
}


//====================================================================
void Fopr_Clover_eo::Meo_gm5(Field_F& v,
                             const Field_F& f, const int ieo)
{
  m_fopr_w->Meo(v, f, ieo);
  m_fopr_csw->mult_csw_inv(m_vF1, v, ieo);

  mult_gm5(v, m_vF1);
}


//====================================================================
void Fopr_Clover_eo::mult_isigma(Field_F& w, const Field_F& f,
                                 const int mu, const int nu)
{
  m_fopr_csw->mult_isigma(w, f, mu, nu);
}


//====================================================================
void Fopr_Clover_eo::prePropD(Field& Be, Field& bo,
                              const Field& b)
{
  ThreadManager_OpenMP::assert_single_thread(class_name);

  Field_F be(m_Nvol2, 1);
  Field_F tmp(m_Nvol2, 1);

  m_idx.convertField(tmp, b, 1);
  m_fopr_csw->mult_csw_inv(be, tmp, 1);
  copy(bo, be);

  m_idx.convertField(tmp, b, 0);
  m_fopr_csw->mult_csw_inv(be, tmp, 0);

  copy(Be, be);
  Meo(tmp, bo, 0);
  axpy(Be, -1.0, tmp);

#pragma omp barrier
}


//====================================================================
void Fopr_Clover_eo::postPropD(Field& x,
                               const Field& xe, const Field& bo)
{
  ThreadManager_OpenMP::assert_single_thread(class_name);

  assert(x.nin() == m_NinF);
  assert(x.nvol() == m_Nvol);
  assert(x.nex() == 1);

  Field_F xo(m_Nvol2, 1);

  Meo(xo, xe, 1);
  aypx(-1.0, xo, bo);

  m_idx.reverseField(x, xe, 0);
  m_idx.reverseField(x, xo, 1);

#pragma omp barrier
}


//====================================================================
void Fopr_Clover_eo::prePropDag(Field& Be, Field& bo,
                                const Field& b)
{
  ThreadManager_OpenMP::assert_single_thread(class_name);

  Field_F be(m_Nvol2, 1);
  Field_F tmp(m_Nvol2, 1);

  m_idx.convertField(tmp, b, 1);
  m_fopr_csw->mult_csw_inv(be, tmp, 1);
  copy(bo, be);

  m_idx.convertField(tmp, b, 0);
  m_fopr_csw->mult_csw_inv(be, tmp, 0);

  copy(Be, be);
  Mdageo(tmp, bo, 0);
  axpy(Be, -1.0, tmp);

#pragma omp barrier
}


//====================================================================
void Fopr_Clover_eo::postPropDag(Field& x,
                                 const Field& xe, const Field& bo)
{
  ThreadManager_OpenMP::assert_single_thread(class_name);

  assert(x.nin() == m_NinF);
  assert(x.nvol() == m_Nvol);
  assert(x.nex() == 1);

  Field_F xo(m_Nvol2, 1);

  Mdageo(xo, xe, 1);
  aypx(-1.0, xo, bo);

  m_idx.reverseField(x, xe, 0);
  m_idx.reverseField(x, xo, 1);

#pragma omp barrier
}


//====================================================================
//============================================================END=====
