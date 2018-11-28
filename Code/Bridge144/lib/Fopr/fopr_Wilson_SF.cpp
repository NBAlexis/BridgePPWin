#include "BridgeLib_Private.h"

/*!
        @file    $Id:: fopr_Wilson_SF.cpp #$

        @brief

        @author  Yusuke Taniguchi
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 00:49:46 #$

        @version $LastChangedRevision: 1561 $
*/

#include "fopr_Wilson_SF.h"



#ifdef USE_FACTORY
namespace {
  Fopr *create_object()
  {
    return new Fopr_Wilson_SF();
  }


  bool init = Fopr::Factory_noarg::Register("Wilson_SF", create_object);
}
#endif



const std::string Fopr_Wilson_SF::class_name = "Fopr_Wilson_SF";

//====================================================================
void Fopr_Wilson_SF::set_parameters(const Parameters& params)
{
  const string str_vlevel = params.get_string("verbose_level");

  m_vl = vout.set_verbose_level(str_vlevel);

  //- fetch and check input parameters
  double           kappa;
  std::vector<int> bc;

  int err = 0;
  err += params.fetch_double("hopping_parameter", kappa);
  err += params.fetch_int_vector("boundary_condition", bc);

  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameter not found.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }

  set_parameters(kappa, bc);
}


//====================================================================
void Fopr_Wilson_SF::set_parameters(const double kappa, const std::vector<int> bc)
{
  //- print input parameters
  vout.general(m_vl, "%s:\n", class_name.c_str());
  vout.general(m_vl, "  kappa = %12.8f\n", kappa);
  for (int dir = 0; dir < m_Ndim; ++dir) {
    vout.general(m_vl, "  boundary[%d] = %2d\n", dir, bc[dir]);
  }

  //- range check
  // NB. kappa = 0 is allowed.
  assert(bc.size() == m_Ndim);

  //- store values
  m_kappa = kappa;

  m_boundary.resize(m_Ndim);
  for (int mu = 0; mu < m_Ndim; ++mu) {
    m_boundary[mu] = bc[mu];
  }

  //- propagate parameters
  m_fopr_w->set_parameters(m_kappa, m_boundary);
}


//====================================================================
void Fopr_Wilson_SF::DdagD(Field& w, const Field& f)
{
  Field w2(f.nin(), f.nvol(), f.nex());

  D(w2, f);
  mult_gm5(w, w2);
  D(w2, w);
  mult_gm5(w, w2);
}


//====================================================================
void Fopr_Wilson_SF::Ddag(Field& w, const Field& f)
{
  Field w2(f.nin(), f.nvol(), f.nex());

  mult_gm5(w, f);
  D(w2, w);
  mult_gm5(w, w2);
}


//====================================================================
void Fopr_Wilson_SF::H(Field& w, const Field& f)
{
  Field w2(f.nin(), f.nvol(), f.nex());

  D(w2, f);
  mult_gm5(w, w2);
}


//====================================================================
void Fopr_Wilson_SF::D(Field& w, const Field& f)
{
  Field w2(f);

  set_zero.set_boundary_zero(w2);
  m_fopr_w->D(w, w2);
  set_zero.set_boundary_zero(w);
}


//====================================================================
void Fopr_Wilson_SF::mult_gm5p(int mu, Field_F& v, const Field_F& w)
{
  Field_F w2(w);

  set_zero.set_boundary_zero(w2);
  m_fopr_w->mult_gm5p(mu, v, w2);
  set_zero.set_boundary_zero(v);
}


//====================================================================
double Fopr_Wilson_SF::flop_count()
{
  //- Counting of floating point operations.
  //  not implemented, yet.

  double flop = 0.0;

  return flop;
}


//====================================================================
//============================================================END=====
