/*!
        @file    fopr_Wilson_SF.cpp

        @brief

        @author  Yusuke Taniguchi
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/
#include "BridgeLib_Private.h"
#include "fopr_Wilson_SF.h"

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init = Fopr_Wilson_SF::register_factory();
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
  m_boundary = bc;

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
void Fopr_Wilson_SF::mult_gm5p(const int mu, Field_F& v, const Field_F& w)
{
  Field_F w2(w);

  set_zero.set_boundary_zero(w2);
  m_fopr_w->mult_gm5p(mu, v, w2);
  set_zero.set_boundary_zero(v);
}


//====================================================================
double Fopr_Wilson_SF::flop_count()
{
  //- Counting of floating point operations in giga unit.
  //  not implemented, yet.

  const double gflop = 0;

  return gflop;
}


//====================================================================
//============================================================END=====
