#include "BridgeLib_Private.h"

/*!
        @file    $Id:: force_F_Wilson_eo.cpp #$

        @brief

        @author  UEDA, Satoru  (sueda)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 00:49:46 #$

        @version $LastChangedRevision: 1561 $
*/

#include "force_F_Wilson_eo.h"


using Bridge::vout;



const std::string Force_F_Wilson_eo::class_name = "Force_F_Wilson_eo";

//====================================================================
void Force_F_Wilson_eo::set_parameters(const Parameters& params)
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
void Force_F_Wilson_eo::set_parameters(double kappa, const std::vector<int> bc)
{
  int Ndim = CommonParameters::Ndim();

  //- print input parameters
  vout.general(m_vl, "%s:\n", class_name.c_str());
  vout.general(m_vl, "  kappa = %12.8f\n", kappa);
  for (int mu = 0; mu < Ndim; ++mu) {
    vout.general(m_vl, "  boundary[%d] = %2d\n", mu, bc[mu]);
  }

  //- range check
  // NB. kappa == 0 is allowed.
  assert(bc.size() == Ndim);

  //- store values
  m_kappa = kappa;

  m_boundary.resize(Ndim);
  for (int mu = 0; mu < Ndim; ++mu) {
    m_boundary[mu] = bc[mu];
  }

  //- post-process
  m_fopr_w->set_parameters(m_kappa, m_boundary);
}


//====================================================================
void Force_F_Wilson_eo::force_udiv(Field& force_, const Field& eta_e)
{
  //int Nc    = CommonParameters::Nc();
  //int Nd    = CommonParameters::Nd();
  int Nvol  = CommonParameters::Nvol();
  int Ndim  = CommonParameters::Ndim();
  int Nvol2 = Nvol / 2;

  Field_G force1(Nvol, Ndim);

  Field_F eta_o(Nvol2, 1);
  Field_F eta(Nvol, 1);

  Field_F zeta_e(Nvol2, 1);
  Field_F zeta_o(Nvol2, 1);
  Field_F zeta(Nvol, 1);

  m_fopr_w->Meo(eta_o, eta_e, 1);
  m_index.mergeField(eta, eta_e, eta_o);

  m_fopr_w->set_mode("H");
  m_fopr_w->mult(zeta_e, eta_e);
  m_fopr_w->Meo(zeta_o, zeta_e, 1);
  m_index.mergeField(zeta, zeta_e, zeta_o);

  force_udiv1_impl(force1, zeta, eta);
  copy(force_, force1); // force_ = force1;

  force_udiv1_impl(force1, eta, zeta);
  axpy(force_, 1.0, force1); // force_ += force1;
}


//====================================================================
void Force_F_Wilson_eo::force_udiv1(Field& force_, const Field& zeta_, const Field& eta_)
{
  int Nvol = CommonParameters::Nvol();
  int Ndim = CommonParameters::Ndim();

  Field_G force(Nvol, Ndim);
  Field_F zeta(zeta_);
  Field_F eta(eta_);

  force_udiv1_impl(force, zeta, eta);

  copy(force_, force); // force_ = force;
}


//====================================================================
void Force_F_Wilson_eo::force_udiv1_impl(Field_G& force, const Field_F& zeta, const Field_F& eta)
{
  //int Nc   = CommonParameters::Nc();
  //int Nd   = CommonParameters::Nd();
  int Nvol = CommonParameters::Nvol();
  int Ndim = CommonParameters::Ndim();

  Field_F eta2(Nvol, 1);

  force.set(0.0);

  for (int mu = 0; mu < Ndim; ++mu) {
    m_fopr_w->gm5p(mu, eta2, eta);
    tensorProd_Field_F(force, mu, zeta, eta2);
  }

  scal(force, -m_kappa); // force *= -m_kappa;
}


//====================================================================
//============================================================END=====
