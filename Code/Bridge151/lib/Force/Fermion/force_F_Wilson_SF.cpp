/*!
        @file    force_F_Wilson_SF.cpp

        @brief

        @author  Yusuke Taniguchi
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/
#include "BridgeLib_Private.h"
#include "force_F_Wilson_SF.h"

const std::string Force_F_Wilson_SF::class_name = "Force_F_Wilson_SF";

//====================================================================
void Force_F_Wilson_SF::set_parameters(const Parameters& params)
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
void Force_F_Wilson_SF::set_parameters(const double kappa, const std::vector<int> bc)
{
  const int Ndim = CommonParameters::Ndim();

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
  m_boundary = bc;

  //- propagate parameters
  m_fopr_w->set_parameters(m_kappa, m_boundary);
}


//====================================================================
void Force_F_Wilson_SF::force_udiv(Field& force_, const Field& eta_)
{
  const int Nvol = CommonParameters::Nvol();
  const int Ndim = CommonParameters::Ndim();

  Field_F zeta;
  Field_F eta(eta_);

  m_fopr_w->set_mode("H");
  m_fopr_w->mult(zeta, eta);

  Field_G force(Nvol, Ndim);
  force_udiv1_impl(force, zeta, eta);

  Field_G force1(Nvol, Ndim);
  force_udiv1_impl(force1, eta, zeta);
  axpy(force, 1.0, force1); // force += force1;

  Field_G_SF set_wk;
  set_wk.set_boundary_spatial_link_zero(force);

  copy(force_, force); // force_ = force;
}


//====================================================================
void Force_F_Wilson_SF::force_udiv1(Field& force_, const Field& zeta_, const Field& eta_)
{
  const int Nvol = CommonParameters::Nvol();
  const int Ndim = CommonParameters::Ndim();

  Field_G force(Nvol, Ndim);
  Field_F zeta(zeta_);
  Field_F eta(eta_);

  force_udiv1_impl(force, zeta, eta);

  copy(force_, force); // force_ = force;
}


//====================================================================
void Force_F_Wilson_SF::force_udiv1_impl(Field_G& force, const Field_F& zeta, const Field_F& eta)
{
  const int Ndim = CommonParameters::Ndim();

  force.set(0.0);

  for (int mu = 0; mu < Ndim; ++mu) {
    Field_F eta2;
    m_fopr_w->mult_gm5p(mu, eta2, eta);

    Field_F eta3;
    mult_Field_Gd(eta3, 0, *m_U, mu, eta2, 0);

    tensorProd_Field_F(force, mu, zeta, eta3);
  }

  scal(force, -m_kappa); // force *= -m_kappa;
}


//====================================================================
//============================================================END=====
