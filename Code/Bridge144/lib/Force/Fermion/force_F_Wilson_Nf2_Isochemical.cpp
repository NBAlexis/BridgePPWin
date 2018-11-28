#include "BridgeLib_Private.h"

/*!
        @file    $Id:: force_F_Wilson_Nf2_Isochemical.cpp #$

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 00:49:46 #$

        @version $LastChangedRevision: 1561 $
*/

#include "force_F_Wilson_Nf2_Isochemical.h"



const std::string Force_F_Wilson_Nf2_Isochemical::class_name = "Force_F_Wilson_Nf2_Isochemical";

//====================================================================
void Force_F_Wilson_Nf2_Isochemical::set_parameters(const Parameters& params)
{
  const string str_vlevel = params.get_string("verbose_level");

  m_vl = vout.set_verbose_level(str_vlevel);

  //- fetch and check input parameters
  double           kappa, mu;
  std::vector<int> bc;

  int err = 0;
  err += params.fetch_double("hopping_parameter", kappa);
  err += params.fetch_double("isospin_chemical_potential", mu);
  err += params.fetch_int_vector("boundary_condition", bc);

  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameter not found.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }


  set_parameters(kappa, mu, bc);
}


//====================================================================
void Force_F_Wilson_Nf2_Isochemical::set_parameters(const double kappa, const double mu,
                                                    const std::vector<int> bc)
{
  int Ndim = CommonParameters::Ndim();

  //- print input parameters
  vout.general(m_vl, "%s:\n", class_name.c_str());
  vout.general(m_vl, "  kappa = %12.8f\n", kappa);
  vout.general(m_vl, "  mu    = %12.8f\n", mu);
  for (int imu = 0; imu < Ndim; ++imu) {
    vout.general(m_vl, "  boundary[%d] = %2d\n", imu, bc[imu]);
  }

  //- range check
  // NB. kappa,mu == 0 is allowed.
  assert(bc.size() == Ndim);

  //- store values
  m_kappa  = kappa;
  m_mu     = mu;
  m_exp_mu = exp(mu);

  m_boundary.resize(Ndim);
  for (int dir = 0; dir < Ndim; ++dir) {
    m_boundary[dir] = bc[dir];
  }

  //- post-process
  m_fopr_w->set_parameters(m_kappa, m_mu, m_boundary);
}


//====================================================================
void Force_F_Wilson_Nf2_Isochemical::force_udiv(Field& force_, const Field& eta_)
{
  //int Nc   = CommonParameters::Nc();
  //int Nd   = CommonParameters::Nd();
  int Nvol = CommonParameters::Nvol();
  int Ndim = CommonParameters::Ndim();

  Field_G force(Nvol, Ndim);
  Field_F zeta(Nvol, 1);
  Field_F eta(eta_);

  m_fopr_w->set_mode("H");
  m_fopr_w->mult(zeta, eta);

  set_mode("H");
  force_udiv1_impl(force, zeta, eta);
  copy(force_, force); // force_ = force;

  set_mode("Hdag");
  force_udiv1_impl(force, eta, zeta);
  axpy(force_, 1.0, force); // force_ += force;
}


//====================================================================
void Force_F_Wilson_Nf2_Isochemical::force_udiv1(Field& force_, const Field& zeta_, const Field& eta_)
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
void Force_F_Wilson_Nf2_Isochemical::force_udiv1_impl(Field_G& force, const Field_F& zeta, const Field_F& eta)
{
  //int Nc   = CommonParameters::Nc();
  //int Nd   = CommonParameters::Nd();
  int Nvol = CommonParameters::Nvol();
  int Ndim = CommonParameters::Ndim();

  Field_F eta2(Nvol, 1), eta3(Nvol, 1);

  force.set(0.0);

  for (int dir = 0; dir < Ndim - 1; ++dir) {
    m_fopr_w->mult_gm5p(dir, eta2, eta);
    mult_Field_Gd(eta3, 0, *m_U, dir, eta2, 0);
    scal(eta3, -m_kappa); // eta3 *= -m_kappa;
    tensorProd_Field_F(force, dir, zeta, eta3);
  }

  {
    int dir = Ndim - 1;
    m_fopr_w->mult_gm5p(dir, eta2, eta);
    mult_Field_Gd(eta3, 0, *m_U, dir, eta2, 0);
    if (m_mode == "H") {
      scal(eta3, -(m_kappa * m_exp_mu)); // eta3 *= -(m_kappa * m_exp_mu);
    } else if (m_mode == "Hdag") {
      scal(eta3, -(m_kappa / m_exp_mu)); // eta3 *= -(m_kappa / m_exp_mu);
    } else {
      vout.crucial(m_vl, "Error at %s: illegal mode.\n", class_name.c_str());
      exit(EXIT_FAILURE);
    }

    tensorProd_Field_F(force, dir, zeta, eta3);
  }
}


//====================================================================
//============================================================END=====
