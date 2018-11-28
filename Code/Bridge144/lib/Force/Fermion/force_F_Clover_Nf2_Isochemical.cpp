#include "BridgeLib_Private.h"

/*!
        @file    $Id:: force_F_Clover_Nf2_Isochemical.cpp #$

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 00:49:46 #$

        @version $LastChangedRevision: 1561 $
*/

#include "force_F_Clover_Nf2_Isochemical.h"



const std::string Force_F_Clover_Nf2_Isochemical::class_name = "Force_F_Clover_Nf2_Isochemical";

//====================================================================
void Force_F_Clover_Nf2_Isochemical::set_parameters(const Parameters& params)
{
  const string str_vlevel = params.get_string("verbose_level");

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
void Force_F_Clover_Nf2_Isochemical::set_parameters(const double kappa, const double cSW,
                                                    const double mu, const std::vector<int> bc)
{
  int Ndim = CommonParameters::Ndim();

  //- print input parameters
  vout.general(m_vl, "%s:\n", class_name.c_str());
  vout.general(m_vl, "  kappa = %12.8f\n", kappa);
  vout.general(m_vl, "  cSW   = %12.8f\n", cSW);
  vout.general(m_vl, "  mu    = %12.8f\n", mu);
  for (int imu = 0; imu < Ndim; ++imu) {
    vout.general(m_vl, "  boundary[%d] = %2d\n", imu, bc[imu]);
  }

  //- range check
  // NB. kappa,cSW,mu == 0 is allowed.
  assert(bc.size() == Ndim);

  //- store values
  m_kappa = kappa;
  m_cSW   = cSW;
  m_mu    = mu;

  m_boundary.resize(Ndim);
  for (int dir = 0; dir < Ndim; ++dir) {
    m_boundary[dir] = bc[dir];
  }

  //- propagate parameters
  m_fopr_c->set_parameters(m_kappa, m_cSW, m_mu, m_boundary);

  m_force_w->set_parameters(m_kappa, m_mu, m_boundary);
  m_force_csw->set_parameters(m_kappa, m_cSW, m_boundary);
}


//====================================================================
void Force_F_Clover_Nf2_Isochemical::init(std::string repr)
{
  m_repr = repr;

  m_Ndim = CommonParameters::Ndim();

  m_fopr_c    = new Fopr_Clover_Isochemical(m_repr);
  m_force_w   = new Force_F_Wilson_Nf2_Isochemical(m_repr);
  m_force_csw = new Force_F_CloverTerm(m_repr);

  m_boundary.resize(CommonParameters::Ndim());
}


//====================================================================
void Force_F_Clover_Nf2_Isochemical::tidyup()
{
  delete m_force_csw;
  delete m_force_w;
  delete m_fopr_c;
}


//====================================================================
void Force_F_Clover_Nf2_Isochemical::force_udiv(Field& force_, const Field& eta_)
{
  //int Nc   = CommonParameters::Nc();
  int Nvol = CommonParameters::Nvol();
  int Ndim = CommonParameters::Ndim();

  Field_G force(Nvol, Ndim);
  Field_F zeta(Nvol, 1);
  Field_F eta(eta_);

  m_fopr_c->H(zeta, eta);

  set_mode("H");
  force_udiv1_impl(force, zeta, eta);
  copy(force_, force); // force_ = force;

  set_mode("Hdag");
  force_udiv1_impl(force, eta, zeta);
  axpy(force_, 1.0, force); // force_ += force;
}


//====================================================================
void Force_F_Clover_Nf2_Isochemical::force_udiv1(Field& force_, const Field& zeta_, const Field& eta_)
{
  //int Nc   = CommonParameters::Nc();
  int Nvol = CommonParameters::Nvol();
  int Ndim = CommonParameters::Ndim();

  Field_G force(Nvol, Ndim);
  Field_F zeta(zeta_);
  Field_F eta(eta_);

  force_udiv1_impl(force, zeta, eta);

  copy(force_, force); // force_ = force;
}


//====================================================================
void Force_F_Clover_Nf2_Isochemical::force_udiv1_impl(Field_G& force, const Field_F& zeta, const Field_F& eta)
{
  //int Nc   = CommonParameters::Nc();
  //int Nd   = CommonParameters::Nd();
  int Nvol = CommonParameters::Nvol();
  int Ndim = CommonParameters::Ndim();

  ShiftField_lex shift;
  Field_G        force2(Nvol, Ndim);

  force.set(0.0);

  m_force_w->set_mode(m_mode);
  m_force_w->force_udiv1(force, zeta, eta);

  m_force_csw->force_udiv1(force2, zeta, eta);

  axpy(force, 1.0, force2); // force += force2;
}


//====================================================================
//============================================================END=====
