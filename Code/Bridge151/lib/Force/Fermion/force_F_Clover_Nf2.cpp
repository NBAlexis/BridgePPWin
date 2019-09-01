/*!
        @file    force_F_Clover_Nf2.cpp

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/
#include "BridgeLib_Private.h"
#include "force_F_Clover_Nf2.h"

const std::string Force_F_Clover_Nf2::class_name = "Force_F_Clover_Nf2";

//====================================================================
void Force_F_Clover_Nf2::set_parameters(const Parameters& params)
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
void Force_F_Clover_Nf2::set_parameters(const double kappa, const double cSW,
                                        const std::vector<int> bc)
{
  const int Ndim = CommonParameters::Ndim();

  //- print input parameters
  vout.general(m_vl, "%s:\n", class_name.c_str());
  vout.general(m_vl, "  kappa = %12.8f\n", kappa);
  vout.general(m_vl, "  cSW   = %12.8f\n", cSW);
  for (int mu = 0; mu < Ndim; ++mu) {
    vout.general(m_vl, "  boundary[%d] = %2d\n", mu, bc[mu]);
  }

  //- range check
  // NB. kappa,cSW == 0 is allowed.
  assert(bc.size() == Ndim);

  //- store values
  m_kappa = kappa;
  m_cSW   = cSW;

  m_boundary.resize(Ndim);
  m_boundary = bc;

  //- propagate parameters
  m_fopr_c->set_parameters(m_kappa, m_cSW, m_boundary);

  m_force_w->set_parameters(m_kappa, m_boundary);
  m_force_csw->set_parameters(m_kappa, m_cSW, m_boundary);
}


//====================================================================
void Force_F_Clover_Nf2::init(std::string repr)
{
  m_repr = repr;

  m_fopr_c    = new Fopr_Clover(repr);
  m_force_w   = new Force_F_Wilson_Nf2(repr);
  m_force_csw = new Force_F_CloverTerm(repr);

  m_boundary.resize(CommonParameters::Ndim());

  const int Nvol = CommonParameters::Nvol();
  m_Ndim = CommonParameters::Ndim();

  m_Cud = new Field_G(Nvol, m_Ndim * m_Ndim);
}


//====================================================================
void Force_F_Clover_Nf2::tidyup()
{
  delete m_Cud;
  delete m_force_csw;
  delete m_force_w;
  delete m_fopr_c;
}


//====================================================================
void Force_F_Clover_Nf2::force_udiv(Field& force_, const Field& eta_)
{
  const int Nvol = CommonParameters::Nvol();
  const int Ndim = CommonParameters::Ndim();

  Field_F zeta;
  Field_F eta(eta_);

  m_fopr_c->H(zeta, eta);

  Field_G force1(Nvol, Ndim);
  force_udiv1(force1, eta, zeta);
  copy(force_, force1); // force_ = force1;

  force_udiv1(force1, zeta, eta);
  axpy(force_, 1.0, force1); // force_ += force1;
}


//====================================================================
void Force_F_Clover_Nf2::force_udiv1(Field& force_, const Field& zeta_, const Field& eta_)
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
void Force_F_Clover_Nf2::force_udiv1_impl(Field_G& force, const Field_F& zeta, const Field_F& eta)
{
  const int Nvol = CommonParameters::Nvol();
  const int Ndim = CommonParameters::Ndim();

  force.set(0.0);

  m_force_w->force_udiv1(force, zeta, eta);

  Field_G force2(Nvol, Ndim);
  m_force_csw->force_udiv1(force2, zeta, eta);

  axpy(force, 1.0, force2); // force += force2;
}


//====================================================================
void Force_F_Clover_Nf2::set_component()
{
  for (int mu = 0; mu < m_Ndim; ++mu) {
    for (int nu = 0; nu < m_Ndim; ++nu) {
      if (nu == mu) continue;

      Staple_lex staple;

      Field_G Cmu_ud1;
      staple.upper(Cmu_ud1, *m_U, mu, nu);

      Field_G Cmu_ud2;
      staple.lower(Cmu_ud2, *m_U, mu, nu);
      //Cmu_ud -= staple.lower(*m_U, mu, nu);

      axpy(Cmu_ud1, -1.0, Cmu_ud2);
      m_Cud->setpart_ex(index_dir(mu, nu), Cmu_ud1, 0);
    }
  }
}


//====================================================================
//============================================================END=====
