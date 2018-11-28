#include "BridgeLib_Private.h"

/*!
        @file    $Id:: smear_HYP.cpp #$

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 00:49:46 #$

        @version $LastChangedRevision: 1561 $
*/

#include "smear_HYP.h"


#ifdef USE_FACTORY
namespace {
  Smear *create_object(Projection *proj)
  {
    return new Smear_HYP(proj);
  }


  bool init = Smear::Factory::Register("HYP", create_object);
}
#endif



const std::string Smear_HYP::class_name = "Smear_HYP";

//====================================================================
void Smear_HYP::set_parameters(const Parameters& params)
{
  const string str_vlevel = params.get_string("verbose_level");

  m_vl = vout.set_verbose_level(str_vlevel);

  //- fetch and check input parameters
  double alpha1, alpha2, alpha3;

  int err = 0;
  err += params.fetch_double("alpha1", alpha1);
  err += params.fetch_double("alpha2", alpha2);
  err += params.fetch_double("alpha3", alpha3);

  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameter not found.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }

  set_parameters(alpha1, alpha2, alpha3);
}


//====================================================================
void Smear_HYP::set_parameters(double alpha1, double alpha2, double alpha3)
{
  //- print input parameters
  vout.general(m_vl, "%s:\n", class_name.c_str());
  vout.general(m_vl, "  alpha1 = %10.6F\n", alpha1);
  vout.general(m_vl, "  alpha2 = %10.6F\n", alpha2);
  vout.general(m_vl, "  alpha3 = %10.6F\n", alpha3);

  //- range check
  // NB. alpha1,alpha2,alpha3 == 0 is allowed.

  //- store values
  m_alpha1 = alpha1;
  m_alpha2 = alpha2;
  m_alpha3 = alpha3;
}


//====================================================================
void Smear_HYP::init()
{
  m_Ndim = CommonParameters::Ndim();
  m_Nvol = CommonParameters::Nvol();

  m_U.resize(m_Ndim);
  m_v1.resize(size_v1());
  m_v2.resize(size_v2());
}


//====================================================================
void Smear_HYP::smear(Field_G& Usmear, const Field_G& U)
{
  assert(U.nvol() == m_Nvol);
  assert(U.nex() == m_Ndim);

  assert(Usmear.nvol() == m_Nvol);
  assert(Usmear.nex() == m_Ndim);

  for (int mu = 0; mu < m_Ndim; ++mu) {
    m_U[mu].setpart_ex(0, U, mu);
  }

  step1();
  //  vout.general(m_vl, "level-1 step finished.\n");
  step2();
  //  vout.general(m_vl, "level-2 step finished.\n");
  step3(Usmear);
  //  vout.general(m_vl, "level-3 step finished.\n");
}


//====================================================================
void Smear_HYP::step1()
{
  Field_G c1(m_Nvol, 1);

  for (int mu = 0; mu < m_Ndim; ++mu) {
    for (int nu = 0; nu < m_Ndim; ++nu) {
      if (nu == mu) continue;

      for (int rho = nu + 1; rho < m_Ndim; ++rho) {
        if (rho == mu) continue;

        int sig = 6 - mu - nu - rho;
        staple(c1, m_U[mu], m_U[sig], mu, sig);
        //c1 *= m_alpha3 / 2.0;
        scal(c1, m_alpha3 * 0.5);
        m_proj->project(m_v1[index_v1(mu, nu, rho)], m_alpha3, c1, m_U[mu]);
      }
    }
  }
}


//====================================================================
void Smear_HYP::step2()
{
  Field_G c2(m_Nvol, 1), u_tmp(m_Nvol, 1);

  for (int mu = 0; mu < m_Ndim; ++mu) {
    for (int nu = 0; nu < m_Ndim; ++nu) {
      if (nu == mu) continue;
      c2.set(0.0);

      for (int rho = 0; rho < m_Ndim; ++rho) {
        if ((rho != mu) && (rho != nu)) {
          staple(u_tmp, m_v1[index_v1(mu, nu, rho)],
                 m_v1[index_v1(rho, nu, mu)], mu, rho);
          c2.addpart_ex(0, u_tmp, 0);
        }
      }

      //c2 *= m_alpha2 / 4.0;
      scal(c2, m_alpha2 * 0.25);
      m_proj->project(m_v2[index_v2(mu, nu)], m_alpha2, c2, m_U[mu]);
    }
  }
}


//====================================================================
void Smear_HYP::step3(Field_G& v)
{
  Field_G c3(m_Nvol, 1), u_tmp(m_Nvol, 1);

  for (int mu = 0; mu < m_Ndim; ++mu) {
    c3.set(0.0);

    for (int nu = 0; nu < m_Ndim; ++nu) {
      if (nu != mu) {
        staple(u_tmp, m_v2[index_v2(mu, nu)],
               m_v2[index_v2(nu, mu)], mu, nu);
        c3.addpart_ex(0, u_tmp, 0);
      }
    }

    //c3 *= m_alpha1 / 6.0;
    scal(c3, m_alpha1 / 6.0);
    m_proj->project(u_tmp, m_alpha1, c3, m_U[mu]);
    v.setpart_ex(mu, u_tmp, 0);
  }
}


//====================================================================
void Smear_HYP::staple(Field_G& c,
                       const Field_G& u_mu, const Field_G& u_nu,
                       int mu, int nu)
{
  Field_G v1(m_Nvol, 1), v2(m_Nvol, 1);

  //- upper direction
  m_shift.backward(v1, u_mu, nu);
  mult_Field_Gnn(v2, 0, u_nu, 0, v1, 0);

  m_shift.backward(v1, u_nu, mu);
  mult_Field_Gnd(c, 0, v2, 0, v1, 0);

  //- lower direction
  m_shift.backward(v2, u_nu, mu);
  mult_Field_Gnn(v1, 0, u_mu, 0, v2, 0);
  mult_Field_Gdn(v2, 0, u_nu, 0, v1, 0);
  m_shift.forward(v1, v2, nu);
  c.addpart_ex(0, v1, 0);
}


//====================================================================
//============================================================END=====
