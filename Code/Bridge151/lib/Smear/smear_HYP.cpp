/*!
        @file    smear_HYP.cpp

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/
#include "BridgeLib_Private.h"
#include "smear_HYP.h"

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init = Smear_HYP::register_factory();
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
void Smear_HYP::set_parameters(const double alpha1, const double alpha2, const double alpha3)
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
  for (int mu = 0; mu < m_Ndim; ++mu) {
    for (int nu = 0; nu < m_Ndim; ++nu) {
      if (nu == mu) continue;

      for (int rho = nu + 1; rho < m_Ndim; ++rho) {
        if (rho == mu) continue;

        int sig = 6 - mu - nu - rho;

        Field_G c1;
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
  for (int mu = 0; mu < m_Ndim; ++mu) {
    for (int nu = 0; nu < m_Ndim; ++nu) {
      if (nu == mu) continue;

      Field_G c2;
      c2.set(0.0);

      for (int rho = 0; rho < m_Ndim; ++rho) {
        if ((rho != mu) && (rho != nu)) {
          Field_G u_tmp;
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
  for (int mu = 0; mu < m_Ndim; ++mu) {
    Field_G c3;
    c3.set(0.0);

    Field_G u_tmp;

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
                       const int mu, const int nu)
{
  //- upper direction
  Field_G v1;

  m_shift.backward(v1, u_mu, nu);

  Field_G v2;
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
