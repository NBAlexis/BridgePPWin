/*!
        @file    smear_APE.cpp

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/
#include "BridgeLib_Private.h"
#include "smear_APE.h"

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init = Smear_APE::register_factory();
}
#endif

const std::string Smear_APE::class_name = "Smear_APE";

//====================================================================
void Smear_APE::set_parameters(const Parameters& params)
{
  const string str_vlevel = params.get_string("verbose_level");

  m_vl = vout.set_verbose_level(str_vlevel);

  //- fetch and check input parameters
  double rho1;

  int err = 0;
  err += params.fetch_double("rho_uniform", rho1);

  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameter not found.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }

  set_parameters(rho1);
}


//====================================================================
void Smear_APE::set_parameters(const double rho1)
{
  //- print input parameters
  vout.general(m_vl, "%s:\n", class_name.c_str());
  vout.general(m_vl, "  rho = %8.4f\n", rho1);

  //- range check
  // NB. rho == 0 is allowed.

  //- store values
  // m_rho.resize(m_Ndim * m_Ndim);  // already resized in init.
  for (int mu = 0; mu < m_Ndim; ++mu) {
    for (int nu = 0; nu < m_Ndim; ++nu) {
      m_rho[mu + nu * m_Ndim] = rho1;
    }
  }
}


//====================================================================
void Smear_APE::set_parameters(const std::vector<double>& rho)
{
  //- print input parameters
  vout.general(m_vl, "%s:\n", class_name.c_str());
  for (int mu = 0; mu < m_Ndim; ++mu) {
    vout.general(m_vl, "  rho[%d] = %8.4f\n", mu, rho[mu]);
  }

  // range check
  // NB. rho == 0 is allowed.
  assert(rho.size() == m_Ndim * m_Ndim);

  //- store values
  // m_rho.resize(m_Ndim * m_Ndim);  // already resized in init.
  for (int mu = 0; mu < m_Ndim; ++mu) {
    for (int nu = 0; nu < m_Ndim; ++nu) {
      m_rho[mu + nu * m_Ndim] = rho[mu + nu * m_Ndim];
    }
  }
}


//====================================================================
void Smear_APE::smear(Field_G& Usmear, const Field_G& U)
{
  const int Nvol = CommonParameters::Nvol();

  assert(U.nvol() == Nvol);
  assert(U.nex() == m_Ndim);
  assert(Usmear.nvol() == Nvol);
  assert(Usmear.nex() == m_Ndim);

  Usmear.set(0.0);

  for (int mu = 0; mu < m_Ndim; ++mu) {
    Field_G c_tmp;
    c_tmp.set(0.0);

    Field_G u_tmp;
    u_tmp.setpart_ex(0, U, mu);

    Field_G u_tmp2;

    for (int nu = 0; nu < m_Ndim; ++nu) {
      if (nu != mu) {
        Staple_lex staple;

        double rho = m_rho[mu + m_Ndim * nu];
        staple.upper(u_tmp2, U, mu, nu);
        c_tmp.addpart_ex(0, u_tmp2, 0, rho);

        staple.lower(u_tmp2, U, mu, nu);
        c_tmp.addpart_ex(0, u_tmp2, 0, rho);
      }
    }

    double rho0 = m_rho[mu + m_Ndim * mu];
    m_proj->project(u_tmp2, rho0, c_tmp, u_tmp);
    Usmear.setpart_ex(mu, u_tmp2, 0);
  }
}


//====================================================================
//============================================================END=====
