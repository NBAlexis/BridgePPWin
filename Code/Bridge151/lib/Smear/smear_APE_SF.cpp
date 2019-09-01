/*!
        @file    smear_APE_SF.cpp

        @brief

        @author  Yusuke Taniguchi
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/
#include "BridgeLib_Private.h"
#include "smear_APE_SF.h"

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init = Smear_APE_SF::register_factory();
}
#endif

const std::string Smear_APE_SF::class_name = "Smear_APE_SF";

//====================================================================
void Smear_APE_SF::set_parameters(const Parameters& params)
{
  const string str_vlevel = params.get_string("verbose_level");

  m_vl = vout.set_verbose_level(str_vlevel);

  //- fetch and check input parameters
  double              rho1;
  std::vector<double> phi, phipr;

  int err = 0;
  err += params.fetch_double("rho_uniform", rho1);
  err += params.fetch_double_vector("phi", phi);
  err += params.fetch_double_vector("phipr", phipr);

  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameter not found.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }

  set_parameters(rho1, &phi[0], &phipr[0]);
}


//====================================================================
void Smear_APE_SF::set_parameters(const double rho1,
                                  double *phi, double *phipr)
{
  //- print input parameters
  vout.general(m_vl, "%s:\n", class_name.c_str());
  vout.general(m_vl, "  rho = %8.4f\n", rho1);

  vout.general(m_vl, "  phi1  = %12.6f\n", phi[0]);
  vout.general(m_vl, "  phi2  = %12.6f\n", phi[1]);
  vout.general(m_vl, "  phi3  = %12.6f\n", phi[2]);
  vout.general(m_vl, "  phipr1= %12.6f\n", phipr[0]);
  vout.general(m_vl, "  phipr2= %12.6f\n", phipr[1]);
  vout.general(m_vl, "  phipr3= %12.6f\n", phipr[2]);

  //- range check
  // NB. rho == 0 is allowed.
  // NB. phi,phipr == 0 is allowed.

  //- store values
  // m_rho.resize(m_Ndim * m_Ndim);  // already resized in init.
  for (int mu = 0; mu < m_Ndim; ++mu) {
    for (int nu = 0; nu < m_Ndim; ++nu) {
      m_rho[mu + nu * m_Ndim] = rho1;
    }
  }

  for (int i = 0; i < 3; ++i) {
    m_phi[i]   = phi[i];
    m_phipr[i] = phipr[i];
  }
}


//====================================================================
void Smear_APE_SF::set_parameters(const std::vector<double>& rho, double *phi, double *phipr)
{
  //- print input parameters
  vout.general(m_vl, "%s:\n", class_name.c_str());
  for (int mu = 0; mu < m_Ndim; ++mu) {
    vout.general(m_vl, "  rho[%d] = %8.4f\n", mu, rho[mu]);
  }

  vout.general(m_vl, "  phi1  = %12.6f\n", phi[0]);
  vout.general(m_vl, "  phi2  = %12.6f\n", phi[1]);
  vout.general(m_vl, "  phi3  = %12.6f\n", phi[2]);
  vout.general(m_vl, "  phipr1= %12.6f\n", phipr[0]);
  vout.general(m_vl, "  phipr2= %12.6f\n", phipr[1]);
  vout.general(m_vl, "  phipr3= %12.6f\n", phipr[2]);

  //- range check
  // NB. rho == 0 is allowed.
  // NB. phi,phipr == 0 is allowed.
  assert(rho.size() == m_Ndim * m_Ndim);

  // store values
  // m_rho.resize(m_Ndim * m_Ndim);  // already resized in init.
  for (int mu = 0; mu < m_Ndim; ++mu) {
    for (int nu = 0; nu < m_Ndim; ++nu) {
      m_rho[mu + nu * m_Ndim] = rho[mu + nu * m_Ndim];
    }
  }

  for (int i = 0; i < 3; ++i) {
    m_phi[i]   = phi[i];
    m_phipr[i] = phipr[i];
  }
}


//====================================================================
void Smear_APE_SF::smear(Field_G& Usmear, const Field_G& U)
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

    Field_G_SF u_tmp2;

    for (int nu = 0; nu < m_Ndim; ++nu) {
      if (nu != mu) {
        Staple_SF staple;
        staple.set_parameters(m_phi, m_phipr);

        double rho = m_rho[mu + m_Ndim * nu];
        staple.upper(u_tmp2, U, mu, nu);
        c_tmp.addpart_ex(0, u_tmp2, 0, rho);

        staple.lower(u_tmp2, U, mu, nu);
        c_tmp.addpart_ex(0, u_tmp2, 0, rho);
      }
    }

    double rho0 = m_rho[mu + m_Ndim * mu];
    //    vout.general(m_vl,"mu=%d\n",mu);
    m_proj->project(u_tmp2, rho0, c_tmp, u_tmp);

    Field_G_SF set_wk(m_phi, m_phipr);
    if (mu != 3) set_wk.set_boundary_wk(u_tmp2);
    Usmear.setpart_ex(mu, u_tmp2, 0);

    /* For a debugging
    for(int site = 0; site < Nvol; ++site){
      vout.general(m_vl,"site smeared=%d\n",site);
      for(int ll=0; ll<9; ++ll){
    vout.general(m_vl,"(%lf %lf)\n",u_tmp2.cmp_r(ll,site,0),u_tmp2.cmp_i(ll,site,0));
      }
    }
    */
  }
}


//====================================================================
//============================================================END=====
