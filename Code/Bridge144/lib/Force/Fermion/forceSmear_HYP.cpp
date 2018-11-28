#include "BridgeLib_Private.h"

/*!
        @file    $Id:: forceSmear_HYP.cpp #$

        @brief

        @author  Hideo Matsufuru  (matsufuru)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 00:49:46 #$

        @version $LastChangedRevision: 1561 $
*/

#include "forceSmear_HYP.h"



#ifdef USE_FACTORY
namespace {
  ForceSmear *create_object(Projection *proj)
  {
    return new ForceSmear_HYP(proj);
  }


  bool init = ForceSmear::Factory::Register("HYP", create_object);
}
#endif



const std::string ForceSmear_HYP::class_name = "ForceSmear_HYP";

//====================================================================
void ForceSmear_HYP::set_parameters(const Parameters& params)
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
void ForceSmear_HYP::set_parameters(double alpha1, double alpha2, double alpha3)
{
  //- print input parameters
  vout.general(m_vl, "Force of %s:\n", class_name.c_str());
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
void ForceSmear_HYP::init()
{
  m_Ndim = CommonParameters::Ndim();
  m_Nvol = CommonParameters::Nvol();

  m_U.resize(m_Ndim);

  m_v1.resize(size1());
  m_v2.resize(size2());

  m_Sigma3.resize(size2());
  m_Sigma2.resize(size1b());

  m_iTheta3.resize(m_Ndim);
  m_iTheta2.resize(size2());
  m_iTheta1.resize(size1b());
}


//====================================================================
void ForceSmear_HYP::force_udiv(Field_G& Sigma, const Field_G& Sigmap, const Field_G& U)
{
  assert(U.nvol() == m_Nvol);
  assert(U.nex() == m_Ndim);
  assert(Sigmap.nvol() == m_Nvol);
  assert(Sigmap.nex() == m_Ndim);

  for (int mu = 0; mu < m_Ndim; ++mu) {
    m_U[mu].setpart_ex(0, U, mu);
  }

  // vout.general(m_vl, " smearing step-1\n");
  smear_step1();
  // vout.general(m_vl, " smearing step-2\n");
  smear_step2();

  Sigma.set(0.0);

  // vout.general(m_vl, " smeared force step-3\n");
  force_step3(Sigma, Sigmap);

  // vout.general(m_vl, " smeared force step-2\n");
  force_step2(Sigma);

  // vout.general(m_vl, " smeared force step-1\n");
  force_step1(Sigma);

  // vout.general(m_vl, " smeared force finished\n");
}


//====================================================================
void ForceSmear_HYP::force_step3(Field_G& Sigma, const Field_G& Sigmap)
{
  Field_G Sigmap_tmp(m_Nvol, 1), C(m_Nvol, 1), c_tmp(m_Nvol, 1);
  Field_G Xi(m_Nvol, 1);

  for (int mu = 0; mu < m_Ndim; ++mu) {
    C.set(0.0);

    for (int nu = 0; nu < m_Ndim; ++nu) {
      if (nu == mu) continue;
      staple(c_tmp, m_v2[idx2(mu, nu)], m_v2[idx2(nu, mu)], mu, nu);
      C.addpart_ex(0, c_tmp, 0);
    }
    scal(C, m_alpha1 / 6.0); // C *= m_alpha1 / 6.0;

    Sigmap_tmp.setpart_ex(0, Sigmap, mu);
    m_proj->force_recursive(Xi, m_iTheta3[mu],
                            m_alpha1, Sigmap_tmp, C, m_U[mu]);
    Sigma.addpart_ex(mu, Xi, 0);
  }

  for (int mu = 0; mu < m_Ndim; ++mu) {
    for (int nu = 0; nu < m_Ndim; ++nu) {
      if (nu == mu) continue;

      force_each(m_Sigma3[idx2(mu, nu)],
                 m_v2[idx2(mu, nu)], m_v2[idx2(nu, mu)],
                 m_iTheta3[mu], m_iTheta3[nu], mu, nu);

      scal(m_Sigma3[idx2(mu, nu)], m_alpha1 / 6.0); // m_Sigma3[idx2(mu, nu)] *= m_alpha1 / 6.0;
    }
  }
}


//====================================================================
void ForceSmear_HYP::force_step2(Field_G& Sigma)
{
  Field_G C(m_Nvol, 1), c_tmp(m_Nvol, 1), Xi(m_Nvol, 1);

  for (int mu = 0; mu < m_Ndim; ++mu) {
    for (int nu = 0; nu < m_Ndim; ++nu) {
      if (nu == mu) continue;

      C.set(0.0);

      for (int rho = 0; rho < m_Ndim; ++rho) {
        if ((rho == mu) || (rho == nu)) continue;
        staple(c_tmp, m_v1[idx1(mu, nu, rho)],
               m_v1[idx1(rho, nu, mu)], mu, rho);
        C.addpart_ex(0, c_tmp, 0);
      }
      scal(C, m_alpha2 / 4.0); // C *= m_alpha2 / 4.0;

      m_proj->force_recursive(Xi, m_iTheta2[idx2(mu, nu)],
                              m_alpha2, m_Sigma3[idx2(mu, nu)], C, m_U[mu]);
      Sigma.addpart_ex(mu, Xi, 0);
    }
  }

  for (int mu = 0; mu < m_Ndim; ++mu) {
    for (int nu = 0; nu < m_Ndim; ++nu) {
      if (nu == mu) continue;
      for (int rho = 0; rho < m_Ndim; ++rho) {
        if ((rho == mu) || (rho == nu)) continue;
        force_each(m_Sigma2[idx1b(mu, nu, rho)],
                   m_v1[idx1(mu, nu, rho)], m_v1[idx1(rho, nu, mu)],
                   m_iTheta2[idx2(mu, nu)], m_iTheta2[idx2(rho, nu)], mu, rho);
        scal(m_Sigma2[idx1b(mu, nu, rho)], m_alpha2 / 4.0); // m_Sigma2[idx1b(mu, nu, rho)] *= m_alpha2 / 4.0;
      }
    }
  }
}


//====================================================================
void ForceSmear_HYP::force_step1(Field_G& Sigma)
{
  Field_G Sigma_tmp(m_Nvol, 1), C(m_Nvol, 1), Xi(m_Nvol, 1);

  for (int mu = 0; mu < m_Ndim; ++mu) {
    for (int nu = 0; nu < m_Ndim; ++nu) {
      if (nu == mu) continue;
      for (int rho = 0; rho < m_Ndim; ++rho) {
        if ((rho == mu) || (rho == nu)) continue;

        int sig = 6 - mu - nu - rho;
        staple(C, m_U[mu], m_U[sig], mu, sig);

        scal(C, m_alpha3 / 2.0); // C *= m_alpha3 / 2.0;

        m_proj->force_recursive(Xi, m_iTheta1[idx1b(mu, nu, rho)],
                                m_alpha3, m_Sigma2[idx1b(mu, nu, rho)], C, m_U[mu]);
        Sigma.addpart_ex(mu, Xi, 0);
      }
    }
  }

  for (int mu = 0; mu < m_Ndim; ++mu) {
    for (int nu = 0; nu < m_Ndim; ++nu) {
      if (nu == mu) continue;
      for (int rho = 0; rho < m_Ndim; ++rho) {
        if ((rho == mu) || (rho == nu)) continue;
        int sig = 6 - mu - nu - rho;
        force_each(Sigma_tmp, m_U[mu], m_U[sig],
                   m_iTheta1[idx1b(mu, nu, rho)], m_iTheta1[idx1b(sig, nu, rho)],
                   mu, sig);
        scal(Sigma_tmp, m_alpha3 / 2.0); // Sigma_tmp *= m_alpha3 / 2.0;
        Sigma.addpart_ex(mu, Sigma_tmp, 0);
      }
    }
  }
}


//====================================================================
void ForceSmear_HYP::force_each(Field_G& Sigma_mu,
                                const Field_G& V_mu, const Field_G& V_nu,
                                const Field_G& iTheta_mu,
                                const Field_G& iTheta_nu,
                                int mu, int nu)
{
  Field_G vt1(m_Nvol, 1), vt2(m_Nvol, 1), vt3(m_Nvol, 1);

  Sigma_mu.set(0.0);

  m_shift.backward(vt1, V_nu, mu);
  m_shift.backward(vt2, V_mu, nu);
  mult_Field_Gnd(vt3, 0, vt1, 0, vt2, 0);
  multadd_Field_Gnd(Sigma_mu, 0, vt3, 0, iTheta_nu, 0, 1.0);

  mult_Field_Gdn(vt3, 0, iTheta_mu, 0, V_nu, 0);
  mult_Field_Gdn(vt2, 0, vt1, 0, vt3, 0);
  m_shift.forward(vt3, vt2, nu);
  axpy(Sigma_mu, 1.0, vt3); // Sigma_mu += vt3;

  mult_Field_Gdn(vt3, 0, V_mu, 0, iTheta_nu, 0);
  mult_Field_Gdn(vt2, 0, vt1, 0, vt3, 0);
  m_shift.forward(vt3, vt2, nu);
  axpy(Sigma_mu, 1.0, vt3); // Sigma_mu += vt3;

  m_shift.backward(vt1, iTheta_nu, mu);
  m_shift.backward(vt2, V_mu, nu);
  mult_Field_Gnd(vt3, 0, vt1, 0, vt2, 0);
  multadd_Field_Gnd(Sigma_mu, 0, vt3, 0, V_nu, 0, 1.0);

  mult_Field_Gdd(vt2, 0, vt1, 0, V_mu, 0);
  mult_Field_Gnn(vt3, 0, vt2, 0, V_nu, 0);
  m_shift.forward(vt2, vt3, nu);
  axpy(Sigma_mu, 1.0, vt2); // Sigma_mu += vt2;

  m_shift.backward(vt1, V_nu, mu);
  m_shift.backward(vt2, iTheta_mu, nu);
  mult_Field_Gnd(vt3, 0, vt1, 0, vt2, 0);
  multadd_Field_Gnd(Sigma_mu, 0, vt3, 0, V_nu, 0, 1.0);
}


//====================================================================
void ForceSmear_HYP::smear_step1()
{
  Field_G c1(m_Nvol, 1);

  for (int mu = 0; mu < m_Ndim; ++mu) {
    for (int nu = 0; nu < m_Ndim; ++nu) {
      if (nu == mu) continue;
      for (int rho = nu + 1; rho < m_Ndim; ++rho) {
        if (rho == mu) continue;
        int sig = 6 - mu - nu - rho;
        staple(c1, m_U[mu], m_U[sig], mu, sig);
        scal(c1, m_alpha3 / 2.0); // c1 *= m_alpha3 / 2.0;
        m_proj->project(m_v1[idx1(mu, nu, rho)], m_alpha3, c1, m_U[mu]);
      }
    }
  }
}


//====================================================================
void ForceSmear_HYP::smear_step2()
{
  Field_G c2(m_Nvol, 1), u_tmp(m_Nvol, 1);

  for (int mu = 0; mu < m_Ndim; ++mu) {
    for (int nu = 0; nu < m_Ndim; ++nu) {
      if (nu == mu) continue;

      c2.set(0.0);

      for (int rho = 0; rho < m_Ndim; ++rho) {
        if ((rho != mu) && (rho != nu)) {
          staple(u_tmp,
                 m_v1[idx1(mu, nu, rho)], m_v1[idx1(rho, nu, mu)], mu, rho);
          c2.addpart_ex(0, u_tmp, 0);
        }
      }
      scal(c2, m_alpha2 / 4.0); // c2 *= m_alpha2 / 4.0;
      m_proj->project(m_v2[idx2(mu, nu)], m_alpha2, c2, m_U[mu]);
    }
  }
}


//====================================================================
void ForceSmear_HYP::staple(Field_G& c,
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
