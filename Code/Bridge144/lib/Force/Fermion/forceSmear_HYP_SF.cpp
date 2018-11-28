#include "BridgeLib_Private.h"

/*!
        @file    $Id:: forceSmear_HYP_SF.cpp #$

        @brief

        @author  Yusuke Tanigchi  (tanigchi)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 00:49:46 #$

        @version $LastChangedRevision: 1561 $
*/

#include "forceSmear_HYP_SF.h"


#ifdef USE_FACTORY
namespace {
  ForceSmear *create_object(Projection *proj)
  {
    return new ForceSmear_HYP_SF(proj);
  }


  bool init = ForceSmear::Factory::Register("HYP_SF", create_object);
}
#endif



const std::string ForceSmear_HYP_SF::class_name = "ForceSmear_HYP_SF";

//====================================================================
void ForceSmear_HYP_SF::set_parameters(const Parameters& params)
{
  const string str_vlevel = params.get_string("verbose_level");

  m_vl = vout.set_verbose_level(str_vlevel);

  //- fetch and check input parameters
  double              alpha1, alpha2, alpha3;
  std::vector<double> phi, phipr;

  int err = 0;
  err += params.fetch_double("alpha1", alpha1);
  err += params.fetch_double("alpha2", alpha2);
  err += params.fetch_double("alpha3", alpha3);
  err += params.fetch_double_vector("phi", phi);
  err += params.fetch_double_vector("phipr", phipr);

  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameter not found.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }


  set_parameters(alpha1, alpha2, alpha3, &phi[0], &phipr[0]);
}


//====================================================================
void ForceSmear_HYP_SF::set_parameters(double alpha1, double alpha2, double alpha3,
                                       double *phi, double *phipr)
{
  //- print input parameters
  vout.general(m_vl, "%s:\n", class_name.c_str());
  vout.general(m_vl, "  alpha1 = %10.6F\n", alpha1);
  vout.general(m_vl, "  alpha2 = %10.6F\n", alpha2);
  vout.general(m_vl, "  alpha3 = %10.6F\n", alpha3);

  vout.general(m_vl, "  phi1  = %12.6f\n", phi[0]);
  vout.general(m_vl, "  phi2  = %12.6f\n", phi[1]);
  vout.general(m_vl, "  phi3  = %12.6f\n", phi[2]);
  vout.general(m_vl, "  phipr1= %12.6f\n", phipr[0]);
  vout.general(m_vl, "  phipr2= %12.6f\n", phipr[1]);
  vout.general(m_vl, "  phipr3= %12.6f\n", phipr[2]);

  //- range check
  // NB. alpha1,alpha2,alpha3 == 0 is allowed.
  // NB. phi,phipr == 0 is allowed.

  //- store values
  m_alpha1 = alpha1;
  m_alpha2 = alpha2;
  m_alpha3 = alpha3;

  for (int i = 0; i < 3; ++i) {
    m_phi[i]   = phi[i];
    m_phipr[i] = phipr[i];
  }

  //- propagate parameters
  set_wk.set_parameters(m_phi, m_phipr);
}


//====================================================================
void ForceSmear_HYP_SF::init()
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

/*
Field ForceSmear_HYP_SF::force_udiv(const Field_G& Sigmap, const Field_G& U)
{
  Field_G Sigma(Sigmap.nvol(), Sigmap.nex());

  force_udiv(Sigma, Sigmap, U);

  return Sigma;
}
*/

//====================================================================
void ForceSmear_HYP_SF::force_udiv(Field_G& Sigma, const Field_G& Sigmap, const Field_G& U)
{
  assert(U.nvol() == m_Nvol);
  assert(U.nex() == m_Ndim);
  assert(Sigmap.nvol() == m_Nvol);
  assert(Sigmap.nex() == m_Ndim);

  for (int mu = 0; mu < m_Ndim; ++mu) {
    m_U[mu].setpart_ex(0, U, mu);
    if (mu != 3) set_wk.set_boundary_wk(m_U[mu]);
  }

  // vout.general(m_vl," smearing step-1\n");
  smear_step1();
  // vout.general(m_vl," smearing step-2\n");
  smear_step2();

  Sigma.set(0.0);

  // vout.general(m_vl," smeared force step-3\n");
  force_step3(Sigma, Sigmap);

  // vout.general(m_vl," smeared force step-2\n");
  force_step2(Sigma);

  // vout.general(m_vl," smeared force step-1\n");
  force_step1(Sigma);

  // vout.general(m_vl," smeared force finished\n");
}


//====================================================================

/*!
<ul>
<li>See the implementation note "note_cloverHMC.pdf" (2 Apr 2012) by H.Matsufuru.
<li>Sigmap \f$=\Sigma_\mu'(x)\f$ in eq.(3.17) in terms of (k)-th smearing.
<li>m_U[mu] \f$=U_{\mu}(x)\f$ is (k-1)-th smeared 0-th level link.
<li>m_v2[idx2(mu,nu)] \f$=V_{\mu;\nu}^{(2)}(x)\f$ is (k-1)-th smeared 2nd level link.
<li>Xi is used as \f$\Xi_\mu^{(3)}(x)=\Sigma_\mu'(x)\exp(iQ_\mu^{(3)}(x))+C_\mu^{(3)\dagger} i\Lambda_\mu^{(3)}(x)\f$ and is to be added to Sigma for the (k-1)-th smeared force.
<li>m_iTheta3[mu] is \f$i\Lambda_\mu^{(3)}(x)U_\mu(x)\f$.
<li>The pre-force field m_Sigma3[idx2(mu,nu)] \f$=\Sigma_{\mu;\nu}^{(3)}(x)\f$ is set to zero for the boundary spatial link.
</ul>
*/
void ForceSmear_HYP_SF::force_step3(Field_G& Sigma, const Field_G& Sigmap)
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
      if (mu != 3) set_wk.set_boundary_zero(m_Sigma3[idx2(mu, nu)]);
    }
  }
}


//====================================================================

/*!
<ul>
<li>See the implementation note "note_cloverHMC.pdf" (2 Apr 2012) by H.Matsufuru.
<li>m_U[mu] \f$=U_{\mu}(x)\f$ is (k-1)-th smeared 0-th level link.
<li>m_v1[idx1(mu,nu,rho)] \f$=V_{\mu;\nu\rho}^{(1)}(x)\f$ is (k-1)-th smeared 1st level link.
<li>Xi is used as \f$\Xi_\mu^{(2)}(x)=\Sigma_{\mu;\nu}^{(3)}(x)\exp(iQ_{\mu;\nu}^{(2)}(x))+C_{\mu;\nu}^{(2)\dagger} i\Lambda_{\mu;\nu}^{(2)}(x)\f$ and is to be added to Sigma for the (k-1)-th smeared force.
<li>m_iTheta2[idx2(mu,nu)] is \f$i\Lambda_{\mu;\nu}^{(2)}(x)U_\mu(x)\f$.
<li>The pre-force field m_Sigma2[idx1b(mu,nu,rho)] \f$=\Sigma_{\mu;\nu\rho}^{(2)}(x)\f$ is set to zero for the boundary spatial link.
</ul>
*/
void ForceSmear_HYP_SF::force_step2(Field_G& Sigma)
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
        if (mu != 3) set_wk.set_boundary_zero(m_Sigma2[idx1b(mu, nu, rho)]);
      }
    }
  }
}


//====================================================================

/*!
<ul>
<li>See the implementation note "note_cloverHMC.pdf" (2 Apr 2012) by H.Matsufuru.
<li>m_U[mu] \f$=U_{\mu}(x)\f$ is (k-1)-th smeared 0-th level link.
<li>Xi is used as \f$\Xi_\mu^{(1)}(x)=\Sigma_{\mu;\nu\rho}^{(2)}(x)\exp(iQ_{\mu;\nu\rho}^{(1)}(x))+C_{\mu;\nu\rho}^{(1)\dagger}i\Lambda_{\mu;\nu\rho}^{(1)}(x)\f$ and is to be added to Sigma for the (k-1)-th smeared force.
<li>m_iTheta1[idx1b(mu,nu,rho)] is \f$i\Lambda_{\mu;\nu\rho}^{(1)}(x)U_\mu(x)\f$.
<li>The force field Sigma is set to zero for the boundary spatial link.
</ul>
*/
void ForceSmear_HYP_SF::force_step1(Field_G& Sigma)
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

  set_wk.set_boundary_spatial_link_zero(Sigma);
}


//====================================================================

/*!
  <ul>
  <li>Each block evaluates:
  <li>The 1st block: \f$-U_\nu(x+\mu)U_\mu^\dagger(x+\nu)U_\nu^\dagger(x)\Lambda_\nu(x)\f$
  <li>The 2nd block: \f$-U_\nu^\dagger(x+\mu-\nu)U_\mu^\dagger(x-\nu)\Lambda_\mu(x-\nu)U_\nu(x-\nu)\f$
  <li>The 3rd block: \f$U_\nu^\dagger(x+\mu-\nu)U_\mu^\dagger(x-\nu)\Lambda_\nu(x-\nu)U_\nu(x-\nu)\f$
  <li>The 4th block: \f$\Lambda_\nu(x+\mu)U_\nu(x+\mu)U_\mu^\dagger(x+\nu)U_\nu^\dagger(x)\f$
  <li>The 5th block: \f$-U_\nu^\dagger(x+\mu-\nu)\Lambda_\mu(x+\mu-\nu)U_\mu^\dagger(x-\nu)U_\nu(x-\nu)\f$
  <li>The 6th block: \f$-U_\nu(x+\mu)U_\mu^\dagger(x+\nu)\Lambda_\mu(x+\nu)U_\nu^\dagger(x)\f$
  </ul>
*/
void ForceSmear_HYP_SF::force_each(Field_G& Sigma_mu,
                                   const Field_G& V_mu, const Field_G& V_nu,
                                   const Field_G& iTheta_mu,
                                   const Field_G& iTheta_nu,
                                   int mu, int nu)
{
  Field_G vt1(m_Nvol, 1), vt2(m_Nvol, 1), vt3(m_Nvol, 1);

  Sigma_mu.set(0.0);

  //- The 1st block
  m_shift.backward(vt1, V_nu, mu);
  if (mu == 3) set_wk.set_boundary_wkpr(vt1);
  m_shift.backward(vt2, V_mu, nu);
  if (nu == 3) set_wk.set_boundary_wkpr(vt2);
  mult_Field_Gnd(vt3, 0, vt1, 0, vt2, 0);
  multadd_Field_Gnd(Sigma_mu, 0, vt3, 0, iTheta_nu, 0, 1.0);

  //- The 2nd block
  mult_Field_Gdn(vt3, 0, iTheta_mu, 0, V_nu, 0);
  mult_Field_Gdn(vt2, 0, vt1, 0, vt3, 0);
  m_shift.forward(vt3, vt2, nu);
  axpy(Sigma_mu, 1.0, vt3); // Sigma_mu += vt3;

  //- The 3rd block
  mult_Field_Gdn(vt3, 0, V_mu, 0, iTheta_nu, 0);
  mult_Field_Gdn(vt2, 0, vt1, 0, vt3, 0);
  m_shift.forward(vt3, vt2, nu);
  axpy(Sigma_mu, 1.0, vt3); // Sigma_mu += vt3;

  //- The 4th block
  m_shift.backward(vt1, iTheta_nu, mu);
  m_shift.backward(vt2, V_mu, nu);
  if (nu == 3) set_wk.set_boundary_wkpr(vt2);
  mult_Field_Gnd(vt3, 0, vt1, 0, vt2, 0);
  multadd_Field_Gnd(Sigma_mu, 0, vt3, 0, V_nu, 0, 1.0);

  //- The 5th block
  mult_Field_Gdd(vt2, 0, vt1, 0, V_mu, 0);
  mult_Field_Gnn(vt3, 0, vt2, 0, V_nu, 0);
  m_shift.forward(vt2, vt3, nu);
  axpy(Sigma_mu, 1.0, vt2); // Sigma_mu += vt2;

  //- The 6th block
  m_shift.backward(vt1, V_nu, mu);
  if (mu == 3) set_wk.set_boundary_wkpr(vt1);
  m_shift.backward(vt2, iTheta_mu, nu);
  mult_Field_Gnd(vt3, 0, vt1, 0, vt2, 0);
  multadd_Field_Gnd(Sigma_mu, 0, vt3, 0, V_nu, 0, 1.0);
}


//====================================================================
void ForceSmear_HYP_SF::smear_step1()
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
        if (mu != 3) set_wk.set_boundary_wk(m_v1[idx1(mu, nu, rho)]);
      }
    }
  }
}


//====================================================================
void ForceSmear_HYP_SF::smear_step2()
{
  Field_G c2(m_Nvol, 1), u_tmp(m_Nvol, 1);

  for (int mu = 0; mu < m_Ndim; ++mu) {
    for (int nu = 0; nu < m_Ndim; ++nu) {
      if (nu == mu) continue;

      c2.set(0.0);

      for (int rho = 0; rho < m_Ndim; ++rho) {
        if ((rho != mu) && (rho != nu)) {
          staple(u_tmp, m_v1[idx1(mu, nu, rho)],
                 m_v1[idx1(rho, nu, mu)], mu, rho);
          c2.addpart_ex(0, u_tmp, 0);
        }
      }
      scal(c2, m_alpha2 / 4.0); // c2 *= m_alpha2 / 4.0;
      m_proj->project(m_v2[idx2(mu, nu)], m_alpha2, c2, m_U[mu]);
      if (mu != 3) set_wk.set_boundary_wk(m_v2[idx2(mu, nu)]);
    }
  }
}


//====================================================================
void ForceSmear_HYP_SF::staple(Field_G& c,
                               const Field_G& u_mu, const Field_G& u_nu,
                               int mu, int nu)
{
  Field_G v1(m_Nvol, 1), v2(m_Nvol, 1);

  // upper direction
  m_shift.backward(v1, u_mu, nu);
  if (nu == 3) set_wk.set_boundary_wkpr(v1);
  mult_Field_Gnn(v2, 0, u_nu, 0, v1, 0);

  m_shift.backward(v1, u_nu, mu);
  if (mu == 3) set_wk.set_boundary_wkpr(v1);
  mult_Field_Gnd(c, 0, v2, 0, v1, 0);

  // lower direction
  m_shift.backward(v2, u_nu, mu);
  if (mu == 3) set_wk.set_boundary_wkpr(v2);
  mult_Field_Gnn(v1, 0, u_mu, 0, v2, 0);
  mult_Field_Gdn(v2, 0, u_nu, 0, v1, 0);
  m_shift.forward(v1, v2, nu);
  c.addpart_ex(0, v1, 0);
}


//====================================================================
//============================================================END=====
