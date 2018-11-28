#include "BridgeLib_Private.h"

/*!
        @file    $Id:: corr2pt_Wilson_SF.cpp #$

        @brief

        @author  Yusuke Taniguchi
                 $LastChangedBy: namekawa $

        @date    $LastChangedDate:: 2017-03-03 15:12:38 #$

        @version $LastChangedRevision: 1582 $
*/

#include "corr2pt_Wilson_SF.h"

//====================================================================

/*!
  The axial current:
  \f[
  A_\mu^a(x)={\overline{\psi}}(x)\gamma_\mu\gamma_5\frac{\tau^a}{2}\psi(x)
  \f]
  The pseudo scalar density:
  \f[
  P^a(x)={\overline{\psi}}(x)\gamma_5\frac{\tau^a}{2}\psi(x)
  \f]
  The boundary pseudo scalar density at t=0
  \f[
  {\cal O}^a=a^6\sum_{\vec{x},\vec{y}}\overline{\zeta}\left(\vec{x}\right)
  \gamma_5\frac{\tau^a}{2}\zeta\left(\vec{y}\right)
  \f]
  The boundary pseudo scalar density at t=T
  \f[
  {\cal O'}^a=a^6\sum_{\vec{x},\vec{y}}\overline{\zeta}'\left(\vec{x}\right)
  \gamma_5\frac{\tau^a}{2}\zeta'\left(\vec{y}\right)
  \f]
  The two point functions are
  \f[
  f_A(x_0)=-\frac{1}{N_f^2-1}\left\langle{A_0^a(x_0){\cal O}^a}\right\rangle
  =\frac{1}{2}\sum_{\vec{v},\vec{y}}{\rm tr}\left(
  \left\langle{\zeta(\vec{v}){\overline{\psi}}(x)}\right\rangle\gamma_0\gamma_5
  \left\langle{\psi(x)\overline{\zeta}(\vec{y})}\right\rangle\gamma_5\right)
  \f]
  \f[
  f_P(x_0)=-\frac{1}{N_f^2-1}\left\langle{P^a(x_0){\cal O}^a}\right\rangle
  =\frac{1}{2}\sum_{\vec{v},\vec{y}}{\rm tr}\left(
  \left\langle{\zeta(\vec{v}){\overline{\psi}}(x)}\right\rangle\gamma_5
  \left\langle{\psi(x)\overline{\zeta}(\vec{y})}\right\rangle\gamma_5\right)
  \f]
  \f[
  f_A'(x_0)=+\frac{1}{N_f^2-1}\left\langle{A_0^a(T-x_0){\cal O'}^a}\right\rangle
  =-\frac{1}{2}\sum_{\vec{v},\vec{y}}{\rm tr}\left(
  \left\langle{\zeta'(\vec{v}){\overline{\psi}}(T-x_0)}\right\rangle\gamma_0\gamma_5
  \left\langle{\psi(T-x_0)\overline{\zeta}'(\vec{y})}\right\rangle\gamma_5\right)
  \f]
  \f[
  f_P'(x_0)=-\frac{1}{N_f^2-1}\left\langle{P^a(T-x_0){\cal O'}^a}\right\rangle
  =\frac{1}{2}\sum_{\vec{v},\vec{y}}{\rm tr}\left(
  \left\langle{\zeta'(\vec{v}){\overline{\psi}}(T-x_0)}\right\rangle\gamma_5
  \left\langle{\psi(T-x_0)\overline{\zeta}'(\vec{y})}\right\rangle\gamma_5\right)
  \f]
*/

const std::string Corr2pt_Wilson_SF::class_name = "Corr2pt_Wilson_SF";

//====================================================================
double Corr2pt_Wilson_SF::fAfP(const std::vector<Field_F>& sq1,
                               const std::vector<Field_F>& sq2)
{
  int Lt = CommonParameters::Lt();

  std::vector<dcomplex> mcorr(Lt);
  GammaMatrix           qn_src, qn_sink;

  int    Lvol = CommonParameters::Lvol();
  //int    m_Nc = CommonParameters::Nc();
  double norm = 0.5 / Lvol * Lt;

  //  double norm=0.5/Lvol*Lt*(2*0.130*2*0.130);
  //  vout.general(m_vl,"norm=%lf\n",norm);

  vout.general(m_vl, "AO correlator:\n");

  qn_src  = m_gmset->get_GM(m_gmset->GAMMA5);
  qn_sink = m_gmset->get_GM(m_gmset->GAMMA54);

  meson_corr(mcorr, qn_sink, qn_src, sq1, sq1);

  for (int t = 0; t < mcorr.size(); ++t) {
    vout.general(m_vl, "fA  %4d  %20.12e  %20.12e\n",
                 t, -norm * real(mcorr[t]), -norm * imag(mcorr[t]));
  }


  vout.general(m_vl, "PO correlator:\n");

  qn_src  = m_gmset->get_GM(m_gmset->GAMMA5);
  qn_sink = m_gmset->get_GM(m_gmset->GAMMA5);
  //  qn_src  = m_gmset->get_GM(m_gmset->UNITY);
  //  qn_sink = m_gmset->get_GM(m_gmset->UNITY);

  meson_corr(mcorr, qn_sink, qn_src, sq1, sq1);

  for (int t = 0; t < mcorr.size(); ++t) {
    vout.general(m_vl, "fP  %4d  %20.12e  %20.12e\n",
                 t, norm * real(mcorr[t]), norm * imag(mcorr[t]));
  }


  vout.general(m_vl, "AO' correlator:\n");

  qn_src  = m_gmset->get_GM(m_gmset->GAMMA5);
  qn_sink = m_gmset->get_GM(m_gmset->GAMMA54);

  meson_corr(mcorr, qn_sink, qn_src, sq2, sq2);

  for (int t = 0; t < mcorr.size(); ++t) {
    vout.general(m_vl, "fA'  %4d  %20.12e  %20.12e\n",
                 t, norm * real(mcorr[Lt - 1 - t]), norm * imag(mcorr[Lt - 1 - t]));
  }


  vout.general(m_vl, "PO' correlator:\n");

  qn_src  = m_gmset->get_GM(m_gmset->GAMMA5);
  qn_sink = m_gmset->get_GM(m_gmset->GAMMA5);
  //  qn_src  = m_gmset->get_GM(m_gmset->UNITY);
  //  qn_sink = m_gmset->get_GM(m_gmset->UNITY);

  meson_corr(mcorr, qn_sink, qn_src, sq2, sq2);

  for (int t = 0; t < mcorr.size(); ++t) {
    vout.general(m_vl, "fP'  %4d  %20.12e  %20.12e\n",
                 t, norm * real(mcorr[Lt - 1 - t]), norm * imag(mcorr[Lt - 1 - t]));
  }


  //- NB. use meson correlator at t=1 for a non-trivial check.
  // double result = real(mcorr[0]);
  double result = real(mcorr[1]);

  return result;
}


//====================================================================
double Corr2pt_Wilson_SF::meson_corr(std::vector<dcomplex>& meson,
                                     const GammaMatrix& qn_sink,
                                     const GammaMatrix& qn_src,
                                     const std::vector<Field_F>& sq1,
                                     const std::vector<Field_F>& sq2)
{
  int Nc   = CommonParameters::Nc();
  int Nd   = CommonParameters::Nd();
  int Nvol = CommonParameters::Nvol();
  int Lt   = CommonParameters::Lt();
  int Nx   = CommonParameters::Nx();
  int Ny   = CommonParameters::Ny();
  int Nz   = CommonParameters::Nz();
  int Nt   = CommonParameters::Nt();

  assert(meson.size() == Lt);

  GammaMatrix gm_src, gm_sink, gm5;
  gm5     = m_gmset->get_GM(m_gmset->GAMMA5);
  gm_src  = qn_src.mult(gm5);
  gm_sink = gm5.mult(qn_sink);

  Index_lex          index;
  valarray<dcomplex> corr_local(Nt);
  valarray<double>   corr_r(Nd), corr_i(Nd);
//  std::vector<int>      s2(Nd);

  Field corrF(2, Nvol, 2);

  corr_local = cmplx(0.0, 0.0);

  for (int c0 = 0; c0 < Nc; ++c0) {
    for (int d0 = 0; d0 < Nd; ++d0) {
      int d1 = gm_src.index(d0);
      //      int d1 = qn_src.index(d0);

      for (int t = 0; t < Nt; ++t) {
        corr_r = 0.0;
        corr_i = 0.0;
        for (int z = 0; z < Nz; ++z) {
          for (int y = 0; y < Ny; ++y) {
            for (int x = 0; x < Nx; ++x) {
              int site = index.site(x, y, z, t);

              for (int s0 = 0; s0 < Nd; ++s0) {
                int s1 = gm_sink.index(s0);
                //int s1 = qn_sink.index(s0);

                for (int c1 = 0; c1 < Nc; ++c1) {
                  corr_r[s0] += sq1[c0 + Nc * d0].cmp_r(c1, s1, site)
                                * sq2[c0 + Nc * d1].cmp_r(c1, s0, site)
                                + sq1[c0 + Nc * d0].cmp_i(c1, s1, site)
                                * sq2[c0 + Nc * d1].cmp_i(c1, s0, site);

                  corr_i[s0] += sq1[c0 + Nc * d0].cmp_r(c1, s1, site)
                                * sq2[c0 + Nc * d1].cmp_i(c1, s0, site)
                                - sq1[c0 + Nc * d0].cmp_i(c1, s1, site)
                                * sq2[c0 + Nc * d1].cmp_r(c1, s0, site);
                }
              }
            }
          }
        }
        for (int s0 = 0; s0 < Nd; ++s0) {
          dcomplex gmf = gm_src.value(d0) * gm_sink.value(s0);
          //dcomplex gmf = qn_src.value(d0)*qn_sink.value(s0);
          dcomplex corr = cmplx(corr_r[s0], corr_i[s0]);
          corr_local[t] += gmf * corr;
          //corr_local[t] += corr;
        }
        //vout.general(m_vl,"%d %lf %lf\n",t,real(corr_local[t]),imag(corr_local[t]));
      }
    }
  }

  std::vector<dcomplex> corr_tmp(Lt);

  int ipet = Communicator::ipe(3);
  for (int t = 0; t < Lt; ++t) {
    corr_tmp[t] = cmplx(0.0, 0.0);
  }
  for (int t = 0; t < Nt; ++t) {
    int t2 = t + ipet * Nt;
    //    vout.general(m_vl,"%d %d\n",t,t2);
    corr_tmp[t2] = corr_local[t];
  }
  for (int t = 0; t < Lt; ++t) {
    double crr = real(corr_tmp[t]);
    double cri = imag(corr_tmp[t]);
    crr      = Communicator::reduce_sum(crr);
    cri      = Communicator::reduce_sum(cri);
    meson[t] = cmplx(crr, cri);
    //    vout.general(m_vl,"%d %lf %lf\n",t,crr,cri);
  }


  //- NB. use meson correlator at t=1 for a non-trivial check.
  // double result = real(mcorr[0]);
  double result = real(meson[1]);


  return result;
}


//====================================================================
//============================================================END=====
