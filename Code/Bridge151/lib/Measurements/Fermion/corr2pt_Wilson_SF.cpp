/*!
        @file    corr2pt_Wilson_SF.cpp

        @brief

        @author  Yusuke Taniguchi
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/
#include "BridgeLib_Private.h"
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
  const int Lt = CommonParameters::Lt();

  const int Nvol = CommonParameters::Nvol();
  const int NPE  = CommonParameters::NPE();

  const double norm = (0.5 / Nvol / NPE) * Lt;

  std::vector<dcomplex> m_corr(Lt);

  //  double norm=0.5/Lvol*Lt*(2*0.130*2*0.130);
  //  vout.general(m_vl,"norm=%lf\n",norm);

  vout.general(m_vl, "AO correlator:\n");

  GammaMatrix qn_src  = m_gmset->get_GM(m_gmset->GAMMA5);
  GammaMatrix qn_sink = m_gmset->get_GM(m_gmset->GAMMA54);

  meson_corr(m_corr, qn_sink, qn_src, sq1, sq1);

  for (int t = 0; t < m_corr.size(); ++t) {
    vout.general(m_vl, "fA  %4d  %20.12e  %20.12e\n",
                 t, -norm * real(m_corr[t]), -norm * imag(m_corr[t]));
  }


  vout.general(m_vl, "PO correlator:\n");

  qn_src  = m_gmset->get_GM(m_gmset->GAMMA5);
  qn_sink = m_gmset->get_GM(m_gmset->GAMMA5);
  //  qn_src  = m_gmset->get_GM(m_gmset->UNITY);
  //  qn_sink = m_gmset->get_GM(m_gmset->UNITY);

  meson_corr(m_corr, qn_sink, qn_src, sq1, sq1);

  for (int t = 0; t < m_corr.size(); ++t) {
    vout.general(m_vl, "fP  %4d  %20.12e  %20.12e\n",
                 t, norm * real(m_corr[t]), norm * imag(m_corr[t]));
  }


  vout.general(m_vl, "AO' correlator:\n");

  qn_src  = m_gmset->get_GM(m_gmset->GAMMA5);
  qn_sink = m_gmset->get_GM(m_gmset->GAMMA54);

  meson_corr(m_corr, qn_sink, qn_src, sq2, sq2);

  for (int t = 0; t < m_corr.size(); ++t) {
    vout.general(m_vl, "fA'  %4d  %20.12e  %20.12e\n",
                 t, norm * real(m_corr[Lt - 1 - t]), norm * imag(m_corr[Lt - 1 - t]));
  }


  vout.general(m_vl, "PO' correlator:\n");

  qn_src  = m_gmset->get_GM(m_gmset->GAMMA5);
  qn_sink = m_gmset->get_GM(m_gmset->GAMMA5);
  //  qn_src  = m_gmset->get_GM(m_gmset->UNITY);
  //  qn_sink = m_gmset->get_GM(m_gmset->UNITY);

  meson_corr(m_corr, qn_sink, qn_src, sq2, sq2);

  for (int t = 0; t < m_corr.size(); ++t) {
    vout.general(m_vl, "fP'  %4d  %20.12e  %20.12e\n",
                 t, norm * real(m_corr[Lt - 1 - t]), norm * imag(m_corr[Lt - 1 - t]));
  }


  //- NB. use meson correlator at t=1 for a non-trivial check.
  // double result = real(m_corr[0]);
  const double result = real(m_corr[1]);

  return result;
}


//====================================================================
double Corr2pt_Wilson_SF::meson_corr(std::vector<dcomplex>& meson,
                                     const GammaMatrix& qn_sink,
                                     const GammaMatrix& qn_src,
                                     const std::vector<Field_F>& sq1,
                                     const std::vector<Field_F>& sq2)
{
  const int Nc   = CommonParameters::Nc();
  const int Nd   = CommonParameters::Nd();
  const int Nvol = CommonParameters::Nvol();
  const int Lt   = CommonParameters::Lt();
  const int Nx   = CommonParameters::Nx();
  const int Ny   = CommonParameters::Ny();
  const int Nz   = CommonParameters::Nz();
  const int Nt   = CommonParameters::Nt();

  assert(meson.size() == Lt);

  const GammaMatrix gm5     = m_gmset->get_GM(m_gmset->GAMMA5);
  const GammaMatrix gm_src  = qn_src.mult(gm5);
  const GammaMatrix gm_sink = gm5.mult(qn_sink);

  const Index_lex index;

  std::valarray<double> corr_r(Nd), corr_i(Nd);
  //  std::vector<int>      s2(Nd);
  Field corrF(2, Nvol, 2);

  std::valarray<dcomplex> corr_local(Nt);
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


  const int ipe_t = Communicator::ipe(3);

  std::vector<dcomplex> corr_global(Lt);
  for (int t_global = 0; t_global < Lt; ++t_global) {
    corr_global[t_global] = cmplx(0.0, 0.0);
  }
  for (int t = 0; t < Nt; ++t) {
    int t_global = t + ipe_t * Nt;
    corr_global[t_global] = corr_local[t];
  }
  for (int t_global = 0; t_global < Lt; ++t_global) {
    double cr_r = real(corr_global[t_global]);
    double cr_i = imag(corr_global[t_global]);

    cr_r = Communicator::reduce_sum(cr_r);
    cr_i = Communicator::reduce_sum(cr_i);

    meson[t_global] = cmplx(cr_r, cr_i);
    //    vout.general(m_vl,"%d %lf %lf\n",t,crr,cri);
  }


  //- NB. use meson correlator at t=1 for a non-trivial check.
  // double result = real(m_corr[0]);
  const double result = real(meson[1]);


  return result;
}


//====================================================================
//============================================================END=====
