/*!
        @file    corr2pt_4spinor.cpp

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/
#include "BridgeLib_Private.h"
#include "corr2pt_4spinor.h"

#include "Tools/epsilonTensor.h"

const std::string Corr2pt_4spinor::class_name = "Corr2pt_4spinor";

//====================================================================
void Corr2pt_4spinor::set_parameters(const Parameters& params)
{
  m_filename_output = params.get_string("filename_output");
  if (m_filename_output.empty()) {
    m_filename_output = "stdout";
  }

  const string str_vlevel = params.get_string("verbose_level");
  m_vl = vout.set_verbose_level(str_vlevel);
}


//====================================================================
void Corr2pt_4spinor::init()
{
  assert(CommonParameters::Nc() == 3);

  m_filename_output = "stdout";
}


//====================================================================
double Corr2pt_4spinor::meson_all(const std::vector<Field_F>& sq1,
                                  const std::vector<Field_F>& sq2)
{
  const int Lt = CommonParameters::Lt();

  std::ostream& log_file_previous = vout.getStream();
  std::ofstream log_file;

  if (m_filename_output != "stdout") {
    log_file.open(m_filename_output.c_str(), std::ios::app);
    vout.init(log_file);
  }

  std::vector<dcomplex> corr(Lt);

  GammaMatrix gm_src  = m_gmset->get_GM(m_gmset->GAMMA5);
  GammaMatrix gm_sink = m_gmset->get_GM(m_gmset->GAMMA5);
  vout.general(m_vl, "PS <-- PS correlator:\n");
  meson_correlator(corr, gm_sink, gm_src, sq1, sq2);
  for (int t = 0; t < corr.size(); ++t) {
    vout.general(m_vl, "  %4d  %20.12e  %20.12e\n",
                 t, real(corr[t]), imag(corr[t]));
  }
  const double result = real(corr[0]);

  gm_src  = m_gmset->get_GM(m_gmset->GAMMA1);
  gm_sink = m_gmset->get_GM(m_gmset->GAMMA1);
  vout.general(m_vl, "V1 <-- V1 correlator:\n");
  meson_correlator(corr, gm_sink, gm_src, sq1, sq2);
  for (int t = 0; t < corr.size(); ++t) {
    vout.general(m_vl, "  %4d  %20.12e  %20.12e\n",
                 t, real(corr[t]), imag(corr[t]));
  }

  gm_src  = m_gmset->get_GM(m_gmset->GAMMA2);
  gm_sink = m_gmset->get_GM(m_gmset->GAMMA2);
  vout.general(m_vl, "V2 <-- V2 correlator:\n");
  meson_correlator(corr, gm_sink, gm_src, sq1, sq2);
  for (int t = 0; t < corr.size(); ++t) {
    vout.general(m_vl, "  %4d  %20.12e  %20.12e\n",
                 t, real(corr[t]), imag(corr[t]));
  }

  gm_src  = m_gmset->get_GM(m_gmset->GAMMA3);
  gm_sink = m_gmset->get_GM(m_gmset->GAMMA3);
  vout.general(m_vl, "V3 <-- V3 correlator:\n");
  meson_correlator(corr, gm_sink, gm_src, sq1, sq2);
  for (int t = 0; t < corr.size(); ++t) {
    vout.general(m_vl, "  %4d  %20.12e  %20.12e\n",
                 t, real(corr[t]), imag(corr[t]));
  }

  gm_src  = m_gmset->get_GM(m_gmset->GAMMA5);
  gm_sink = m_gmset->get_GM(m_gmset->GAMMA45);
  vout.general(m_vl, "A4 <-- PS correlator:\n");
  meson_correlator(corr, gm_sink, gm_src, sq1, sq2);
  for (int t = 0; t < corr.size(); ++t) {
    vout.general(m_vl, "  %4d  %20.12e  %20.12e\n",
                 t, real(corr[t]), imag(corr[t]));
  }

  gm_src  = m_gmset->get_GM(m_gmset->GAMMA54);
  gm_sink = m_gmset->get_GM(m_gmset->GAMMA5);
  vout.general(m_vl, "PS <-- A4 correlator:\n");
  meson_correlator(corr, gm_sink, gm_src, sq1, sq2);
  for (int t = 0; t < corr.size(); ++t) {
    vout.general(m_vl, "  %4d  %20.12e  %20.12e\n",
                 t, real(corr[t]), imag(corr[t]));
  }

  gm_src  = m_gmset->get_GM(m_gmset->UNITY);
  gm_sink = m_gmset->get_GM(m_gmset->UNITY);
  vout.general(m_vl, "S <-- S correlator:\n");
  meson_correlator(corr, gm_sink, gm_src, sq1, sq2);
  for (int t = 0; t < corr.size(); ++t) {
    vout.general(m_vl, "  %4d  %20.12e  %20.12e\n",
                 t, real(corr[t]), imag(corr[t]));
  }

  gm_src  = m_gmset->get_GM(m_gmset->GAMMA51);
  gm_sink = m_gmset->get_GM(m_gmset->GAMMA51);
  vout.general(m_vl, "GAMMA51 <-- GAMMA51 correlator:\n");
  meson_correlator(corr, gm_sink, gm_src, sq1, sq2);
  for (int t = 0; t < corr.size(); ++t) {
    vout.general(m_vl, "  %4d  %20.12e  %20.12e\n",
                 t, real(corr[t]), imag(corr[t]));
  }

  gm_src  = m_gmset->get_GM(m_gmset->GAMMA52);
  gm_sink = m_gmset->get_GM(m_gmset->GAMMA52);
  vout.general(m_vl, "GAMMA52 <-- GAMMA52 correlator:\n");
  meson_correlator(corr, gm_sink, gm_src, sq1, sq2);
  for (int t = 0; t < corr.size(); ++t) {
    vout.general(m_vl, "  %4d  %20.12e  %20.12e\n",
                 t, real(corr[t]), imag(corr[t]));
  }

  gm_src  = m_gmset->get_GM(m_gmset->GAMMA53);
  gm_sink = m_gmset->get_GM(m_gmset->GAMMA53);
  vout.general(m_vl, "GAMMA53 <-- GAMMA53 correlator:\n");
  meson_correlator(corr, gm_sink, gm_src, sq1, sq2);
  for (int t = 0; t < corr.size(); ++t) {
    vout.general(m_vl, "  %4d  %20.12e  %20.12e\n",
                 t, real(corr[t]), imag(corr[t]));
  }

  gm_src  = m_gmset->get_GM(m_gmset->SIGMA12);
  gm_sink = m_gmset->get_GM(m_gmset->SIGMA12);
  vout.general(m_vl, "SIGMA12 <-- SIGMA12 correlator:\n");
  meson_correlator(corr, gm_sink, gm_src, sq1, sq2);
  for (int t = 0; t < corr.size(); ++t) {
    vout.general(m_vl, "  %4d  %20.12e  %20.12e\n",
                 t, real(corr[t]), imag(corr[t]));
  }

  gm_src  = m_gmset->get_GM(m_gmset->SIGMA23);
  gm_sink = m_gmset->get_GM(m_gmset->SIGMA23);
  vout.general(m_vl, "SIGMA23 <-- SIGMA23 correlator:\n");
  meson_correlator(corr, gm_sink, gm_src, sq1, sq2);
  for (int t = 0; t < corr.size(); ++t) {
    vout.general(m_vl, "  %4d  %20.12e  %20.12e\n",
                 t, real(corr[t]), imag(corr[t]));
  }

  gm_src  = m_gmset->get_GM(m_gmset->SIGMA31);
  gm_sink = m_gmset->get_GM(m_gmset->SIGMA31);
  vout.general(m_vl, "SIGMA31 <-- SIGMA31 correlator:\n");
  meson_correlator(corr, gm_sink, gm_src, sq1, sq2);
  for (int t = 0; t < corr.size(); ++t) {
    vout.general(m_vl, "  %4d  %20.12e  %20.12e\n",
                 t, real(corr[t]), imag(corr[t]));
  }


  if (m_filename_output != "stdout") {
    log_file.close();
    vout.init(log_file_previous);
  }

  return result;
}


//====================================================================
void Corr2pt_4spinor::meson_correlator(std::vector<dcomplex>& corr_global,
                                       const GammaMatrix& gm_sink,
                                       const GammaMatrix& gm_src,
                                       const std::vector<Field_F>& sq1,
                                       const std::vector<Field_F>& sq2)
{
  const int Nc = CommonParameters::Nc();
  const int Nd = CommonParameters::Nd();
  const int Lt = CommonParameters::Lt();
  const int Nt = CommonParameters::Nt();

  assert(corr_global.size() == Lt);

  const GammaMatrix gm5         = m_gmset->get_GM(m_gmset->GAMMA5);
  const GammaMatrix gm_gm5_src  = gm_src.mult(gm5);
  const GammaMatrix gm5_gm_sink = gm5.mult(gm_sink);

  std::vector<dcomplex> corr_local(Nt, cmplx(0.0, 0.0));
  for (int c0 = 0; c0 < Nc; ++c0) {
    for (int d0 = 0; d0 < Nd; ++d0) {
      int d1 = gm_gm5_src.index(d0);

      for (int t = 0; t < Nt; ++t) {
        dcomplex corr_t;

        contract_at_t(corr_t, gm5_gm_sink,
                      sq1[c0 + Nc * d0], sq2[c0 + Nc * d1], t);

        corr_local[t] += gm_gm5_src.value(d0) * corr_t;
      }
    }
  }
  global_corr_t(corr_global, corr_local);
}


//====================================================================
double Corr2pt_4spinor::meson_momentum_all(const std::vector<Field_F>& sq1,
                                           const std::vector<Field_F>& sq2,
                                           const std::vector<int>& source_position)
{
  const int Ndim = CommonParameters::Ndim();
  const int Lt   = CommonParameters::Lt();

  const int N_momentum = 10;

  typedef std::vector<int> MomentumSet;
  std::vector<MomentumSet> momentum_sink(N_momentum);
  for (int i_momentum = 0; i_momentum < N_momentum; i_momentum++) {
    momentum_sink[i_momentum].resize(Ndim - 1);
  }

  //- momentum_sink[0] = (1,0,0)
  int i_momentum = 0;
  momentum_sink[i_momentum][0] = 1;
  momentum_sink[i_momentum][1] = 0;
  momentum_sink[i_momentum][2] = 0;

  //- momentum_sink[1] = (0,1,0)
  i_momentum = 1;
  momentum_sink[i_momentum][0] = 0;
  momentum_sink[i_momentum][1] = 1;
  momentum_sink[i_momentum][2] = 0;

  //- momentum_sink[2] = (0,0,1)
  i_momentum = 2;
  momentum_sink[i_momentum][0] = 0;
  momentum_sink[i_momentum][1] = 0;
  momentum_sink[i_momentum][2] = 1;

  //- momentum_sink[3] = (1,1,0)
  i_momentum = 3;
  momentum_sink[i_momentum][0] = 1;
  momentum_sink[i_momentum][1] = 1;
  momentum_sink[i_momentum][2] = 0;

  //- momentum_sink[4] = (0,1,1)
  i_momentum = 4;
  momentum_sink[i_momentum][0] = 0;
  momentum_sink[i_momentum][1] = 1;
  momentum_sink[i_momentum][2] = 1;

  //- momentum_sink[5] = (1,0,1)
  i_momentum = 5;
  momentum_sink[i_momentum][0] = 1;
  momentum_sink[i_momentum][1] = 0;
  momentum_sink[i_momentum][2] = 1;

  //- momentum_sink[6] = (1,1,1)
  i_momentum = 6;
  momentum_sink[i_momentum][0] = 1;
  momentum_sink[i_momentum][1] = 1;
  momentum_sink[i_momentum][2] = 1;

  //- momentum_sink[7] = (2,0,0)
  i_momentum = 7;
  momentum_sink[i_momentum][0] = 2;
  momentum_sink[i_momentum][1] = 0;
  momentum_sink[i_momentum][2] = 0;

  //- momentum_sink[8] = (0,2,0)
  i_momentum = 8;
  momentum_sink[i_momentum][0] = 0;
  momentum_sink[i_momentum][1] = 2;
  momentum_sink[i_momentum][2] = 0;

  //- momentum_sink[9] = (0,0,2)
  i_momentum = 9;
  momentum_sink[i_momentum][0] = 0;
  momentum_sink[i_momentum][1] = 0;
  momentum_sink[i_momentum][2] = 2;


  std::ostream& log_file_previous = vout.getStream();
  std::ofstream log_file;

  if (m_filename_output != "stdout") {
    log_file.open(m_filename_output.c_str(), std::ios::app);
    vout.init(log_file);
  }

  std::vector<dcomplex> corr(Lt);

  GammaMatrix gm_src  = m_gmset->get_GM(m_gmset->GAMMA5);
  GammaMatrix gm_sink = m_gmset->get_GM(m_gmset->GAMMA5);
  for (int momentum = 0; momentum < N_momentum; momentum++) {
    vout.general(m_vl, "PS_momentum(%d %d %d) <-- PS correlator:\n",
                 momentum_sink[momentum][0],
                 momentum_sink[momentum][1],
                 momentum_sink[momentum][2]);
    meson_momentum_correlator(corr, momentum_sink[momentum], gm_sink, gm_src,
                              sq1, sq2, source_position);
    for (int t = 0; t < corr.size(); ++t) {
      vout.general(m_vl, "  %4d  %20.12e  %20.12e\n",
                   t, real(corr[t]), imag(corr[t]));
    }
  }

  gm_src  = m_gmset->get_GM(m_gmset->GAMMA1);
  gm_sink = m_gmset->get_GM(m_gmset->GAMMA1);
  for (int momentum = 0; momentum < N_momentum; momentum++) {
    vout.general(m_vl, "V1_momentum(%d %d %d) <-- V1 correlator:\n",
                 momentum_sink[momentum][0],
                 momentum_sink[momentum][1],
                 momentum_sink[momentum][2]);
    meson_momentum_correlator(corr, momentum_sink[momentum], gm_sink, gm_src,
                              sq1, sq2, source_position);
    for (int t = 0; t < corr.size(); ++t) {
      vout.general(m_vl, "  %4d  %20.12e  %20.12e\n",
                   t, real(corr[t]), imag(corr[t]));
    }
  }

  gm_src  = m_gmset->get_GM(m_gmset->GAMMA2);
  gm_sink = m_gmset->get_GM(m_gmset->GAMMA2);
  for (int momentum = 0; momentum < N_momentum; momentum++) {
    vout.general(m_vl, "V2_momentum(%d %d %d) <-- V2 correlator:\n",
                 momentum_sink[momentum][0],
                 momentum_sink[momentum][1],
                 momentum_sink[momentum][2]);
    meson_momentum_correlator(corr, momentum_sink[momentum], gm_sink, gm_src,
                              sq1, sq2, source_position);
    for (int t = 0; t < corr.size(); ++t) {
      vout.general(m_vl, "  %4d  %20.12e  %20.12e\n",
                   t, real(corr[t]), imag(corr[t]));
    }
  }

  gm_src  = m_gmset->get_GM(m_gmset->GAMMA3);
  gm_sink = m_gmset->get_GM(m_gmset->GAMMA3);
  for (int momentum = 0; momentum < N_momentum; momentum++) {
    vout.general(m_vl, "V3_momentum(%d %d %d) <-- V3 correlator:\n",
                 momentum_sink[momentum][0],
                 momentum_sink[momentum][1],
                 momentum_sink[momentum][2]);
    meson_momentum_correlator(corr, momentum_sink[momentum], gm_sink, gm_src,
                              sq1, sq2, source_position);
    for (int t = 0; t < corr.size(); ++t) {
      vout.general(m_vl, "  %4d  %20.12e  %20.12e\n",
                   t, real(corr[t]), imag(corr[t]));
    }
  }


  if (m_filename_output != "stdout") {
    log_file.close();
    vout.init(log_file_previous);
  }

  return EXIT_SUCCESS;
}


//====================================================================
void Corr2pt_4spinor::meson_momentum_correlator(std::vector<dcomplex>& corr_global,
                                                const std::vector<int>& momentum_sink,
                                                const GammaMatrix& gm_sink,
                                                const GammaMatrix& gm_src,
                                                const std::vector<Field_F>& sq1,
                                                const std::vector<Field_F>& sq2,
                                                const std::vector<int>& source_position)
{
  const int Nc = CommonParameters::Nc();
  const int Nd = CommonParameters::Nd();
  const int Lt = CommonParameters::Lt();
  const int Nt = CommonParameters::Nt();

  assert(corr_global.size() == Lt);

  const GammaMatrix gm5         = m_gmset->get_GM(m_gmset->GAMMA5);
  const GammaMatrix gm_gm5_src  = gm_src.mult(gm5);
  const GammaMatrix gm5_gm_sink = gm5.mult(gm_sink);

  std::vector<dcomplex> corr_local(Nt, cmplx(0.0, 0.0));
  for (int c0 = 0; c0 < Nc; ++c0) {
    for (int d0 = 0; d0 < Nd; ++d0) {
      int d1 = gm_gm5_src.index(d0);

      for (int t = 0; t < Nt; ++t) {
        dcomplex corr_t;

        contract_at_t(corr_t, momentum_sink, gm5_gm_sink, source_position,
                      sq1[c0 + Nc * d0], sq2[c0 + Nc * d1], t);

        corr_local[t] += gm_gm5_src.value(d0) * corr_t;
      }
    }
  }
  global_corr_t(corr_global, corr_local);
}


//====================================================================
double Corr2pt_4spinor::proton_test(const std::vector<Field_F>& sq_u,
                                    const std::vector<Field_F>& sq_d)
{
  const int Lt = CommonParameters::Lt();

  std::ostream& log_file_previous = vout.getStream();
  std::ofstream log_file;

  if (m_filename_output != "stdout") {
    log_file.open(m_filename_output.c_str(), std::ios::app);
    vout.init(log_file);
  }


  vout.general(m_vl, "proton <-- proton correlator(UNITY):\n");

  const GammaMatrix gm_unit   = m_gmset->get_GM(m_gmset->UNITY);
  const GammaMatrix gm_gamma0 = m_gmset->get_GM(m_gmset->GAMMA4);

  std::vector<dcomplex> p_corr_unity(Lt);
  proton_correlator(p_corr_unity, gm_unit, sq_u, sq_d);

  for (int t = 0; t < p_corr_unity.size(); t++) {
    vout.general(m_vl, "  %4d  %20.12e  %20.12e\n",
                 t, real(p_corr_unity[t]), imag(p_corr_unity[t]));
  }

  vout.general(m_vl, "proton <-- proton correlator(UPPER):\n");

  std::vector<dcomplex> p_corr_gamma0(Lt);
  proton_correlator(p_corr_gamma0, gm_gamma0, sq_u, sq_d);

  std::vector<dcomplex> p_corr_upper(Lt);
  for (int t = 0; t < p_corr_upper.size(); t++) {
    p_corr_upper[t] = (p_corr_unity[t] + p_corr_gamma0[t]) * 0.5;
    vout.general(m_vl, "  %4d  %20.12e  %20.12e\n",
                 t, real(p_corr_upper[t]), imag(p_corr_upper[t]));
  }

  vout.general(m_vl, "proton <-- proton correlator(GAMMA0):\n");

  for (int t = 0; t < p_corr_gamma0.size(); t++) {
    vout.general(m_vl, "  %4d  %20.12e  %20.12e\n",
                 t, real(p_corr_gamma0[t]), imag(p_corr_gamma0[t]));
  }


  if (m_filename_output != "stdout") {
    log_file.close();
    vout.init(log_file_previous);
  }

  const double result = real(p_corr_gamma0[0]);

  return result;
}


//====================================================================
void Corr2pt_4spinor::proton_correlator(std::vector<dcomplex>& corr_global,
                                        const GammaMatrix& gm,
                                        const std::vector<Field_F>& sq_u,
                                        const std::vector<Field_F>& sq_d)
{
  const int Nc = CommonParameters::Nc();
  const int Nd = CommonParameters::Nd();
  const int Lt = CommonParameters::Lt();
  const int Nt = CommonParameters::Nt();

  assert(Nc == 3);
  assert(corr_global.size() == Lt);

  const GammaMatrix gm5 = m_gmset->get_GM(m_gmset->GAMMA5);
  const GammaMatrix c   = m_gmset->get_GM(m_gmset->CHARGECONJG);
  const GammaMatrix cg5 = c.mult(gm5);

#ifdef DEBUG
  vout.general(m_vl, "i:\tgm5\t\t\t\tc\t\t\t\tcg5\t\t\t\tgm\n");
  for (int i = 0; i < Nd; i++) {
    vout.general(m_vl, "%d:\t %d %e %e \t %d  %e %e \t %d  %e %e \t %d  %e %e \n",
                 i,
                 gm5.index(i), real(gm5.value(i)), imag(gm5.value(i)),
                 c.index(i), real(c.value(i)), imag(c.value(i)),
                 cg5.index(i), real(cg5.value(i)), imag(cg5.value(i)),
                 gm.index(i), real(gm.value(i)), imag(gm.value(i))
                 );
  }
#endif

  const int FactNc = 6;
  // This is valid only when Nc =3, which was already asserted.

  std::vector<dcomplex> corr_local(Nt);

  for (int t = 0; t < Nt; t++) {
    vout.paranoiac(m_vl, "# t= %d\n", t);

    dcomplex sum = cmplx(0.0, 0.0);
    for (int i_alpha = 0; i_alpha < Nd; i_alpha++) {
      int i_alphaP  = gm.index(i_alpha);
      int i_alpha3  = i_alpha;
      int i_alpha3P = i_alphaP;

      for (int i_alpha1P = 0; i_alpha1P < Nd; i_alpha1P++) {
        int i_alpha2P = cg5.index(i_alpha1P);

        for (int ic123P = 0; ic123P < FactNc; ic123P++) {
          EpsilonTensor epsilon_tensor;

          int      ic1P   = epsilon_tensor.epsilon_3_index(ic123P, 0);
          int      ic2P   = epsilon_tensor.epsilon_3_index(ic123P, 1);
          int      ic3P   = epsilon_tensor.epsilon_3_index(ic123P, 2);
          dcomplex factor = gm.value(i_alpha)
                            * cg5.value(i_alpha1P)
                            * static_cast<double>(epsilon_tensor.epsilon_3_value(ic123P));

          dcomplex sum1;
          contract_at_t(sum1, cg5, i_alpha3,
                        sq_u[ic1P + Nc * i_alpha1P],
                        sq_d[ic2P + Nc * i_alpha2P],
                        sq_u[ic3P + Nc * i_alpha3P], t);

          dcomplex sum2;
          contract_at_t(sum2, cg5, i_alpha3,
                        sq_u[ic3P + Nc * i_alpha3P],
                        sq_d[ic2P + Nc * i_alpha2P],
                        sq_u[ic1P + Nc * i_alpha1P], t);

          sum += factor * (sum1 - sum2);
        }
      }
    }

    corr_local[t] = sum;
  } // t loop end.

  global_corr_t(corr_global, corr_local);
}


//====================================================================
// meson =tr(gm5.qn_sink.sq1.qn_src.gm5.(sq2)^\dagger);
void Corr2pt_4spinor::meson_correlator_x(std::vector<dcomplex>& meson,
                                         const GammaMatrix& qn_sink,
                                         const GammaMatrix& qn_src,
                                         const std::vector<Field_F>& sq1,
                                         const std::vector<Field_F>& sq2)
{
  int Nc = CommonParameters::Nc();
  int Nd = CommonParameters::Nd();
#if DEBUG
  int Lx = CommonParameters::Lx();
#endif
  int Nx = CommonParameters::Nx();

#if DEBUG
  assert(meson.size() == Lx);
#endif

  GammaMatrix gm_src, gm_sink, gm5;
  gm5     = m_gmset->get_GM(m_gmset->GAMMA5);
  gm_src  = qn_src.mult(gm5);
  gm_sink = gm5.mult(qn_sink);

  std::vector<dcomplex> corr_local(Nx);
  for (int i = 0; i < Nx; ++i) {
    corr_local[i] = 0.0;
  }

  for (int c0 = 0; c0 < Nc; ++c0) {
    for (int d0 = 0; d0 < Nd; ++d0) {
      int d1 = gm_src.index(d0);

      for (int x = 0; x < Nx; ++x) {
        dcomplex corr_x;
        contract_at_x(corr_x, gm_sink,
                      sq1[c0 + Nc * d0], sq2[c0 + Nc * d1], x);

        corr_local[x] += gm_src.value(d0) * corr_x;
      }
    }
  }

  global_corr_x(meson, corr_local);
}


//====================================================================
void Corr2pt_4spinor::meson_momentum_correlator_x(std::vector<dcomplex>& corr_global,
                                                  const std::vector<int>& momentum_sink,
                                                  const GammaMatrix& gm_sink,
                                                  const GammaMatrix& gm_src,
                                                  const std::vector<Field_F>& sq1,
                                                  const std::vector<Field_F>& sq2,
                                                  const std::vector<int>& source_position)
{
  int Nc = CommonParameters::Nc();
  int Nd = CommonParameters::Nd();
#if DEBUG
  int Lx = CommonParameters::Lx();
#endif
  int Nx = CommonParameters::Nx();
#if DEBUG
  assert(corr_global.size() == Lx);
#endif

  GammaMatrix gm_gm5_src, gm5_gm_sink, gm5;
  gm5         = m_gmset->get_GM(m_gmset->GAMMA5);
  gm_gm5_src  = gm_src.mult(gm5);
  gm5_gm_sink = gm5.mult(gm_sink);

  std::vector<dcomplex> corr_local(Nx);

  for (int c0 = 0; c0 < Nc; ++c0) {
    for (int d0 = 0; d0 < Nd; ++d0) {
      int d1 = gm_gm5_src.index(d0);

      for (int x = 0; x < Nx; ++x) {
        dcomplex corr_x;

        contract_at_x(corr_x, momentum_sink, gm5_gm_sink, source_position,
                      sq1[c0 + Nc * d0], sq2[c0 + Nc * d1], x);

        corr_local[x] += gm_gm5_src.value(d0) * corr_x;
      }
    }
  }

  global_corr_x(corr_global, corr_local);
}


//====================================================================
void Corr2pt_4spinor::proton_correlator_x(std::vector<dcomplex>& proton,
                                          const GammaMatrix& gm,
                                          const std::vector<Field_F>& squ,
                                          const std::vector<Field_F>& sqd)
{
  int Nc = CommonParameters::Nc();
  int Nd = CommonParameters::Nd();
#if DEBUG
  int Lx = CommonParameters::Lx();
#endif
  int Nx = CommonParameters::Nx();

  assert(Nc == 3);
#if DEBUG
  assert(proton.size() == Lx);
#endif

  EpsilonTensor epsilon_tensor;

  GammaMatrix cg5, c, gm5;
  gm5 = m_gmset->get_GM(m_gmset->GAMMA5);
  c   = m_gmset->get_GM(m_gmset->CHARGECONJG);
  cg5 = c.mult(gm5);

  int FactNc = 6;
  // This is valid only when Nc =3, which was already asserted.

  std::vector<dcomplex> corr_local(Nx);

  for (int ix = 0; ix < Nx; ix++) {
    dcomplex sum = 0.0;
    dcomplex sum1, sum2;

    for (int ialph = 0; ialph < Nd; ialph++) {
      int ialphP  = gm.index(ialph);
      int ialph3  = ialph;
      int ialph3P = ialphP;

      for (int ialph1P = 0; ialph1P < Nd; ialph1P++) {
        int ialph2P = cg5.index(ialph1P);

        for (int ic123P = 0; ic123P < FactNc; ic123P++) {
          int      ic1P   = epsilon_tensor.epsilon_3_index(ic123P, 0);
          int      ic2P   = epsilon_tensor.epsilon_3_index(ic123P, 1);
          int      ic3P   = epsilon_tensor.epsilon_3_index(ic123P, 2);
          dcomplex factor = gm.value(ialph)
                            * cg5.value(ialph1P)
                            * static_cast<double>(epsilon_tensor.epsilon_3_value(ic123P));

          contract_at_x(sum1, cg5, ialph3,
                        squ[ic1P + Nc * ialph1P],
                        sqd[ic2P + Nc * ialph2P],
                        squ[ic3P + Nc * ialph3P], ix);
          contract_at_x(sum2, cg5, ialph3,
                        squ[ic3P + Nc * ialph3P],
                        sqd[ic2P + Nc * ialph2P],
                        squ[ic1P + Nc * ialph1P], ix);
          sum += factor * (sum1 - sum2);
        }
      }
    }

    corr_local[ix] = sum;
  } // it loop end.

  global_corr_x(proton, corr_local);
}


/* moved by tanigchi
//====================================================================
void Corr2pt_4spinor::global_corr_t(std::vector<dcomplex>& corr_global,
                                    const std::vector<dcomplex>& corr_local)
{
  const int Lt = CommonParameters::Lt();
  const int Nt = CommonParameters::Nt();

  assert(corr_global.size() == Lt);
  assert(corr_local.size() == Nt);

  const int ipe_t = Communicator::ipe(3);

  std::vector<dcomplex> corr_tmp(Lt, cmplx(0.0, 0.0));

  for (int t = 0; t < Nt; ++t) {
    int t_global = t + ipe_t * Nt;
    corr_tmp[t_global] = corr_local[t];
  }

  for (int t_global = 0; t_global < Lt; ++t_global) {
    double cr_r = Communicator::reduce_sum(real(corr_tmp[t_global]));
    double cr_i = Communicator::reduce_sum(imag(corr_tmp[t_global]));

    corr_global[t_global] = cmplx(cr_r, cr_i);
  }
}
*/

//====================================================================
//============================================================END=====
