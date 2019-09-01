/*!
        @file    corr4pt_4spinor.cpp

        @brief

        @author  Yusuke Namekawa (namekawa)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/
#include "BridgeLib_Private.h"
#include "corr4pt_4spinor.h"

const std::string Corr4pt_4spinor::class_name = "Corr4pt_4spinor";

//====================================================================
void Corr4pt_4spinor::set_parameters(const Parameters& params)
{
  m_filename_output = params.get_string("filename_output");
  if (m_filename_output.empty()) {
    m_filename_output = "stdout";
  }

  const string str_vlevel = params.get_string("verbose_level");
  m_vl = vout.set_verbose_level(str_vlevel);
}


//====================================================================
double Corr4pt_4spinor::meson_all(const std::vector<Field_F>& sq1,
                                  const std::vector<Field_F>& sq2,
                                  const std::vector<Field_F>& sq3,
                                  const std::vector<Field_F>& sq4)
{
  const int Lt = CommonParameters::Lt();

  std::vector<dcomplex> corr_global(Lt);
  GammaMatrix           gm_sink_12, gm_sink_34, gm_src_12, gm_src_34;

  std::ostream& log_file_previous = vout.getStream();
  std::ofstream log_file;

  if (m_filename_output != "stdout") {
    log_file.open(m_filename_output.c_str(), std::ios::app);
    vout.init(log_file);
  }


  gm_src_12  = m_gmset->get_GM(m_gmset->GAMMA5);
  gm_src_34  = m_gmset->get_GM(m_gmset->GAMMA5);
  gm_sink_12 = m_gmset->get_GM(m_gmset->GAMMA5);
  gm_sink_34 = m_gmset->get_GM(m_gmset->GAMMA5);
  vout.general(m_vl, "PS <-- PS 4pt_correlator:\n");
  meson_correlator(corr_global, gm_sink_12, gm_sink_34, gm_src_12, gm_src_34, sq1, sq2, sq3, sq4);
  double result = real(corr_global[0]);

#if 0
  gm_src_12  = m_gmset->get_GM(m_gmset->GAMMA1);
  gm_src_34  = m_gmset->get_GM(m_gmset->GAMMA1);
  gm_sink_12 = m_gmset->get_GM(m_gmset->GAMMA1);
  gm_sink_34 = m_gmset->get_GM(m_gmset->GAMMA1);
  vout.general(m_vl, "V1 <-- V1 4pt_correlator:\n");
  meson_correlator(corr_global, gm_sink_12, gm_sink_34, gm_src_12, gm_src_34, sq1, sq2, sq3, sq4);

  gm_src_12  = m_gmset->get_GM(m_gmset->GAMMA2);
  gm_src_34  = m_gmset->get_GM(m_gmset->GAMMA2);
  gm_sink_12 = m_gmset->get_GM(m_gmset->GAMMA2);
  gm_sink_34 = m_gmset->get_GM(m_gmset->GAMMA2);
  vout.general(m_vl, "V2 <-- V2 4pt_correlator:\n");
  meson_correlator(corr_global, gm_sink_12, gm_sink_34, gm_src_12, gm_src_34, sq1, sq2, sq3, sq4);

  gm_src_12  = m_gmset->get_GM(m_gmset->GAMMA3);
  gm_src_34  = m_gmset->get_GM(m_gmset->GAMMA3);
  gm_sink_12 = m_gmset->get_GM(m_gmset->GAMMA3);
  gm_sink_34 = m_gmset->get_GM(m_gmset->GAMMA3);
  vout.general(m_vl, "V3 <-- V3 4pt_correlator:\n");
  meson_correlator(corr_global, gm_sink_12, gm_sink_34, gm_src_12, gm_src_34, sq1, sq2, sq3, sq4);

  gm_src_12  = m_gmset->get_GM(m_gmset->GAMMA5);
  gm_src_34  = m_gmset->get_GM(m_gmset->GAMMA5);
  gm_sink_12 = m_gmset->get_GM(m_gmset->GAMMA45);
  gm_sink_34 = m_gmset->get_GM(m_gmset->GAMMA45);
  vout.general(m_vl, "A4 <-- PS 4pt_correlator:\n");
  meson_correlator(corr_global, gm_sink_12, gm_sink_34, gm_src_12, gm_src_34, sq1, sq2, sq3, sq4);

  gm_src_12  = m_gmset->get_GM(m_gmset->GAMMA54);
  gm_src_34  = m_gmset->get_GM(m_gmset->GAMMA54);
  gm_sink_12 = m_gmset->get_GM(m_gmset->GAMMA5);
  gm_sink_34 = m_gmset->get_GM(m_gmset->GAMMA5);
  vout.general(m_vl, "PS <-- A4 4pt_correlator:\n");
  meson_correlator(corr_global, gm_sink_12, gm_sink_34, gm_src_12, gm_src_34, sq1, sq2, sq3, sq4);

  gm_src_12  = m_gmset->get_GM(m_gmset->UNITY);
  gm_src_34  = m_gmset->get_GM(m_gmset->UNITY);
  gm_sink_12 = m_gmset->get_GM(m_gmset->UNITY);
  gm_sink_34 = m_gmset->get_GM(m_gmset->UNITY);
  vout.general(m_vl, "S <-- S 4pt_correlator:\n");
  meson_correlator(corr_global, gm_sink_12, gm_sink_34, gm_src_12, gm_src_34, sq1, sq2, sq3, sq4);

  gm_src_12  = m_gmset->get_GM(m_gmset->GAMMA51);
  gm_src_34  = m_gmset->get_GM(m_gmset->GAMMA51);
  gm_sink_12 = m_gmset->get_GM(m_gmset->GAMMA51);
  gm_sink_34 = m_gmset->get_GM(m_gmset->GAMMA51);
  vout.general(m_vl, "GAMMA51 <-- GAMMA51 4pt_correlator:\n");
  meson_correlator(corr_global, gm_sink_12, gm_sink_34, gm_src_12, gm_src_34, sq1, sq2, sq3, sq4);

  gm_src_12  = m_gmset->get_GM(m_gmset->GAMMA52);
  gm_src_34  = m_gmset->get_GM(m_gmset->GAMMA52);
  gm_sink_12 = m_gmset->get_GM(m_gmset->GAMMA52);
  gm_sink_34 = m_gmset->get_GM(m_gmset->GAMMA52);
  vout.general(m_vl, "GAMMA52 <-- GAMMA52 4pt_correlator:\n");
  meson_correlator(corr_global, gm_sink_12, gm_sink_34, gm_src_12, gm_src_34, sq1, sq2, sq3, sq4);

  gm_src_12  = m_gmset->get_GM(m_gmset->GAMMA53);
  gm_src_34  = m_gmset->get_GM(m_gmset->GAMMA53);
  gm_sink_12 = m_gmset->get_GM(m_gmset->GAMMA53);
  gm_sink_34 = m_gmset->get_GM(m_gmset->GAMMA53);
  vout.general(m_vl, "GAMMA53 <-- GAMMA53 4pt_correlator:\n");
  meson_correlator(corr_global, gm_sink_12, gm_sink_34, gm_src_12, gm_src_34, sq1, sq2, sq3, sq4);

  gm_src_12  = m_gmset->get_GM(m_gmset->SIGMA12);
  gm_src_34  = m_gmset->get_GM(m_gmset->SIGMA12);
  gm_sink_12 = m_gmset->get_GM(m_gmset->SIGMA12);
  gm_sink_34 = m_gmset->get_GM(m_gmset->SIGMA12);
  vout.general(m_vl, "SIGMA12 <-- SIGMA12 4pt_correlator:\n");
  meson_correlator(corr_global, gm_sink_12, gm_sink_34, gm_src_12, gm_src_34, sq1, sq2, sq3, sq4);

  gm_src_12  = m_gmset->get_GM(m_gmset->SIGMA23);
  gm_src_34  = m_gmset->get_GM(m_gmset->SIGMA23);
  gm_sink_12 = m_gmset->get_GM(m_gmset->SIGMA23);
  gm_sink_34 = m_gmset->get_GM(m_gmset->SIGMA23);
  vout.general(m_vl, "SIGMA23 <-- SIGMA23 4pt_correlator:\n");
  meson_correlator(corr_global, gm_sink_12, gm_sink_34, gm_src_12, gm_src_34, sq1, sq2, sq3, sq4);

  gm_src_12  = m_gmset->get_GM(m_gmset->SIGMA31);
  gm_src_34  = m_gmset->get_GM(m_gmset->SIGMA31);
  gm_sink_12 = m_gmset->get_GM(m_gmset->SIGMA31);
  gm_sink_34 = m_gmset->get_GM(m_gmset->SIGMA31);
  vout.general(m_vl, "SIGMA31 <-- SIGMA31 4pt_correlator:\n");
  meson_correlator(corr_global, gm_sink_12, gm_sink_34, gm_src_12, gm_src_34, sq1, sq2, sq3, sq4);
#endif

  if (m_filename_output != "stdout") {
    log_file.close();
    vout.init(log_file_previous);
  }

  return result;
}


//====================================================================
void Corr4pt_4spinor::meson_correlator(std::vector<dcomplex>& corr_global,
                                       const GammaMatrix& gm_sink_12,
                                       const GammaMatrix& gm_sink_34,
                                       const GammaMatrix& gm_src_12,
                                       const GammaMatrix& gm_src_34,
                                       const std::vector<Field_F>& sq1,
                                       const std::vector<Field_F>& sq2,
                                       const std::vector<Field_F>& sq3,
                                       const std::vector<Field_F>& sq4)
{
  const int Nc   = CommonParameters::Nc();
  const int Nd   = CommonParameters::Nd();
  const int Nx   = CommonParameters::Nx();
  const int Ny   = CommonParameters::Ny();
  const int Nz   = CommonParameters::Nz();
  const int Nt   = CommonParameters::Nt();
  const int Nvol = CommonParameters::Nvol();

  const int Lt = CommonParameters::Lt();

  assert(corr_global.size() == Lt);

  //- M_{ij}(x,t) := \bar{q_i}^c \Gamma_ij q_j^c
  //  Prop_ij(t)  := q_i(t_sink) \bar{q_j}(t_src)
  //
  //- corr_direct = tr(Prop_11(t) Prop_22^{\dagger}(t)) tr(Prop_33(t) Prop_44^{\dagger}(t))
  //               +tr(Prop_13(t) Prop_24^{\dagger}(t)) tr(Prop_31(t) Prop_42^{\dagger}(t))
  //             x   x'      x   x'
  //    t_sink  1 2 3 4     1 2 3 4 1 2
  //            | | | |    / / / / / /
  //            ^ v ^ v  +  v ^ v ^ v
  //            | | | |      / / / /
  //    t_src   1 2 3 4     1 2 3 4
  //            c0  c1      c0  c1
  //
  //- corr_cross = Prop_11(t)_{c0,c2} Prop_24^{\dagger}(t)_{c2,c1} Prop_33(t)_{c1,c3} Prop_42^{\dagger}(t)_{c3,c0}
  //              +Prop_31(t)_{c0,c3} Prop_44^{\dagger}(t)_{c3,c1} Prop_13(t)_{c1,c2} Prop_22^{\dagger}(t)_{c2,c0}
  //             = Prop_1124(t)_{c0,c1} Prop_3342(t)_{c1,c0}
  //              +Prop_3144(t)_{c0,c1} Prop_1322(t)_{c1,c0}
  //             x   x'      x   x'
  //    t_sink  1 2 3 4     1 2 3 4
  //            |  \|/       \|/  |
  //            ^   *   +     *   v
  //            |  /|\       /|\  |
  //    t_src   1 2 3 4     1 2 3 4
  //            c0  c1      c0  c1

  const GammaMatrix gm5 = m_gmset->get_GM(m_gmset->GAMMA5);

  const GammaMatrix gm_gm5_src_12 = gm_src_12.mult(gm5);
  const GammaMatrix gm_gm5_src_34 = gm_src_34.mult(gm5);

  const GammaMatrix gm5_gm_sink_12 = gm5.mult(gm_sink_12);
  const GammaMatrix gm5_gm_sink_34 = gm5.mult(gm_sink_34);

  //- corr_direct = tr(Prop_11(t) Prop_22^{\dagger}(t)) tr(Prop_33(t) Prop_44^{\dagger}(t))
  //               +tr(Prop_13(t) Prop_24^{\dagger}(t)) tr(Prop_31(t) Prop_42^{\dagger}(t))

  std::vector<dcomplex> corr_direct_1122_global(Lt);
  std::vector<dcomplex> corr_direct_3344_global(Lt);
  std::vector<dcomplex> corr_direct_1324_global(Lt);
  std::vector<dcomplex> corr_direct_3142_global(Lt);

  corr_direct(corr_direct_1122_global,
              gm5_gm_sink_12, gm_gm5_src_12,
              sq1, sq2);
  corr_direct(corr_direct_3344_global,
              gm5_gm_sink_34, gm_gm5_src_34,
              sq3, sq4);
  corr_direct(corr_direct_1324_global,
              gm5_gm_sink_12, gm_gm5_src_34,
              sq3, sq4);
  corr_direct(corr_direct_3142_global,
              gm5_gm_sink_34, gm_gm5_src_12,
              sq1, sq2);

  std::vector<dcomplex> corr_direct_global(Lt);
  for (int t_global = 0; t_global < Lt; ++t_global) {
    corr_direct_global[t_global]
      = corr_direct_1122_global[t_global] * corr_direct_3344_global[t_global]
        + corr_direct_1324_global[t_global] * corr_direct_3142_global[t_global];
  }

  //- corr_cross = Prop_11(t)_{c0,c2} Prop_24^{\dagger}(t)_{c2,c1} Prop_33(t)_{c1,c3} Prop_42^{\dagger}(t)_{c3,c0}
  //              +Prop_31(t)_{c0,c3} Prop_44^{\dagger}(t)_{c3,c1} Prop_13(t)_{c1,c2} Prop_22^{\dagger}(t)_{c2,c0}
  //             = Prop_1124(t)_{c0,c1} Prop_3342(t)_{c1,c0}
  //              +Prop_3144(t)_{c0,c1} Prop_1322(t)_{c1,c0}

  typedef std::vector<dcomplex> CorrSet;
  std::vector<CorrSet> corr_cross_1124_global(Lt);
  std::vector<CorrSet> corr_cross_3342_global(Lt);
  std::vector<CorrSet> corr_cross_1322_global(Lt);
  std::vector<CorrSet> corr_cross_3144_global(Lt);
  for (int t_global = 0; t_global < Lt; ++t_global) {
    corr_cross_1124_global[t_global].resize(Nc * Nd * Nc * Nd);
    corr_cross_3342_global[t_global].resize(Nc * Nd * Nc * Nd);
    corr_cross_1322_global[t_global].resize(Nc * Nd * Nc * Nd);
    corr_cross_3144_global[t_global].resize(Nc * Nd * Nc * Nd);
  }

  corr_cross_sink(corr_cross_1124_global,
                  gm5_gm_sink_12, gm_gm5_src_12,
                  sq1, sq4);
  corr_cross_sink(corr_cross_3342_global,
                  gm5_gm_sink_34, gm_gm5_src_34,
                  sq3, sq2);
  corr_cross_sink(corr_cross_3144_global,
                  gm5_gm_sink_34, gm_gm5_src_12,
                  sq1, sq4);
  corr_cross_sink(corr_cross_1322_global,
                  gm5_gm_sink_12, gm_gm5_src_34,
                  sq3, sq2);

  //- sum up at src
  std::vector<dcomplex> corr_cross_global(Lt, cmplx(0.0, 0.0));
  for (int t_global = 0; t_global < Lt; ++t_global) {
    for (int c0 = 0; c0 < Nc; ++c0) {
      for (int d0 = 0; d0 < Nd; ++d0) {
        for (int c1 = 0; c1 < Nc; ++c1) {
          for (int d1 = 0; d1 < Nd; ++d1) {
            int i_cd2_01 = c0 + Nc * d0 + Nc * Nd * c1 + Nc * Nd * Nc * d1;
            int i_cd2_10 = c1 + Nc * d1 + Nc * Nd * c0 + Nc * Nd * Nc * d0;

            dcomplex corr_1124 = corr_cross_1124_global[t_global][i_cd2_01];
            dcomplex corr_3342 = corr_cross_3342_global[t_global][i_cd2_10];

            corr_cross_global[t_global] += corr_1124 * corr_3342;

            dcomplex corr_3144 = corr_cross_3144_global[t_global][i_cd2_01];
            dcomplex corr_1322 = corr_cross_1322_global[t_global][i_cd2_10];

            corr_cross_global[t_global] += corr_3144 * corr_1322;
          }
        }
      }
    }
  }


  //- write outputs
  for (int t_global = 0; t_global < Lt; ++t_global) {
    corr_global[t_global] = corr_direct_global[t_global] - corr_cross_global[t_global];

    vout.general(m_vl, "  %4d  %20.12e  %20.12e\n",
                 t_global, real(corr_global[t_global]), imag(corr_global[t_global]));
  }
  vout.general(m_vl, "\n");


  vout.detailed(m_vl, "  4pt_correlator_direct:\n");
  for (int t_global = 0; t_global < Lt; ++t_global) {
    vout.detailed(m_vl, "  %4d  %20.12e  %20.12e\n",
                  t_global, real(corr_direct_global[t_global]), imag(corr_direct_global[t_global]));
  }
  vout.detailed(m_vl, "\n");

  vout.detailed(m_vl, "  4pt_correlator_cross:\n");
  for (int t_global = 0; t_global < Lt; ++t_global) {
    vout.detailed(m_vl, "  %4d  %20.12e  %20.12e\n",
                  t_global, real(corr_cross_global[t_global]), imag(corr_cross_global[t_global]));
  }
  vout.detailed(m_vl, "\n");
}


//====================================================================
double Corr4pt_4spinor::meson_momentum_all(const std::vector<Field_F>& sq1,
                                           const std::vector<Field_F>& sq2,
                                           const std::vector<Field_F>& sq3,
                                           const std::vector<Field_F>& sq4,
                                           const std::vector<int>& source_position)
{
  const int Ndim = CommonParameters::Ndim();
  const int Lt   = CommonParameters::Lt();

  std::vector<dcomplex> corr_global(Lt);
  GammaMatrix           gm_sink_12, gm_sink_34, gm_src_12, gm_src_34;

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


  gm_src_12  = m_gmset->get_GM(m_gmset->GAMMA5);
  gm_src_34  = m_gmset->get_GM(m_gmset->GAMMA5);
  gm_sink_12 = m_gmset->get_GM(m_gmset->GAMMA5);
  gm_sink_34 = m_gmset->get_GM(m_gmset->GAMMA5);
  for (int momentum = 0; momentum < N_momentum; momentum++) {
    vout.general(m_vl, "PS_momentum(%d %d %d) <-- PS correlator:\n",
                 momentum_sink[momentum][0],
                 momentum_sink[momentum][1],
                 momentum_sink[momentum][2]);
    meson_momentum_correlator(corr_global, momentum_sink[momentum],
                              gm_sink_12, gm_sink_34, gm_src_12, gm_src_34,
                              sq1, sq2, sq3, sq4, source_position);
  }

  gm_src_12  = m_gmset->get_GM(m_gmset->GAMMA1);
  gm_src_34  = m_gmset->get_GM(m_gmset->GAMMA1);
  gm_sink_12 = m_gmset->get_GM(m_gmset->GAMMA1);
  gm_sink_34 = m_gmset->get_GM(m_gmset->GAMMA1);
  for (int momentum = 0; momentum < N_momentum; momentum++) {
    vout.general(m_vl, "V1_momentum(%d %d %d) <-- V1 correlator:\n",
                 momentum_sink[momentum][0],
                 momentum_sink[momentum][1],
                 momentum_sink[momentum][2]);
    meson_momentum_correlator(corr_global, momentum_sink[momentum],
                              gm_sink_12, gm_sink_34, gm_src_12, gm_src_34,
                              sq1, sq2, sq3, sq4, source_position);
  }

  gm_src_12  = m_gmset->get_GM(m_gmset->GAMMA2);
  gm_src_34  = m_gmset->get_GM(m_gmset->GAMMA2);
  gm_sink_12 = m_gmset->get_GM(m_gmset->GAMMA2);
  gm_sink_34 = m_gmset->get_GM(m_gmset->GAMMA2);
  for (int momentum = 0; momentum < N_momentum; momentum++) {
    vout.general(m_vl, "V2_momentum(%d %d %d) <-- V2 correlator:\n",
                 momentum_sink[momentum][0],
                 momentum_sink[momentum][1],
                 momentum_sink[momentum][2]);
    meson_momentum_correlator(corr_global, momentum_sink[momentum],
                              gm_sink_12, gm_sink_34, gm_src_12, gm_src_34,
                              sq1, sq2, sq3, sq4, source_position);
  }

  gm_src_12  = m_gmset->get_GM(m_gmset->GAMMA3);
  gm_src_34  = m_gmset->get_GM(m_gmset->GAMMA3);
  gm_sink_12 = m_gmset->get_GM(m_gmset->GAMMA3);
  gm_sink_34 = m_gmset->get_GM(m_gmset->GAMMA3);
  for (int momentum = 0; momentum < N_momentum; momentum++) {
    vout.general(m_vl, "V3_momentum(%d %d %d) <-- V3 correlator:\n",
                 momentum_sink[momentum][0],
                 momentum_sink[momentum][1],
                 momentum_sink[momentum][2]);
    meson_momentum_correlator(corr_global, momentum_sink[momentum],
                              gm_sink_12, gm_sink_34, gm_src_12, gm_src_34,
                              sq1, sq2, sq3, sq4, source_position);
  }


  if (m_filename_output != "stdout") {
    log_file.close();
    vout.init(log_file_previous);
  }

  return EXIT_SUCCESS;
}


//====================================================================
void Corr4pt_4spinor::meson_momentum_correlator(std::vector<dcomplex>& corr_global,
                                                const std::vector<int>& momentum_sink,
                                                const GammaMatrix& gm_sink_12,
                                                const GammaMatrix& gm_sink_34,
                                                const GammaMatrix& gm_src_12,
                                                const GammaMatrix& gm_src_34,
                                                const std::vector<Field_F>& sq1,
                                                const std::vector<Field_F>& sq2,
                                                const std::vector<Field_F>& sq3,
                                                const std::vector<Field_F>& sq4,
                                                const std::vector<int>& source_position)
{
  const int Nc = CommonParameters::Nc();
  const int Nd = CommonParameters::Nd();
  const int Lt = CommonParameters::Lt();
  const int Nt = CommonParameters::Nt();

  assert(corr_global.size() == Lt);

#if 0
  GammaMatrix gm_gm5_src, gm5_gm_sink, gm5;
  gm5         = m_gmset->get_GM(m_gmset->GAMMA5);
  gm_gm5_src  = gm_src.mult(gm5);
  gm5_gm_sink = gm5.mult(gm_sink);

  std::vector<dcomplex> corr_local(Nt);

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

  for (int t = 0; t < corr_global.size(); ++t) {
    vout.general(m_vl, "  %4d  %20.12e  %20.12e\n",
                 t, real(corr_global[t]), imag(corr_global[t]));
  }
#endif
}


//====================================================================
void Corr4pt_4spinor::corr_direct(std::vector<dcomplex>& corr_direct_global,
                                  const GammaMatrix& gm5_gm_sink,
                                  const GammaMatrix& gm_gm5_src,
                                  const std::vector<Field_F>& sq1,
                                  const std::vector<Field_F>& sq2)
{
  const int Nc = CommonParameters::Nc();
  const int Nd = CommonParameters::Nd();
  const int Nt = CommonParameters::Nt();

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
  global_corr_t(corr_direct_global, corr_local);
}


//====================================================================
void Corr4pt_4spinor::corr_cross_sink(std::vector<std::vector<dcomplex> >& corr_cross_global,
                                      const GammaMatrix& gm5_gm_sink,
                                      const GammaMatrix& gm_gm5_src,
                                      const std::vector<Field_F>& sq1,
                                      const std::vector<Field_F>& sq2)
{
  const int Nc = CommonParameters::Nc();
  const int Nd = CommonParameters::Nd();
  const int Nt = CommonParameters::Nt();
  const int Lt = CommonParameters::Lt();

  typedef std::vector<dcomplex> CorrSet;
  std::vector<CorrSet> corr_local_idx(Nt);
  for (int t = 0; t < Nt; ++t) {
    corr_local_idx[t].resize(Nc * Nd * Nc * Nd);
  }
  for (int t = 0; t < Nt; ++t) {
    for (int i_cd2 = 0; i_cd2 < Nc * Nd * Nc * Nd; ++i_cd2) {
      corr_local_idx[t][i_cd2] = cmplx(0.0, 0.0);
    }
  }

  for (int c0 = 0; c0 < Nc; ++c0) {
    for (int d0 = 0; d0 < Nd; ++d0) {
      for (int c1 = 0; c1 < Nc; ++c1) {
        for (int d1 = 0; d1 < Nd; ++d1) {
          int d51 = gm_gm5_src.index(d1);

          for (int t = 0; t < Nt; ++t) {
            dcomplex corr_t;

            contract_at_t(corr_t, gm5_gm_sink,
                          sq1[c0 + Nc * d0], sq2[c1 + Nc * d51], t);

            int i_cd2 = c0 + Nc * d0 + Nc * Nd * c1 + Nc * Nd * Nc * d51;

            corr_local_idx[t][i_cd2] = gm_gm5_src.value(d1) * corr_t;
          }
        }
      }
    }
  }


  for (int i_cd2 = 0; i_cd2 < Nc * Nd * Nc * Nd; ++i_cd2) {
    std::vector<dcomplex> corr_local(Nt);
    std::vector<dcomplex> corr_global(Lt);

    for (int t = 0; t < Nt; ++t) {
      corr_local[t] = corr_local_idx[t][i_cd2];
    }

    global_corr_t(corr_global, corr_local);

    for (int t_global = 0; t_global < Lt; ++t_global) {
      corr_cross_global[t_global][i_cd2] = corr_global[t_global];
    }
  }
}


//====================================================================
void Corr4pt_4spinor::global_corr_t(std::vector<dcomplex>& corr_global,
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


//====================================================================
//============================================================END=====
