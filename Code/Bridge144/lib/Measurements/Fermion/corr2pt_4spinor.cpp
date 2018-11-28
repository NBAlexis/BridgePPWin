#include "BridgeLib_Private.h"

/*!
        @file    $Id:: corr2pt_4spinor.cpp #$

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-03-14 06:41:34 #$

        @version $LastChangedRevision: 1593 $
*/

#include "corr2pt_4spinor.h"

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

  int Nc = 3;
  int n;

  m_filename_output = "stdout";

  m_epsilon_index.resize(Nc * 6);

  n = 0;
  m_epsilon_index[Nc * n]     = 0;
  m_epsilon_index[1 + Nc * n] = 1;
  m_epsilon_index[2 + Nc * n] = 2;

  n = 1;
  m_epsilon_index[Nc * n]     = 1;
  m_epsilon_index[1 + Nc * n] = 2;
  m_epsilon_index[2 + Nc * n] = 0;

  n = 2;
  m_epsilon_index[Nc * n]     = 2;
  m_epsilon_index[1 + Nc * n] = 0;
  m_epsilon_index[2 + Nc * n] = 1;

  n = 3;
  m_epsilon_index[Nc * n]     = 2;
  m_epsilon_index[1 + Nc * n] = 1;
  m_epsilon_index[2 + Nc * n] = 0;

  n = 4;
  m_epsilon_index[Nc * n]     = 1;
  m_epsilon_index[1 + Nc * n] = 0;
  m_epsilon_index[2 + Nc * n] = 2;

  n = 5;
  m_epsilon_index[Nc * n]     = 0;
  m_epsilon_index[1 + Nc * n] = 2;
  m_epsilon_index[2 + Nc * n] = 1;
}


//====================================================================
double Corr2pt_4spinor::meson_all(const std::vector<Field_F>& sq1,
                                  const std::vector<Field_F>& sq2)
{
  int Lt = CommonParameters::Lt();

  std::vector<dcomplex> corr(Lt);
  GammaMatrix           gm_src, gm_sink;

  std::ostream& log_file_previous = vout.getStream();
  std::ofstream log_file;

  if (m_filename_output != "stdout") {
    log_file.open(m_filename_output.c_str(), std::ios::app);
    vout.init(log_file);
  }


  gm_src  = m_gmset->get_GM(m_gmset->GAMMA5);
  gm_sink = m_gmset->get_GM(m_gmset->GAMMA5);
  vout.general(m_vl, "PS <-- PS correlator:\n");
  meson_correlator(corr, gm_sink, gm_src, sq1, sq2);
  double result = real(corr[0]);

  gm_src  = m_gmset->get_GM(m_gmset->GAMMA1);
  gm_sink = m_gmset->get_GM(m_gmset->GAMMA1);
  vout.general(m_vl, "V1 <-- V1 correlator:\n");
  meson_correlator(corr, gm_sink, gm_src, sq1, sq2);

  gm_src  = m_gmset->get_GM(m_gmset->GAMMA2);
  gm_sink = m_gmset->get_GM(m_gmset->GAMMA2);
  vout.general(m_vl, "V2 <-- V2 correlator:\n");
  meson_correlator(corr, gm_sink, gm_src, sq1, sq2);

  gm_src  = m_gmset->get_GM(m_gmset->GAMMA3);
  gm_sink = m_gmset->get_GM(m_gmset->GAMMA3);
  vout.general(m_vl, "V3 <-- V3 correlator:\n");
  meson_correlator(corr, gm_sink, gm_src, sq1, sq2);

  gm_src  = m_gmset->get_GM(m_gmset->GAMMA5);
  gm_sink = m_gmset->get_GM(m_gmset->GAMMA45);
  vout.general(m_vl, "A4 <-- PS correlator:\n");
  meson_correlator(corr, gm_sink, gm_src, sq1, sq2);

  gm_src  = m_gmset->get_GM(m_gmset->GAMMA54);
  gm_sink = m_gmset->get_GM(m_gmset->GAMMA5);
  vout.general(m_vl, "PS <-- A4 correlator:\n");
  meson_correlator(corr, gm_sink, gm_src, sq1, sq2);

  gm_src  = m_gmset->get_GM(m_gmset->UNITY);
  gm_sink = m_gmset->get_GM(m_gmset->UNITY);
  vout.general(m_vl, "S <-- S correlator:\n");
  meson_correlator(corr, gm_sink, gm_src, sq1, sq2);

  gm_src  = m_gmset->get_GM(m_gmset->GAMMA51);
  gm_sink = m_gmset->get_GM(m_gmset->GAMMA51);
  vout.general(m_vl, "GAMMA51 <-- GAMMA51 correlator:\n");
  meson_correlator(corr, gm_sink, gm_src, sq1, sq2);

  gm_src  = m_gmset->get_GM(m_gmset->GAMMA52);
  gm_sink = m_gmset->get_GM(m_gmset->GAMMA52);
  vout.general(m_vl, "GAMMA52 <-- GAMMA52 correlator:\n");
  meson_correlator(corr, gm_sink, gm_src, sq1, sq2);

  gm_src  = m_gmset->get_GM(m_gmset->GAMMA53);
  gm_sink = m_gmset->get_GM(m_gmset->GAMMA53);
  vout.general(m_vl, "GAMMA53 <-- GAMMA53 correlator:\n");
  meson_correlator(corr, gm_sink, gm_src, sq1, sq2);

  gm_src  = m_gmset->get_GM(m_gmset->SIGMA12);
  gm_sink = m_gmset->get_GM(m_gmset->SIGMA12);
  vout.general(m_vl, "SIGMA12 <-- SIGMA12 correlator:\n");
  meson_correlator(corr, gm_sink, gm_src, sq1, sq2);

  gm_src  = m_gmset->get_GM(m_gmset->SIGMA23);
  gm_sink = m_gmset->get_GM(m_gmset->SIGMA23);
  vout.general(m_vl, "SIGMA23 <-- SIGMA23 correlator:\n");
  meson_correlator(corr, gm_sink, gm_src, sq1, sq2);

  gm_src  = m_gmset->get_GM(m_gmset->SIGMA31);
  gm_sink = m_gmset->get_GM(m_gmset->SIGMA31);
  vout.general(m_vl, "SIGMA31 <-- SIGMA31 correlator:\n");
  meson_correlator(corr, gm_sink, gm_src, sq1, sq2);


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
  int Nc = CommonParameters::Nc();
  int Nd = CommonParameters::Nd();
  int Lt = CommonParameters::Lt();
  int Nt = CommonParameters::Nt();

  assert(corr_global.size() == Lt);

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

        contract_at_t(corr_t, gm5_gm_sink,
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
}


//====================================================================
double Corr2pt_4spinor::meson_momentum_all(const std::vector<Field_F>& sq1,
                                           const std::vector<Field_F>& sq2,
                                           const std::vector<int>& source_position)
{
  const int Ndim = CommonParameters::Ndim();
  const int Lt   = CommonParameters::Lt();

  std::vector<dcomplex> corr(Lt);
  GammaMatrix           gm_src, gm_sink;

  const int N_momentum = 10;

  typedef std::vector<int>   MomentumSet;
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


  gm_src  = m_gmset->get_GM(m_gmset->GAMMA5);
  gm_sink = m_gmset->get_GM(m_gmset->GAMMA5);
  for (int i_momentum = 0; i_momentum < N_momentum; i_momentum++) {
    vout.general(m_vl, "PS_momentum(%d %d %d) <-- PS correlator:\n",
                 momentum_sink[i_momentum][0],
                 momentum_sink[i_momentum][1],
                 momentum_sink[i_momentum][2]);
    meson_momentum_correlator(corr, momentum_sink[i_momentum], gm_sink, gm_src,
                              sq1, sq2, source_position);
  }

  gm_src  = m_gmset->get_GM(m_gmset->GAMMA1);
  gm_sink = m_gmset->get_GM(m_gmset->GAMMA1);
  for (int i_momentum = 0; i_momentum < N_momentum; i_momentum++) {
    vout.general(m_vl, "V1_momentum(%d %d %d) <-- V1 correlator:\n",
                 momentum_sink[i_momentum][0],
                 momentum_sink[i_momentum][1],
                 momentum_sink[i_momentum][2]);
    meson_momentum_correlator(corr, momentum_sink[i_momentum], gm_sink, gm_src,
                              sq1, sq2, source_position);
  }

  gm_src  = m_gmset->get_GM(m_gmset->GAMMA2);
  gm_sink = m_gmset->get_GM(m_gmset->GAMMA2);
  for (int i_momentum = 0; i_momentum < N_momentum; i_momentum++) {
    vout.general(m_vl, "V2_momentum(%d %d %d) <-- V2 correlator:\n",
                 momentum_sink[i_momentum][0],
                 momentum_sink[i_momentum][1],
                 momentum_sink[i_momentum][2]);
    meson_momentum_correlator(corr, momentum_sink[i_momentum], gm_sink, gm_src,
                              sq1, sq2, source_position);
  }

  gm_src  = m_gmset->get_GM(m_gmset->GAMMA3);
  gm_sink = m_gmset->get_GM(m_gmset->GAMMA3);
  for (int i_momentum = 0; i_momentum < N_momentum; i_momentum++) {
    vout.general(m_vl, "V3_momentum(%d %d %d) <-- V3 correlator:\n",
                 momentum_sink[i_momentum][0],
                 momentum_sink[i_momentum][1],
                 momentum_sink[i_momentum][2]);
    meson_momentum_correlator(corr, momentum_sink[i_momentum], gm_sink, gm_src,
                              sq1, sq2, source_position);
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
  int Nc = CommonParameters::Nc();
  int Nd = CommonParameters::Nd();
  int Lt = CommonParameters::Lt();
  int Nt = CommonParameters::Nt();

  assert(corr_global.size() == Lt);

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
}


//====================================================================
double Corr2pt_4spinor::proton_test(const std::vector<Field_F>& sq_u,
                                    const std::vector<Field_F>& sq_d)
{
  int Lt = CommonParameters::Lt();

  std::vector<dcomplex> p_corr_unity(Lt), p_corr_gamma0(Lt), p_corr_upper(Lt);
  GammaMatrix           gm_unit, gm_gamma0;

  std::ostream& log_file_previous = vout.getStream();
  std::ofstream log_file;

  if (m_filename_output != "stdout") {
    log_file.open(m_filename_output.c_str(), std::ios::app);
    vout.init(log_file);
  }


  vout.general(m_vl, "proton <-- proton correlator(UNITY):\n");

  gm_unit   = m_gmset->get_GM(m_gmset->UNITY);
  gm_gamma0 = m_gmset->get_GM(m_gmset->GAMMA4);

  proton_correlator(p_corr_unity, gm_unit, sq_u, sq_d);

  for (int it = 0; it < p_corr_unity.size(); it++) {
    vout.general(m_vl, "  %4d  %20.12e  %20.12e\n",
                 it, real(p_corr_unity[it]), imag(p_corr_unity[it]));
  }

  vout.general(m_vl, "proton <-- proton correlator(UPPER):\n");

  proton_correlator(p_corr_gamma0, gm_gamma0, sq_u, sq_d);
  for (int it = 0; it < p_corr_upper.size(); it++) {
    p_corr_upper[it] = (p_corr_unity[it] + p_corr_gamma0[it]) * 0.5;
    vout.general(m_vl, "  %4d  %20.12e  %20.12e\n",
                 it, real(p_corr_upper[it]), imag(p_corr_upper[it]));
  }

  vout.general(m_vl, "proton <-- proton correlator(GAMMA0):\n");

  for (int it = 0; it < p_corr_gamma0.size(); it++) {
    vout.general(m_vl, "  %4d  %20.12e  %20.12e\n",
                 it, real(p_corr_gamma0[it]), imag(p_corr_gamma0[it]));
  }


  if (m_filename_output != "stdout") {
    log_file.close();
    vout.init(log_file_previous);
  }

  double result = real(p_corr_gamma0[0]);

  return result;
}


//====================================================================
void Corr2pt_4spinor::proton_correlator(std::vector<dcomplex>& corr_global,
                                        const GammaMatrix& gm,
                                        const std::vector<Field_F>& sq_u,
                                        const std::vector<Field_F>& sq_d)
{
  int Nc = CommonParameters::Nc();
  int Nd = CommonParameters::Nd();
  int Lt = CommonParameters::Lt();
  int Nt = CommonParameters::Nt();

  assert(Nc == 3);
  assert(corr_global.size() == Lt);

  GammaMatrix cg5, c, gm5;
  gm5 = m_gmset->get_GM(m_gmset->GAMMA5);
  c   = m_gmset->get_GM(m_gmset->CHARGECONJG);
  cg5 = c.mult(gm5);

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

  int FactNc = 6;
  // This is valid only when Nc =3, which was already asserted.

  std::vector<dcomplex> corr_local(Nt);

  for (int it = 0; it < Nt; it++) {
    vout.paranoiac(m_vl, "# it= %d\n", it);

    dcomplex sum = cmplx(0.0, 0.0);
    dcomplex sum1, sum2;

    for (int i_alpha = 0; i_alpha < Nd; i_alpha++) {
      int i_alphaP  = gm.index(i_alpha);
      int i_alpha3  = i_alpha;
      int i_alpha3P = i_alphaP;

      for (int i_alpha1P = 0; i_alpha1P < Nd; i_alpha1P++) {
        int i_alpha2P = cg5.index(i_alpha1P);

        for (int ic123P = 0; ic123P < FactNc; ic123P++) {
          int      ic1P   = epsilon_index(0, ic123P);
          int      ic2P   = epsilon_index(1, ic123P);
          int      ic3P   = epsilon_index(2, ic123P);
          dcomplex factor = gm.value(i_alpha)
                            * cg5.value(i_alpha1P) * epsilon_value(ic123P);

          contract_at_t(sum1, cg5, i_alpha3,
                        sq_u[ic1P + Nc * i_alpha1P],
                        sq_d[ic2P + Nc * i_alpha2P],
                        sq_u[ic3P + Nc * i_alpha3P], it);
          contract_at_t(sum2, cg5, i_alpha3,
                        sq_u[ic3P + Nc * i_alpha3P],
                        sq_d[ic2P + Nc * i_alpha2P],
                        sq_u[ic1P + Nc * i_alpha1P], it);
          sum += factor * (sum1 - sum2);
        }
      }
    }

    corr_local[it] = sum;
  } // it loop end.

  global_corr_t(corr_global, corr_local);
}


//====================================================================
void Corr2pt_4spinor::global_corr_t(std::vector<dcomplex>& corr_global,
                                    std::vector<dcomplex>& corr_local)
{
  int Lt = CommonParameters::Lt();
  int Nt = CommonParameters::Nt();

  assert(corr_global.size() == Lt);
  assert(corr_local.size() == Nt);

  std::vector<dcomplex> corr_tmp(Lt);

  int ipe_t = Communicator::ipe(3);

  for (int t = 0; t < Lt; ++t) {
    corr_tmp[t] = cmplx(0.0, 0.0);
  }

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
