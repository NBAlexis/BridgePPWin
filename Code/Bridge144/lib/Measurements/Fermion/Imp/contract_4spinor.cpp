#include "BridgeLib_Private.h"
#if USE_IMP

/*!
        @file    $Id:: contract_4spinor.cpp #$

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 18:35:38 #$

        @version $LastChangedRevision: 1571 $
*/

#include "Measurements/Fermion/contract_4spinor.h"

#include <cassert>

#if defined USE_GROUP_SU3
#define NC      3
#define NC2     6
#define ND      4
#define NCD2    24
#define C1      0
#define C2      2
#define C3      4
#elif defined USE_GROUP_SU2
#define NC      2
#define NC2     4
#define ND      4
#define NCD2    16
#endif

//====================================================================
void contract_at_t(dcomplex& corr,
                   const GammaMatrix& gm_sink,
                   const Field_F& v1, const Field_F& v2,
                   const int time)
{
#if defined USE_GROUP_SU_N
  int NC   = CommonParameters::Nc();
  int ND   = CommonParameters::Nd();
  int NC2  = 2 * NC;
  int NCD2 = NC2 * ND;
#endif

  const int Nvol   = v1.nvol();
  const int Nvol_s = Nvol / CommonParameters::Nt();

  assert(Nvol == v2.nvol());
  assert(v1.nex() == 1);
  assert(v2.nex() == 1);

  const double *w1 = v1.ptr(0);
  const double *w2 = v2.ptr(0);

  int id1[ND];
  int id2[ND];
  for (int id = 0; id < ND; ++id) {
    id1[id] = id * NC2;
    id2[id] = gm_sink.index(id) * NC2;
  }

  double c_r[ND];
  double c_i[ND];
  for (int id = 0; id < ND; ++id) {
    c_r[id] = 0.0;
    c_i[id] = 0.0;
  }


  for (int ss = 0; ss < Nvol_s; ++ss) {
    int site = NCD2 * (ss + time * Nvol_s);

    for (int cc = 0; cc < NC; ++cc) {
      for (int id = 0; id < ND; ++id) {
        int ic1_r = 2 * cc + id1[id] + site;
        int ic2_r = 2 * cc + id2[id] + site;

        int ic1_i = 2 * cc + 1 + id1[id] + site;
        int ic2_i = 2 * cc + 1 + id2[id] + site;

        c_r[id] += w1[ic2_r] * w2[ic1_r]
                   + w1[ic2_i] * w2[ic1_i];

        c_i[id] += -w1[ic2_r] * w2[ic1_i]
                   + w1[ic2_i] * w2[ic1_r];
      }
    }
  }

  corr = cmplx(0.0, 0.0);
  for (int id = 0; id < ND; ++id) {
    corr += gm_sink.value(id) * cmplx(c_r[id], c_i[id]);
  }
}


//====================================================================
void contract_at_t(dcomplex& corr,
                   const std::vector<int>& momentum_sink,
                   const GammaMatrix& gm_sink,
                   const std::vector<int>& source_position,
                   const Field_F& v1, const Field_F& v2,
                   const int time)
{
#if defined USE_GROUP_SU_N
  int NC   = CommonParameters::Nc();
  int ND   = CommonParameters::Nd();
  int NC2  = 2 * NC;
  int NCD2 = NC2 * ND;
#endif

  const int Nvol   = v1.nvol();
  const int Nvol_s = Nvol / CommonParameters::Nt();
  const int Nx     = CommonParameters::Nx();
  const int Ny     = CommonParameters::Ny();
  const int Nz     = CommonParameters::Nz();
  const int Lx     = CommonParameters::Lx();
  const int Ly     = CommonParameters::Ly();
  const int Lz     = CommonParameters::Lz();

  assert(Nvol == v2.nvol());
  assert(v1.nex() == 1);
  assert(v2.nex() == 1);
  assert(momentum_sink.size() == ND - 1);

  const double *w1 = v1.ptr(0);
  const double *w2 = v2.ptr(0);

  int id1[ND];
  int id2[ND];
  for (int id = 0; id < ND; ++id) {
    id1[id] = id * NC2;
    id2[id] = gm_sink.index(id) * NC2;
  }

  double c_r[ND];
  double c_i[ND];
  for (int id = 0; id < ND; ++id) {
    c_r[id] = 0.0;
    c_i[id] = 0.0;
  }

  static const double PI = 4.0 * atan(1.0);
  std::vector<double> p_unit(ND - 1);
  p_unit[0] = (2.0 * PI / Lx) * momentum_sink[0];
  p_unit[1] = (2.0 * PI / Ly) * momentum_sink[1];
  p_unit[2] = (2.0 * PI / Lz) * momentum_sink[2];

  std::vector<int> ipe(ND - 1);
  ipe[0] = Communicator::ipe(0);
  ipe[1] = Communicator::ipe(1);
  ipe[2] = Communicator::ipe(2);


  for (int ss = 0; ss < Nvol_s; ++ss) {
    int site = NCD2 * (ss + time * Nvol_s);

    int x = ss % Nx;
    int y = ss % (Nx * Ny) / Nx;
    int z = ss % (Nx * Ny * Nz) / (Nx * Ny);

    int x_global = x + ipe[0] * Nx;
    int y_global = y + ipe[1] * Ny;
    int z_global = z + ipe[2] * Nz;

    double p_x = p_unit[0] * (x_global - source_position[0]);
    double p_y = p_unit[1] * (y_global - source_position[1]);
    double p_z = p_unit[2] * (z_global - source_position[2]);

    double cos_p_xyz = cos(p_x + p_y + p_z);
    double sin_p_xyz = sin(p_x + p_y + p_z);

    for (int cc = 0; cc < NC; ++cc) {
      for (int id = 0; id < ND; ++id) {
        int ic1_r = 2 * cc + id1[id] + site;
        int ic2_r = 2 * cc + id2[id] + site;

        int ic1_i = 2 * cc + 1 + id1[id] + site;
        int ic2_i = 2 * cc + 1 + id2[id] + site;

        double w1_w2_r = w1[ic2_r] * w2[ic1_r] + w1[ic2_i] * w2[ic1_i];
        double w1_w2_i = -w1[ic2_r] * w2[ic1_i] + w1[ic2_i] * w2[ic1_r];

        //- corr[s0] += v2^dagger * v1 * exp(i * p_i * x_i)
        c_r[id] += w1_w2_r * cos_p_xyz - w1_w2_i * sin_p_xyz;
        c_i[id] += w1_w2_r * sin_p_xyz + w1_w2_i * cos_p_xyz;
      }
    }
  }

  corr = cmplx(0.0, 0.0);
  for (int id = 0; id < ND; ++id) {
    corr += gm_sink.value(id) * cmplx(c_r[id], c_i[id]);
  }
}


//====================================================================
void contract_at_t(dcomplex& corr,
                   const GammaMatrix& gm_sink,
                   const int i_alpha,
                   const Field_F& v1, const Field_F& v2, const Field_F& v3,
                   const int time)
{
#if defined USE_GROUP_SU3
  const int Nvol   = v1.nvol();
  const int Nvol_s = Nvol / CommonParameters::Nt();

  assert(Nvol == v2.nvol());
  assert(Nvol == v3.nvol());
  assert(v1.nex() == 1);
  assert(v2.nex() == 1);
  assert(v3.nex() == 1);

  const double *w1 = v1.ptr(0);
  const double *w2 = v2.ptr(0);
  const double *w3 = v3.ptr(0);

  int id1[ND];
  int id2[ND];
  for (int id = 0; id < ND; ++id) {
    id1[id] = id * NC2;
    id2[id] = gm_sink.index(id) * NC2;
  }
  int id3 = i_alpha * NC2;

  double c_r[ND];
  double c_i[ND];
  for (int id = 0; id < ND; ++id) {
    c_r[id] = 0.0;
    c_i[id] = 0.0;
  }


  for (int ss = 0; ss < Nvol_s; ++ss) {
    int site = NCD2 * (ss + time * Nvol_s);

    for (int id = 0; id < ND; ++id) {
      int ic11_r = C1 + id1[id] + site;
      int ic22_r = C2 + id2[id] + site;
      int ic33_r = C3 + id3 + site;

      int ic11_i = C1 + 1 + id1[id] + site;
      int ic22_i = C2 + 1 + id2[id] + site;
      int ic33_i = C3 + 1 + id3 + site;

      int ic21_r = C2 + id1[id] + site;
      int ic32_r = C3 + id2[id] + site;
      int ic13_r = C1 + id3 + site;

      int ic21_i = C2 + 1 + id1[id] + site;
      int ic32_i = C3 + 1 + id2[id] + site;
      int ic13_i = C1 + 1 + id3 + site;

      int ic31_r = C3 + id1[id] + site;
      int ic12_r = C1 + id2[id] + site;
      int ic23_r = C2 + id3 + site;

      int ic31_i = C3 + 1 + id1[id] + site;
      int ic12_i = C1 + 1 + id2[id] + site;
      int ic23_i = C2 + 1 + id3 + site;


      c_r[id] += (w1[ic11_r] * w2[ic22_r] - w1[ic11_i] * w2[ic22_i]) * w3[ic33_r]
                 - (w1[ic11_r] * w2[ic22_i] + w1[ic11_i] * w2[ic22_r]) * w3[ic33_i];
      c_i[id] += (w1[ic11_r] * w2[ic22_r] - w1[ic11_i] * w2[ic22_i]) * w3[ic33_i]
                 + (w1[ic11_r] * w2[ic22_i] + w1[ic11_i] * w2[ic22_r]) * w3[ic33_r];

      c_r[id] += (w1[ic21_r] * w2[ic32_r] - w1[ic21_i] * w2[ic32_i]) * w3[ic13_r]
                 - (w1[ic21_r] * w2[ic32_i] + w1[ic21_i] * w2[ic32_r]) * w3[ic13_i];
      c_i[id] += (w1[ic21_r] * w2[ic32_r] - w1[ic21_i] * w2[ic32_i]) * w3[ic13_i]
                 + (w1[ic21_r] * w2[ic32_i] + w1[ic21_i] * w2[ic32_r]) * w3[ic13_r];

      c_r[id] += (w1[ic31_r] * w2[ic12_r] - w1[ic31_i] * w2[ic12_i]) * w3[ic23_r]
                 - (w1[ic31_r] * w2[ic12_i] + w1[ic31_i] * w2[ic12_r]) * w3[ic23_i];
      c_i[id] += (w1[ic31_r] * w2[ic12_r] - w1[ic31_i] * w2[ic12_i]) * w3[ic23_i]
                 + (w1[ic31_r] * w2[ic12_i] + w1[ic31_i] * w2[ic12_r]) * w3[ic23_r];

      c_r[id] -= (w1[ic31_r] * w2[ic22_r] - w1[ic31_i] * w2[ic22_i]) * w3[ic13_r]
                 - (w1[ic31_r] * w2[ic22_i] + w1[ic31_i] * w2[ic22_r]) * w3[ic13_i];
      c_i[id] -= (w1[ic31_r] * w2[ic22_r] - w1[ic31_i] * w2[ic22_i]) * w3[ic13_i]
                 + (w1[ic31_r] * w2[ic22_i] + w1[ic31_i] * w2[ic22_r]) * w3[ic13_r];

      c_r[id] -= (w1[ic21_r] * w2[ic12_r] - w1[ic21_i] * w2[ic12_i]) * w3[ic33_r]
                 - (w1[ic21_r] * w2[ic12_i] + w1[ic21_i] * w2[ic12_r]) * w3[ic33_i];
      c_i[id] -= (w1[ic21_r] * w2[ic12_r] - w1[ic21_i] * w2[ic12_i]) * w3[ic33_i]
                 + (w1[ic21_r] * w2[ic12_i] + w1[ic21_i] * w2[ic12_r]) * w3[ic33_r];

      c_r[id] -= (w1[ic11_r] * w2[ic32_r] - w1[ic11_i] * w2[ic32_i]) * w3[ic23_r]
                 - (w1[ic11_r] * w2[ic32_i] + w1[ic11_i] * w2[ic32_r]) * w3[ic23_i];
      c_i[id] -= (w1[ic11_r] * w2[ic32_r] - w1[ic11_i] * w2[ic32_i]) * w3[ic23_i]
                 + (w1[ic11_r] * w2[ic32_i] + w1[ic11_i] * w2[ic32_r]) * w3[ic23_r];
    }
  }

  corr = cmplx(0.0, 0.0);
  for (int id = 0; id < ND; ++id) {
    corr += gm_sink.value(id) * cmplx(c_r[id], c_i[id]);
  }
#endif  // (USE_GROUP_SU3)
}


//====================================================================
//============================================================END=====
#endif
