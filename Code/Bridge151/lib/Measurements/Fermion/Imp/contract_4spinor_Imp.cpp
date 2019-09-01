/*!
        @file    contract_4spinor.cpp

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/
#include "BridgeLib_Private.h"

#if USE_IMP
#include "Measurements/Fermion/contract_4spinor.h"

#include <cassert>

#define NDIM    4

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
// corr=(v2^*)_alpha (gm)_{alpha,beta} (v1)_beta
void contract_at_t(dcomplex& corr,
                   const GammaMatrix& gm_sink,
                   const Field_F& v1, const Field_F& v2,
                   const int time)
{
#if defined USE_GROUP_SU_N
  const int NC   = CommonParameters::Nc();
  const int ND   = CommonParameters::Nd();
  const int NC2  = 2 * NC;
  const int NCD2 = NC2 * ND;
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
void contract_at_t(std::vector<dcomplex>& corr_global,
                   const GammaMatrix& gm_sink,
                   const Field_F& v1, const Field_F& v2)
{
  const int Lt = CommonParameters::Lt();
  const int Nt = CommonParameters::Nt();

  assert(corr_global.size() == Lt);

  std::vector<dcomplex> corr_local(Nt, 0.0);
  for (int t = 0; t < Nt; ++t) {
    dcomplex corr_t;
    contract_at_t(corr_t, gm_sink, v1, v2, t);
    corr_local[t] += corr_t;
  }
  global_corr_t(corr_global, corr_local);
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
  const int NC   = CommonParameters::Nc();
  const int ND   = CommonParameters::Nd();
  const int NC2  = 2 * NC;
  const int NCD2 = NC2 * ND;
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
  assert(momentum_sink.size() == NDIM - 1);

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
void contract_at_t(std::vector<dcomplex>& corr_global,
                   const std::vector<int>& momentum_sink,
                   const GammaMatrix& gm_sink,
                   const std::vector<int>& source_position,
                   const Field_F& v1, const Field_F& v2)
{
  const int Lt = CommonParameters::Lt();
  const int Nt = CommonParameters::Nt();

  assert(corr_global.size() == Lt);

  std::vector<dcomplex> corr_local(Nt, 0.0);
  for (int t = 0; t < Nt; ++t) {
    dcomplex corr_t;
    contract_at_t(corr_t, momentum_sink, gm_sink, source_position, v1, v2, t);
    corr_local[t] += corr_t;
  }
  global_corr_t(corr_global, corr_local);
}


//====================================================================
void contract_at_t_cos(dcomplex& corr,
                       const std::vector<int>& momentum_sink,
                       const GammaMatrix& gm_sink,
                       const std::vector<int>& source_position,
                       const Field_F& v1, const Field_F& v2,
                       const int time)
{
#if defined USE_GROUP_SU_N
  const int NC   = CommonParameters::Nc();
  const int ND   = CommonParameters::Nd();
  const int NC2  = 2 * NC;
  const int NCD2 = NC2 * ND;
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
  assert(momentum_sink.size() == NDIM - 1);

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
    //double sin_p_xyz = sin(p_x + p_y + p_z);
    double sin_p_xyz = 0;

    for (int cc = 0; cc < NC; ++cc) {
      for (int id = 0; id < ND; ++id) {
        int ic1_r = 2 * cc + id1[id] + site;
        int ic2_r = 2 * cc + id2[id] + site;

        int ic1_i = 2 * cc + 1 + id1[id] + site;
        int ic2_i = 2 * cc + 1 + id2[id] + site;

        double w1_w2_r = w1[ic2_r] * w2[ic1_r] + w1[ic2_i] * w2[ic1_i];
        double w1_w2_i = -w1[ic2_r] * w2[ic1_i] + w1[ic2_i] * w2[ic1_r];
        //double w1_w2_i = 0;

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
void contract_at_t_cos(std::vector<dcomplex>& corr_global,
                       const std::vector<int>& momentum_sink,
                       const GammaMatrix& gm_sink,
                       const std::vector<int>& source_position,
                       const Field_F& v1, const Field_F& v2)
{
  const int Lt = CommonParameters::Lt();
  const int Nt = CommonParameters::Nt();

  assert(corr_global.size() == Lt);

  std::vector<dcomplex> corr_local(Nt, 0.0);
  for (int t = 0; t < Nt; ++t) {
    dcomplex corr_t;
    contract_at_t_cos(corr_t, momentum_sink, gm_sink, source_position, v1, v2, t);
    corr_local[t] += corr_t;
  }
  global_corr_t(corr_global, corr_local);
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
// corr=(v2^*)_alpha (gm)_{alpha,beta} (v1)_beta at x
void contract_at_x(dcomplex& corr, const GammaMatrix& gm_sink,
                   const Field_F& v1, const Field_F& v2,
                   int x)
{
#if defined USE_GROUP_SU_N
  const int NC   = CommonParameters::Nc();
  const int ND   = CommonParameters::Nd();
  const int NC2  = 2 * NC;
  const int NCD2 = NC2 * ND;
#endif

  const int Nvol   = v1.nvol();
  const int Nvol_s = Nvol / CommonParameters::Nt();
  const int Nx     = CommonParameters::Nx();
  const int Ny     = CommonParameters::Ny();
  const int Nz     = CommonParameters::Nz();
  const int Nt     = CommonParameters::Nt();

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

  for (int t = 0; t < Nt; ++t) {
    for (int z = 0; z < Nz; ++z) {
      for (int y = 0; y < Ny; ++y) {
        int site = NCD2 * (x + Nx * (y + Ny * (z + Nz * t)));

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
    }
  }

  corr = cmplx(0.0, 0.0);
  for (int id = 0; id < ND; ++id) {
    corr += gm_sink.value(id) * cmplx(c_r[id], c_i[id]);
  }
}


//====================================================================
void contract_at_x(std::vector<dcomplex>& corr_global,
                   const GammaMatrix& gm_sink,
                   const Field_F& v1, const Field_F& v2)
{
  const int Lx = CommonParameters::Lx();
  const int Nx = CommonParameters::Nx();

  assert(corr_global.size() == Lx);

  std::vector<dcomplex> corr_local(Nx, 0.0);
  for (int x = 0; x < Nx; ++x) {
    dcomplex corr_x;
    contract_at_x(corr_x, gm_sink, v1, v2, x);
    corr_local[x] += corr_x;
  }
  global_corr_x(corr_global, corr_local);
}


//====================================================================
void contract_at_x(dcomplex& corr,
                   const std::vector<int>& momentum_sink,
                   const GammaMatrix& gm_sink,
                   const std::vector<int>& source_position,
                   const Field_F& v1, const Field_F& v2,
                   const int x)
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
  const int Nt     = CommonParameters::Nt();
  const int Lx     = CommonParameters::Lx();
  const int Ly     = CommonParameters::Ly();
  const int Lz     = CommonParameters::Lz();
  const int Lt     = CommonParameters::Lt();

  assert(Nvol == v2.nvol());
  assert(v1.nex() == 1);
  assert(v2.nex() == 1);
  assert(momentum_sink.size() == NDIM - 1);
  assert(source_position.size() == NDIM);

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
  p_unit[0] = (2.0 * PI / Ly) * momentum_sink[0];
  p_unit[1] = (2.0 * PI / Lz) * momentum_sink[1];
  p_unit[2] = (2.0 * PI / Lt) * momentum_sink[2];

  std::vector<int> ipe(ND);
  ipe[0] = Communicator::ipe(0);
  ipe[1] = Communicator::ipe(1);
  ipe[2] = Communicator::ipe(2);
  ipe[3] = Communicator::ipe(3);

  for (int t = 0; t < Nt; ++t) {
    for (int z = 0; z < Nz; ++z) {
      for (int y = 0; y < Ny; ++y) {
        int site = NCD2 * (x + Nx * (y + Ny * (z + Nz * t)));

        int y_global = y + ipe[1] * Ny;
        int z_global = z + ipe[2] * Nz;
        int t_global = t + ipe[3] * Nt;

        double p_y = p_unit[0] * (y_global - source_position[1]);
        double p_z = p_unit[1] * (z_global - source_position[2]);
        double p_t = p_unit[2] * (t_global - source_position[3]);

        double cos_p_yzt = cos(p_t + p_y + p_z);
        double sin_p_yzt = sin(p_t + p_y + p_z);

        for (int cc = 0; cc < NC; ++cc) {
          for (int id = 0; id < ND; ++id) {
            int ic1_r = 2 * cc + id1[id] + site;
            int ic2_r = 2 * cc + id2[id] + site;

            int ic1_i = 2 * cc + 1 + id1[id] + site;
            int ic2_i = 2 * cc + 1 + id2[id] + site;

            double w1_w2_r = w1[ic2_r] * w2[ic1_r] + w1[ic2_i] * w2[ic1_i];
            double w1_w2_i = -w1[ic2_r] * w2[ic1_i] + w1[ic2_i] * w2[ic1_r];

            //- corr[s0] += v2^dagger * v1 * exp(i * p_i * x_i)
            c_r[id] += w1_w2_r * cos_p_yzt - w1_w2_i * sin_p_yzt;
            c_i[id] += w1_w2_r * sin_p_yzt + w1_w2_i * cos_p_yzt;
          }
        }
      }
    }
  }

  corr = cmplx(0.0, 0.0);
  for (int id = 0; id < ND; ++id) {
    corr += gm_sink.value(id) * cmplx(c_r[id], c_i[id]);
  }
}


//====================================================================
void contract_at_x(std::vector<dcomplex>& corr_global,
                   const std::vector<int>& momentum_sink,
                   const GammaMatrix& gm_sink,
                   const std::vector<int>& source_position,
                   const Field_F& v1, const Field_F& v2)
{
  const int Lx = CommonParameters::Lx();
  const int Nx = CommonParameters::Nx();

  assert(corr_global.size() == Lx);

  std::vector<dcomplex> corr_local(Nx, 0.0);
  for (int x = 0; x < Nx; ++x) {
    dcomplex corr_x;
    contract_at_x(corr_x, momentum_sink, gm_sink, source_position, v1, v2, x);
    corr_local[x] += corr_x;
  }
  global_corr_x(corr_global, corr_local);
}


//====================================================================
void contract_at_x_cos(dcomplex& corr,
                       const std::vector<int>& momentum_sink,
                       const GammaMatrix& gm_sink,
                       const std::vector<int>& source_position,
                       const Field_F& v1, const Field_F& v2,
                       const int x)
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
  const int Nt     = CommonParameters::Nt();
  const int Lx     = CommonParameters::Lx();
  const int Ly     = CommonParameters::Ly();
  const int Lz     = CommonParameters::Lz();
  const int Lt     = CommonParameters::Lt();

  assert(Nvol == v2.nvol());
  assert(v1.nex() == 1);
  assert(v2.nex() == 1);
  assert(momentum_sink.size() == NDIM - 1);
  assert(source_position.size() == NDIM);

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
  p_unit[0] = (2.0 * PI / Ly) * momentum_sink[0];
  p_unit[1] = (2.0 * PI / Lz) * momentum_sink[1];
  p_unit[2] = (2.0 * PI / Lt) * momentum_sink[2];

  std::vector<int> ipe(ND);
  ipe[0] = Communicator::ipe(0);
  ipe[1] = Communicator::ipe(1);
  ipe[2] = Communicator::ipe(2);
  ipe[3] = Communicator::ipe(3);


  for (int t = 0; t < Nt; ++t) {
    for (int z = 0; z < Nz; ++z) {
      for (int y = 0; y < Ny; ++y) {
        int site = NCD2 * (x + Nx * (y + Ny * (z + Nz * t)));

        int y_global = y + ipe[1] * Ny;
        int z_global = z + ipe[2] * Nz;
        int t_global = t + ipe[3] * Nt;

        double p_y = p_unit[0] * (y_global - source_position[1]);
        double p_z = p_unit[1] * (z_global - source_position[2]);
        double p_t = p_unit[2] * (t_global - source_position[3]);

        double cos_p_yzt = cos(p_t + p_y + p_z);
        //double sin_p_yzt = sin(p_t + p_y + p_z);
        double sin_p_yzt = 0;

        for (int cc = 0; cc < NC; ++cc) {
          for (int id = 0; id < ND; ++id) {
            int ic1_r = 2 * cc + id1[id] + site;
            int ic2_r = 2 * cc + id2[id] + site;

            int ic1_i = 2 * cc + 1 + id1[id] + site;
            int ic2_i = 2 * cc + 1 + id2[id] + site;

            double w1_w2_r = w1[ic2_r] * w2[ic1_r] + w1[ic2_i] * w2[ic1_i];
            double w1_w2_i = -w1[ic2_r] * w2[ic1_i] + w1[ic2_i] * w2[ic1_r];
            //double w1_w2_i = 0;

            //- corr[s0] += v2^dagger * v1 * exp(i * p_i * x_i)
            c_r[id] += w1_w2_r * cos_p_yzt - w1_w2_i * sin_p_yzt;
            c_i[id] += w1_w2_r * sin_p_yzt + w1_w2_i * cos_p_yzt;
          }
        }
      }
    }
  }

  corr = cmplx(0.0, 0.0);
  for (int id = 0; id < ND; ++id) {
    corr += gm_sink.value(id) * cmplx(c_r[id], c_i[id]);
  }
}


//====================================================================
void contract_at_x_cos(std::vector<dcomplex>& corr_global,
                       const std::vector<int>& momentum_sink,
                       const GammaMatrix& gm_sink,
                       const std::vector<int>& source_position,
                       const Field_F& v1, const Field_F& v2)
{
  const int Lx = CommonParameters::Lx();
  const int Nx = CommonParameters::Nx();

  assert(corr_global.size() == Lx);

  std::vector<dcomplex> corr_local(Nx, 0.0);
  for (int x = 0; x < Nx; ++x) {
    dcomplex corr_x;
    contract_at_x_cos(corr_x, momentum_sink, gm_sink, source_position, v1, v2, x);
    corr_local[x] += corr_x;
  }
  global_corr_x(corr_global, corr_local);
}


//====================================================================
void contract_at_x(dcomplex& corr,
                   const GammaMatrix& gm_sink,
                   const int i_alpha,
                   const Field_F& v1, const Field_F& v2, const Field_F& v3,
                   const int x)
{
#if defined USE_GROUP_SU3
  const int Nvol   = v1.nvol();
  const int Nvol_s = Nvol / CommonParameters::Nt();
  const int Nx     = CommonParameters::Nx();
  const int Ny     = CommonParameters::Ny();
  const int Nz     = CommonParameters::Nz();
  const int Nt     = CommonParameters::Nt();

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


  for (int t = 0; t < Nt; ++t) {
    for (int z = 0; z < Nz; ++z) {
      for (int y = 0; y < Ny; ++y) {
        int site = NCD2 * (x + Nx * (y + Ny * (z + Nz * t)));

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
    }
  }

  corr = cmplx(0.0, 0.0);
  for (int id = 0; id < ND; ++id) {
    corr += gm_sink.value(id) * cmplx(c_r[id], c_i[id]);
  }
#endif  // (USE_GROUP_SU3)
}


//====================================================================
// corr=(v2^*)_alpha (gm)_{alpha,beta} (v1)_beta
void contract_at_y(dcomplex& corr,
                   const GammaMatrix& gm_sink,
                   const Field_F& v1, const Field_F& v2,
                   const int y)
{
#if defined USE_GROUP_SU_N
  const int NC   = CommonParameters::Nc();
  const int ND   = CommonParameters::Nd();
  const int NC2  = 2 * NC;
  const int NCD2 = NC2 * ND;
#endif

  const int Nvol = v1.nvol();
  int       Nx   = CommonParameters::Nx();
  int       Ny   = CommonParameters::Ny();
  int       Nz   = CommonParameters::Nz();
  int       Nt   = CommonParameters::Nt();

  assert(Nvol == v1.nvol());
  assert(Nvol == v2.nvol());
  assert(v1.nex() == 1);
  assert(v2.nex() == 1);
  assert(y < Ny);

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

  for (int t = 0; t < Nt; ++t) {
    for (int z = 0; z < Nz; ++z) {
      for (int x = 0; x < Nx; ++x) {
        int site = NCD2 * (x + Nx * (y + Ny * (z + Nz * t)));
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
    }
  }

  corr = cmplx(0.0, 0.0);
  for (int id = 0; id < ND; ++id) {
    corr += gm_sink.value(id) * cmplx(c_r[id], c_i[id]);
  }
}


//====================================================================
void contract_at_y(std::vector<dcomplex>& corr_global,
                   const GammaMatrix& gm_sink,
                   const Field_F& v1, const Field_F& v2)
{
  const int Ly = CommonParameters::Ly();
  const int Ny = CommonParameters::Ny();

  assert(corr_global.size() == Ly);

  std::vector<dcomplex> corr_local(Ny, 0.0);
  for (int y = 0; y < Ny; ++y) {
    dcomplex corr_y;
    contract_at_y(corr_y, gm_sink, v1, v2, y);
    corr_local[y] += corr_y;
  }
  global_corr_y(corr_global, corr_local);
}


//====================================================================
void contract_at_y(dcomplex& corr,
                   const std::vector<int>& momentum_sink,
                   const GammaMatrix& gm_sink,
                   const std::vector<int>& source_position,
                   const Field_F& v1, const Field_F& v2,
                   const int y)
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
  const int Nt     = CommonParameters::Nt();
  const int Lx     = CommonParameters::Lx();
  const int Ly     = CommonParameters::Ly();
  const int Lz     = CommonParameters::Lz();
  const int Lt     = CommonParameters::Lt();

  assert(Nvol == v2.nvol());
  assert(v1.nex() == 1);
  assert(v2.nex() == 1);
  assert(momentum_sink.size() == NDIM - 1);
  assert(source_position.size() == NDIM);

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
  p_unit[0] = (2.0 * PI / Ly) * momentum_sink[0];
  p_unit[1] = (2.0 * PI / Lz) * momentum_sink[1];
  p_unit[2] = (2.0 * PI / Lt) * momentum_sink[2];

  std::vector<int> ipe(ND);
  ipe[0] = Communicator::ipe(0);
  ipe[1] = Communicator::ipe(1);
  ipe[2] = Communicator::ipe(2);
  ipe[3] = Communicator::ipe(3);

  for (int t = 0; t < Nt; ++t) {
    for (int z = 0; z < Nz; ++z) {
      for (int x = 0; x < Nx; ++x) {
        int site = NCD2 * (x + Nx * (y + Ny * (z + Nz * t)));

        int x_global = x + ipe[0] * Nx;
        int z_global = z + ipe[2] * Nz;
        int t_global = t + ipe[3] * Nt;

        double p_x = p_unit[0] * (x_global - source_position[0]);
        double p_z = p_unit[1] * (z_global - source_position[2]);
        double p_t = p_unit[2] * (t_global - source_position[3]);

        double cos_p_xzt = cos(p_t + p_x + p_z);
        double sin_p_xzt = sin(p_t + p_x + p_z);

        for (int cc = 0; cc < NC; ++cc) {
          for (int id = 0; id < ND; ++id) {
            int ic1_r = 2 * cc + id1[id] + site;
            int ic2_r = 2 * cc + id2[id] + site;

            int ic1_i = 2 * cc + 1 + id1[id] + site;
            int ic2_i = 2 * cc + 1 + id2[id] + site;

            double w1_w2_r = w1[ic2_r] * w2[ic1_r] + w1[ic2_i] * w2[ic1_i];
            double w1_w2_i = -w1[ic2_r] * w2[ic1_i] + w1[ic2_i] * w2[ic1_r];

            //- corr[s0] += v2^dagger * v1 * exp(i * p_i * x_i)
            c_r[id] += w1_w2_r * cos_p_xzt - w1_w2_i * sin_p_xzt;
            c_i[id] += w1_w2_r * sin_p_xzt + w1_w2_i * cos_p_xzt;
          }
        }
      }
    }
  }

  corr = cmplx(0.0, 0.0);
  for (int id = 0; id < ND; ++id) {
    corr += gm_sink.value(id) * cmplx(c_r[id], c_i[id]);
  }
}


//====================================================================
void contract_at_y(std::vector<dcomplex>& corr_global,
                   const std::vector<int>& momentum_sink,
                   const GammaMatrix& gm_sink,
                   const std::vector<int>& source_position,
                   const Field_F& v1, const Field_F& v2)
{
  const int Ly = CommonParameters::Ly();
  const int Ny = CommonParameters::Ny();

  assert(corr_global.size() == Ly);

  std::vector<dcomplex> corr_local(Ny, 0.0);
  for (int y = 0; y < Ny; ++y) {
    dcomplex corr_y;
    contract_at_y(corr_y, momentum_sink, gm_sink, source_position, v1, v2, y);
    corr_local[y] += corr_y;
  }
  global_corr_y(corr_global, corr_local);
}


//====================================================================
void contract_at_y_cos(dcomplex& corr,
                       const std::vector<int>& momentum_sink,
                       const GammaMatrix& gm_sink,
                       const std::vector<int>& source_position,
                       const Field_F& v1, const Field_F& v2,
                       const int y)
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
  const int Nt     = CommonParameters::Nt();
  const int Lx     = CommonParameters::Lx();
  const int Ly     = CommonParameters::Ly();
  const int Lz     = CommonParameters::Lz();
  const int Lt     = CommonParameters::Lt();

  assert(Nvol == v2.nvol());
  assert(v1.nex() == 1);
  assert(v2.nex() == 1);
  assert(momentum_sink.size() == NDIM - 1);
  assert(source_position.size() == NDIM);

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
  p_unit[0] = (2.0 * PI / Ly) * momentum_sink[0];
  p_unit[1] = (2.0 * PI / Lz) * momentum_sink[1];
  p_unit[2] = (2.0 * PI / Lt) * momentum_sink[2];

  std::vector<int> ipe(ND);
  ipe[0] = Communicator::ipe(0);
  ipe[1] = Communicator::ipe(1);
  ipe[2] = Communicator::ipe(2);
  ipe[3] = Communicator::ipe(3);


  for (int t = 0; t < Nt; ++t) {
    for (int z = 0; z < Nz; ++z) {
      for (int x = 0; x < Nx; ++x) {
        int site = NCD2 * (x + Nx * (y + Ny * (z + Nz * t)));

        int x_global = x + ipe[0] * Nx;
        int z_global = z + ipe[2] * Nz;
        int t_global = t + ipe[3] * Nt;

        double p_x = p_unit[0] * (x_global - source_position[0]);
        double p_z = p_unit[1] * (z_global - source_position[2]);
        double p_t = p_unit[2] * (t_global - source_position[3]);

        double cos_p_xzt = cos(p_t + p_x + p_z);
        //double sin_p_xzt = sin(p_t + p_x + p_z);
        double sin_p_xzt = 0;

        for (int cc = 0; cc < NC; ++cc) {
          for (int id = 0; id < ND; ++id) {
            int ic1_r = 2 * cc + id1[id] + site;
            int ic2_r = 2 * cc + id2[id] + site;

            int ic1_i = 2 * cc + 1 + id1[id] + site;
            int ic2_i = 2 * cc + 1 + id2[id] + site;

            double w1_w2_r = w1[ic2_r] * w2[ic1_r] + w1[ic2_i] * w2[ic1_i];
            double w1_w2_i = -w1[ic2_r] * w2[ic1_i] + w1[ic2_i] * w2[ic1_r];
            //double w1_w2_i = 0;

            //- corr[s0] += v2^dagger * v1 * exp(i * p_i * x_i)
            c_r[id] += w1_w2_r * cos_p_xzt - w1_w2_i * sin_p_xzt;
            c_i[id] += w1_w2_r * sin_p_xzt + w1_w2_i * cos_p_xzt;
          }
        }
      }
    }
  }

  corr = cmplx(0.0, 0.0);
  for (int id = 0; id < ND; ++id) {
    corr += gm_sink.value(id) * cmplx(c_r[id], c_i[id]);
  }
}


//====================================================================
void contract_at_y_cos(std::vector<dcomplex>& corr_global,
                       const std::vector<int>& momentum_sink,
                       const GammaMatrix& gm_sink,
                       const std::vector<int>& source_position,
                       const Field_F& v1, const Field_F& v2)
{
  const int Ly = CommonParameters::Ly();
  const int Ny = CommonParameters::Ny();

  assert(corr_global.size() == Ly);

  std::vector<dcomplex> corr_local(Ny, 0.0);
  for (int y = 0; y < Ny; ++y) {
    dcomplex corr_y;
    contract_at_y_cos(corr_y, momentum_sink, gm_sink, source_position, v1, v2, y);
    corr_local[y] += corr_y;
  }
  global_corr_y(corr_global, corr_local);
}


//====================================================================
// corr=(v2^*)_alpha (gm)_{alpha,beta} (v1)_beta
void contract_at_z(dcomplex& corr,
                   const GammaMatrix& gm_sink,
                   const Field_F& v1, const Field_F& v2,
                   const int z)
{
#if defined USE_GROUP_SU_N
  const int NC   = CommonParameters::Nc();
  const int ND   = CommonParameters::Nd();
  const int NC2  = 2 * NC;
  const int NCD2 = NC2 * ND;
#endif

  int       Nx   = CommonParameters::Nx();
  int       Ny   = CommonParameters::Ny();
  int       Nz   = CommonParameters::Nz();
  int       Nt   = CommonParameters::Nt();
  const int Nvol = v1.nvol();

  assert(Nvol == v1.nvol());
  assert(Nvol == v2.nvol());
  assert(v1.nex() == 1);
  assert(v2.nex() == 1);
  assert(z < Nz);

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

  for (int t = 0; t < Nt; ++t) {
    for (int y = 0; y < Ny; ++y) {
      for (int x = 0; x < Nx; ++x) {
        int site = NCD2 * (x + Nx * (y + Ny * (z + Nz * t)));
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
    }
  }

  corr = cmplx(0.0, 0.0);
  for (int id = 0; id < ND; ++id) {
    corr += gm_sink.value(id) * cmplx(c_r[id], c_i[id]);
  }
}


//====================================================================
void contract_at_z(std::vector<dcomplex>& corr_global,
                   const GammaMatrix& gm_sink,
                   const Field_F& v1, const Field_F& v2)
{
  const int Lz = CommonParameters::Lz();
  const int Nz = CommonParameters::Nz();

  assert(corr_global.size() == Lz);

  std::vector<dcomplex> corr_local(Nz, 0.0);
  for (int z = 0; z < Nz; ++z) {
    dcomplex corr_z;
    contract_at_z(corr_z, gm_sink, v1, v2, z);
    corr_local[z] += corr_z;
  }
  global_corr_z(corr_global, corr_local);
}


//====================================================================
void contract_at_z(dcomplex& corr,
                   const std::vector<int>& momentum_sink,
                   const GammaMatrix& gm_sink,
                   const std::vector<int>& source_position,
                   const Field_F& v1, const Field_F& v2,
                   const int z)
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
  const int Nt     = CommonParameters::Nt();
  const int Lx     = CommonParameters::Lx();
  const int Ly     = CommonParameters::Ly();
  const int Lz     = CommonParameters::Lz();
  const int Lt     = CommonParameters::Lt();

  assert(Nvol == v2.nvol());
  assert(v1.nex() == 1);
  assert(v2.nex() == 1);
  assert(momentum_sink.size() == NDIM - 1);
  assert(source_position.size() == NDIM);

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
  p_unit[0] = (2.0 * PI / Ly) * momentum_sink[0];
  p_unit[1] = (2.0 * PI / Lz) * momentum_sink[1];
  p_unit[2] = (2.0 * PI / Lt) * momentum_sink[2];

  std::vector<int> ipe(ND);
  ipe[0] = Communicator::ipe(0);
  ipe[1] = Communicator::ipe(1);
  ipe[2] = Communicator::ipe(2);
  ipe[3] = Communicator::ipe(3);

  for (int t = 0; t < Nt; ++t) {
    for (int y = 0; y < Ny; ++y) {
      for (int x = 0; x < Nx; ++x) {
        int site = NCD2 * (x + Nx * (y + Ny * (z + Nz * t)));

        int x_global = x + ipe[0] * Nx;
        int y_global = y + ipe[1] * Ny;
        int t_global = t + ipe[3] * Nt;

        double p_x = p_unit[0] * (x_global - source_position[0]);
        double p_y = p_unit[1] * (y_global - source_position[1]);
        double p_t = p_unit[2] * (t_global - source_position[3]);

        double cos_p_xyt = cos(p_t + p_x + p_y);
        double sin_p_xyt = sin(p_t + p_x + p_y);

        for (int cc = 0; cc < NC; ++cc) {
          for (int id = 0; id < ND; ++id) {
            int ic1_r = 2 * cc + id1[id] + site;
            int ic2_r = 2 * cc + id2[id] + site;

            int ic1_i = 2 * cc + 1 + id1[id] + site;
            int ic2_i = 2 * cc + 1 + id2[id] + site;

            double w1_w2_r = w1[ic2_r] * w2[ic1_r] + w1[ic2_i] * w2[ic1_i];
            double w1_w2_i = -w1[ic2_r] * w2[ic1_i] + w1[ic2_i] * w2[ic1_r];

            //- corr[s0] += v2^dagger * v1 * exp(i * p_i * x_i)
            c_r[id] += w1_w2_r * cos_p_xyt - w1_w2_i * sin_p_xyt;
            c_i[id] += w1_w2_r * sin_p_xyt + w1_w2_i * cos_p_xyt;
          }
        }
      }
    }
  }

  corr = cmplx(0.0, 0.0);
  for (int id = 0; id < ND; ++id) {
    corr += gm_sink.value(id) * cmplx(c_r[id], c_i[id]);
  }
}


//====================================================================
void contract_at_z(std::vector<dcomplex>& corr_global,
                   const std::vector<int>& momentum_sink,
                   const GammaMatrix& gm_sink,
                   const std::vector<int>& source_position,
                   const Field_F& v1, const Field_F& v2)
{
  const int Lz = CommonParameters::Lz();
  const int Nz = CommonParameters::Nz();

  assert(corr_global.size() == Lz);

  std::vector<dcomplex> corr_local(Nz, 0.0);
  for (int z = 0; z < Nz; ++z) {
    dcomplex corr_z;
    contract_at_z(corr_z, momentum_sink, gm_sink, source_position, v1, v2, z);
    corr_local[z] += corr_z;
  }
  global_corr_z(corr_global, corr_local);
}


//====================================================================
void contract_at_z_cos(dcomplex& corr,
                       const std::vector<int>& momentum_sink,
                       const GammaMatrix& gm_sink,
                       const std::vector<int>& source_position,
                       const Field_F& v1, const Field_F& v2,
                       const int z)
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
  const int Nt     = CommonParameters::Nt();
  const int Lx     = CommonParameters::Lx();
  const int Ly     = CommonParameters::Ly();
  const int Lz     = CommonParameters::Lz();
  const int Lt     = CommonParameters::Lt();

  assert(Nvol == v2.nvol());
  assert(v1.nex() == 1);
  assert(v2.nex() == 1);
  assert(momentum_sink.size() == NDIM - 1);
  assert(source_position.size() == NDIM);

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
  p_unit[0] = (2.0 * PI / Ly) * momentum_sink[0];
  p_unit[1] = (2.0 * PI / Lz) * momentum_sink[1];
  p_unit[2] = (2.0 * PI / Lt) * momentum_sink[2];

  std::vector<int> ipe(ND);
  ipe[0] = Communicator::ipe(0);
  ipe[1] = Communicator::ipe(1);
  ipe[2] = Communicator::ipe(2);
  ipe[3] = Communicator::ipe(3);


  for (int t = 0; t < Nt; ++t) {
    for (int y = 0; y < Ny; ++y) {
      for (int x = 0; x < Nx; ++x) {
        int site = NCD2 * (x + Nx * (y + Ny * (z + Nz * t)));

        int x_global = x + ipe[0] * Nx;
        int y_global = y + ipe[1] * Ny;
        int t_global = t + ipe[3] * Nt;

        double p_x = p_unit[0] * (x_global - source_position[0]);
        double p_y = p_unit[1] * (y_global - source_position[1]);
        double p_t = p_unit[2] * (t_global - source_position[2]);

        double cos_p_xyt = cos(p_t + p_x + p_y);
        //double sin_p_xyt = sin(p_t + p_x + p_y);
        double sin_p_xyt = 0;

        for (int cc = 0; cc < NC; ++cc) {
          for (int id = 0; id < ND; ++id) {
            int ic1_r = 2 * cc + id1[id] + site;
            int ic2_r = 2 * cc + id2[id] + site;

            int ic1_i = 2 * cc + 1 + id1[id] + site;
            int ic2_i = 2 * cc + 1 + id2[id] + site;

            double w1_w2_r = w1[ic2_r] * w2[ic1_r] + w1[ic2_i] * w2[ic1_i];
            double w1_w2_i = -w1[ic2_r] * w2[ic1_i] + w1[ic2_i] * w2[ic1_r];
            //double w1_w2_i = 0;

            //- corr[s0] += v2^dagger * v1 * exp(i * p_i * x_i)
            c_r[id] += w1_w2_r * cos_p_xyt - w1_w2_i * sin_p_xyt;
            c_i[id] += w1_w2_r * sin_p_xyt + w1_w2_i * cos_p_xyt;
          }
        }
      }
    }
  }

  corr = cmplx(0.0, 0.0);
  for (int id = 0; id < ND; ++id) {
    corr += gm_sink.value(id) * cmplx(c_r[id], c_i[id]);
  }
}


//====================================================================
void contract_at_z_cos(std::vector<dcomplex>& corr_global,
                       const std::vector<int>& momentum_sink,
                       const GammaMatrix& gm_sink,
                       const std::vector<int>& source_position,
                       const Field_F& v1, const Field_F& v2)
{
  const int Lz = CommonParameters::Lz();
  const int Nz = CommonParameters::Nz();

  assert(corr_global.size() == Lz);

  std::vector<dcomplex> corr_local(Nz, 0.0);
  for (int z = 0; z < Nz; ++z) {
    dcomplex corr_z;
    contract_at_z_cos(corr_z, momentum_sink, gm_sink, source_position, v1, v2, z);
    corr_local[z] += corr_z;
  }
  global_corr_z(corr_global, corr_local);
}


//====================================================================
void global_corr_x(std::vector<dcomplex>& corr_global,
                   std::vector<dcomplex>& corr_local)
{
  int Lx = CommonParameters::Lx();
  int Nx = CommonParameters::Nx();

  assert(corr_global.size() == Lx);
  assert(corr_local.size() == Nx);

  std::vector<dcomplex> corr_tmp(Lx, 0);

  int ipex = Communicator::ipe(0);

  for (int x = 0; x < Nx; ++x) {
    int x2 = x + ipex * Nx;
    corr_tmp[x2] = corr_local[x];
  }

  for (int x = 0; x < Lx; ++x) {
    double crr = Communicator::reduce_sum(real(corr_tmp[x]));
    double cri = Communicator::reduce_sum(imag(corr_tmp[x]));
    corr_global[x] = cmplx(crr, cri);
  }
}


//====================================================================
void global_corr_y(std::vector<dcomplex>& corr_global,
                   std::vector<dcomplex>& corr_local)
{
  int Ly = CommonParameters::Ly();
  int Ny = CommonParameters::Ny();

  assert(corr_global.size() == Ly);
  assert(corr_local.size() == Ny);

  std::vector<dcomplex> corr_tmp(Ly, 0);

  int ipey = Communicator::ipe(1);

  for (int y = 0; y < Ny; ++y) {
    int y2 = y + ipey * Ny;
    corr_tmp[y2] = corr_local[y];
  }

  for (int y = 0; y < Ly; ++y) {
    double crr = Communicator::reduce_sum(real(corr_tmp[y]));
    double cri = Communicator::reduce_sum(imag(corr_tmp[y]));
    corr_global[y] = cmplx(crr, cri);
  }
}


//====================================================================
void global_corr_z(std::vector<dcomplex>& corr_global,
                   std::vector<dcomplex>& corr_local)
{
  int Lz = CommonParameters::Lz();
  int Nz = CommonParameters::Nz();

  assert(corr_global.size() == Lz);
  assert(corr_local.size() == Nz);

  std::vector<dcomplex> corr_tmp(Lz);

  int ipez = Communicator::ipe(2);

  for (int z = 0; z < Nz; ++z) {
    int z2 = z + ipez * Nz;
    corr_tmp[z2] = corr_local[z];
  }

  for (int z = 0; z < Lz; ++z) {
    double crr = Communicator::reduce_sum(real(corr_tmp[z]));
    double cri = Communicator::reduce_sum(imag(corr_tmp[z]));
    corr_global[z] = cmplx(crr, cri);
  }
}


//====================================================================
void global_corr_t(std::vector<dcomplex>& corr_global,
                   std::vector<dcomplex>& corr_local)
{
  int Lt = CommonParameters::Lt();
  int Nt = CommonParameters::Nt();

  assert(corr_global.size() == Lt);
  assert(corr_local.size() == Nt);

  std::vector<dcomplex> corr_tmp(Lt, 0);

  int ipe_t = Communicator::ipe(3);

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
#endif

//====================================================================
//============================================================END=====
