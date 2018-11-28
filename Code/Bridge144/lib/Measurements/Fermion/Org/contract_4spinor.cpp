#include "BridgeLib_Private.h"
#if USE_ORG

/*!
        @file    $Id:: contract_4spinor.cpp #$

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 18:35:38 #$

        @version $LastChangedRevision: 1571 $
*/

#include "Measurements/Fermion/contract_4spinor.h"
#include "Field/index_lex.h"

#include <cassert>

#if defined USE_GROUP_SU3
#define C1    0
#define C2    1
#define C3    2
#endif

//====================================================================
void contract_at_t(dcomplex& corr,
                   const GammaMatrix& gm_sink,
                   const Field_F& v1, const Field_F& v2,
                   const int time)
{
  int Nc   = CommonParameters::Nc();
  int Nd   = CommonParameters::Nd();
  int Nx   = CommonParameters::Nx();
  int Ny   = CommonParameters::Ny();
  int Nz   = CommonParameters::Nz();
  int Nt   = CommonParameters::Nt();
  int Nvol = CommonParameters::Nvol();

  assert(Nvol == v1.nvol());
  assert(Nvol == v2.nvol());
  assert(time < Nt);

  Index_lex           index;
  std::vector<int>    gm_index(Nd);
  std::vector<double> corr_r(Nd), corr_i(Nd);

  for (int i = 0; i < Nd; ++i) {
    gm_index[i] = gm_sink.index(i);
    corr_r[i]   = 0.0;
    corr_i[i]   = 0.0;
  }

  for (int z = 0; z < Nz; ++z) {
    for (int y = 0; y < Ny; ++y) {
      for (int x = 0; x < Nx; ++x) {
        int site = index.site(x, y, z, time);

        for (int s0 = 0; s0 < Nd; ++s0) {
          int s1 = gm_index[s0];

          for (int c1 = 0; c1 < Nc; ++c1) {
            corr_r[s0] += v1.cmp_r(c1, s1, site) * v2.cmp_r(c1, s0, site)
                          + v1.cmp_i(c1, s1, site) * v2.cmp_i(c1, s0, site);

            corr_i[s0] += -v1.cmp_r(c1, s1, site) * v2.cmp_i(c1, s0, site)
                          + v1.cmp_i(c1, s1, site) * v2.cmp_r(c1, s0, site);
          }
        }
      }
    }
  }

  corr = cmplx(0.0, 0.0);
  for (int s0 = 0; s0 < Nd; ++s0) {
    corr += gm_sink.value(s0) * cmplx(corr_r[s0], corr_i[s0]);
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
  const int Nc   = CommonParameters::Nc();
  const int Nd   = CommonParameters::Nd();
  const int Nx   = CommonParameters::Nx();
  const int Ny   = CommonParameters::Ny();
  const int Nz   = CommonParameters::Nz();
  const int Nt   = CommonParameters::Nt();
  const int Nvol = CommonParameters::Nvol();
  const int Ndim = CommonParameters::Ndim();
  const int Lx   = CommonParameters::Lx();
  const int Ly   = CommonParameters::Ly();
  const int Lz   = CommonParameters::Lz();

  assert(Nvol == v1.nvol());
  assert(Nvol == v2.nvol());
  assert(time < Nt);
  assert(momentum_sink.size() == Ndim - 1);

  Index_lex           index;
  std::vector<int>    gm_index(Nd);
  std::vector<double> corr_r(Nd), corr_i(Nd);

  for (int i = 0; i < Nd; ++i) {
    gm_index[i] = gm_sink.index(i);
    corr_r[i]   = 0.0;
    corr_i[i]   = 0.0;
  }

  static const double PI = 4.0 * atan(1.0);
  std::vector<double> p_unit(Ndim - 1);
  p_unit[0] = (2.0 * PI / Lx) * momentum_sink[0];
  p_unit[1] = (2.0 * PI / Ly) * momentum_sink[1];
  p_unit[2] = (2.0 * PI / Lz) * momentum_sink[2];

  std::vector<int> ipe(Ndim - 1);
  ipe[0] = Communicator::ipe(0);
  ipe[1] = Communicator::ipe(1);
  ipe[2] = Communicator::ipe(2);

  for (int z = 0; z < Nz; ++z) {
    for (int y = 0; y < Ny; ++y) {
      for (int x = 0; x < Nx; ++x) {
        int site = index.site(x, y, z, time);

        int x_global = x + ipe[0] * Nx;
        int y_global = y + ipe[1] * Ny;
        int z_global = z + ipe[2] * Nz;

        double p_x = p_unit[0] * (x_global - source_position[0]);
        double p_y = p_unit[1] * (y_global - source_position[1]);
        double p_z = p_unit[2] * (z_global - source_position[2]);

        double cos_p_xyz = cos(p_x + p_y + p_z);
        double sin_p_xyz = sin(p_x + p_y + p_z);

        for (int s0 = 0; s0 < Nd; ++s0) {
          int s1 = gm_index[s0];

          double v1_v2_r = 0.0;
          double v1_v2_i = 0.0;

          for (int c1 = 0; c1 < Nc; ++c1) {
            v1_v2_r += v1.cmp_r(c1, s1, site) * v2.cmp_r(c1, s0, site)
                       + v1.cmp_i(c1, s1, site) * v2.cmp_i(c1, s0, site);

            v1_v2_i += -v1.cmp_r(c1, s1, site) * v2.cmp_i(c1, s0, site)
                       + v1.cmp_i(c1, s1, site) * v2.cmp_r(c1, s0, site);
          }

          //- corr[s0] += v2^dagger * v1 * exp(i * p_i * x_i)
          corr_r[s0] += v1_v2_r * cos_p_xyz - v1_v2_i * sin_p_xyz;
          corr_i[s0] += v1_v2_r * sin_p_xyz + v1_v2_i * cos_p_xyz;
        }
      }
    }
  }

  corr = cmplx(0.0, 0.0);
  for (int s0 = 0; s0 < Nd; ++s0) {
    corr += gm_sink.value(s0) * cmplx(corr_r[s0], corr_i[s0]);
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
  int Nc   = CommonParameters::Nc();
  int Nd   = CommonParameters::Nd();
  int Nx   = CommonParameters::Nx();
  int Ny   = CommonParameters::Ny();
  int Nz   = CommonParameters::Nz();
  int Nt   = CommonParameters::Nt();
  int Nvol = CommonParameters::Nvol();

  assert(Nc == 3);
  assert(Nvol == v1.nvol());
  assert(Nvol == v2.nvol());
  assert(Nvol == v3.nvol());
  assert(time < Nt);

  Index_lex           index;
  std::vector<int>    gm_index(Nd);
  std::vector<double> c_r(Nd), c_i(Nd);

  for (int i = 0; i < Nd; ++i) {
    gm_index[i] = gm_sink.index(i);
    c_r[i]      = 0.0;
    c_i[i]      = 0.0;
  }

  for (int z = 0; z < Nz; ++z) {
    for (int y = 0; y < Ny; ++y) {
      for (int x = 0; x < Nx; ++x) {
        int site = index.site(x, y, z, time);

        for (int d1 = 0; d1 < Nd; ++d1) {
          int d2 = gm_index[d1];
          int d3 = i_alpha;

          c_r[d1] += (v1.cmp_r(C1, d1, site) * v2.cmp_r(C2, d2, site)
                      - v1.cmp_i(C1, d1, site) * v2.cmp_i(C2, d2, site)) * v3.cmp_r(C3, d3, site)
                     - (v1.cmp_r(C1, d1, site) * v2.cmp_i(C2, d2, site)
                        + v1.cmp_i(C1, d1, site) * v2.cmp_r(C2, d2, site)) * v3.cmp_i(C3, d3, site);
          c_i[d1] += (v1.cmp_r(C1, d1, site) * v2.cmp_r(C2, d2, site)
                      - v1.cmp_i(C1, d1, site) * v2.cmp_i(C2, d2, site)) * v3.cmp_i(C3, d3, site)
                     + (v1.cmp_r(C1, d1, site) * v2.cmp_i(C2, d2, site)
                        + v1.cmp_i(C1, d1, site) * v2.cmp_r(C2, d2, site)) * v3.cmp_r(C3, d3, site);

          c_r[d1] += (v1.cmp_r(C2, d1, site) * v2.cmp_r(C3, d2, site)
                      - v1.cmp_i(C2, d1, site) * v2.cmp_i(C3, d2, site)) * v3.cmp_r(C1, d3, site)
                     - (v1.cmp_r(C2, d1, site) * v2.cmp_i(C3, d2, site)
                        + v1.cmp_i(C2, d1, site) * v2.cmp_r(C3, d2, site)) * v3.cmp_i(C1, d3, site);
          c_i[d1] += (v1.cmp_r(C2, d1, site) * v2.cmp_r(C3, d2, site)
                      - v1.cmp_i(C2, d1, site) * v2.cmp_i(C3, d2, site)) * v3.cmp_i(C1, d3, site)
                     + (v1.cmp_r(C2, d1, site) * v2.cmp_i(C3, d2, site)
                        + v1.cmp_i(C2, d1, site) * v2.cmp_r(C3, d2, site)) * v3.cmp_r(C1, d3, site);

          c_r[d1] += (v1.cmp_r(C3, d1, site) * v2.cmp_r(C1, d2, site)
                      - v1.cmp_i(C3, d1, site) * v2.cmp_i(C1, d2, site)) * v3.cmp_r(C2, d3, site)
                     - (v1.cmp_r(C3, d1, site) * v2.cmp_i(C1, d2, site)
                        + v1.cmp_i(C3, d1, site) * v2.cmp_r(C1, d2, site)) * v3.cmp_i(C2, d3, site);
          c_i[d1] += (v1.cmp_r(C3, d1, site) * v2.cmp_r(C1, d2, site)
                      - v1.cmp_i(C3, d1, site) * v2.cmp_i(C1, d2, site)) * v3.cmp_i(C2, d3, site)
                     + (v1.cmp_r(C3, d1, site) * v2.cmp_i(C1, d2, site)
                        + v1.cmp_i(C3, d1, site) * v2.cmp_r(C1, d2, site)) * v3.cmp_r(C2, d3, site);

          c_r[d1] -= (v1.cmp_r(C3, d1, site) * v2.cmp_r(C2, d2, site)
                      - v1.cmp_i(C3, d1, site) * v2.cmp_i(C2, d2, site)) * v3.cmp_r(C1, d3, site)
                     - (v1.cmp_r(C3, d1, site) * v2.cmp_i(C2, d2, site)
                        + v1.cmp_i(C3, d1, site) * v2.cmp_r(C2, d2, site)) * v3.cmp_i(C1, d3, site);
          c_i[d1] -= (v1.cmp_r(C3, d1, site) * v2.cmp_r(C2, d2, site)
                      - v1.cmp_i(C3, d1, site) * v2.cmp_i(C2, d2, site)) * v3.cmp_i(C1, d3, site)
                     + (v1.cmp_r(C3, d1, site) * v2.cmp_i(C2, d2, site)
                        + v1.cmp_i(C3, d1, site) * v2.cmp_r(C2, d2, site)) * v3.cmp_r(C1, d3, site);

          c_r[d1] -= (v1.cmp_r(C2, d1, site) * v2.cmp_r(C1, d2, site)
                      - v1.cmp_i(C2, d1, site) * v2.cmp_i(C1, d2, site)) * v3.cmp_r(C3, d3, site)
                     - (v1.cmp_r(C2, d1, site) * v2.cmp_i(C1, d2, site)
                        + v1.cmp_i(C2, d1, site) * v2.cmp_r(C1, d2, site)) * v3.cmp_i(C3, d3, site);
          c_i[d1] -= (v1.cmp_r(C2, d1, site) * v2.cmp_r(C1, d2, site)
                      - v1.cmp_i(C2, d1, site) * v2.cmp_i(C1, d2, site)) * v3.cmp_i(C3, d3, site)
                     + (v1.cmp_r(C2, d1, site) * v2.cmp_i(C1, d2, site)
                        + v1.cmp_i(C2, d1, site) * v2.cmp_r(C1, d2, site)) * v3.cmp_r(C3, d3, site);

          c_r[d1] -= (v1.cmp_r(C1, d1, site) * v2.cmp_r(C3, d2, site)
                      - v1.cmp_i(C1, d1, site) * v2.cmp_i(C3, d2, site)) * v3.cmp_r(C2, d3, site)
                     - (v1.cmp_r(C1, d1, site) * v2.cmp_i(C3, d2, site)
                        + v1.cmp_i(C1, d1, site) * v2.cmp_r(C3, d2, site)) * v3.cmp_i(C2, d3, site);
          c_i[d1] -= (v1.cmp_r(C1, d1, site) * v2.cmp_r(C3, d2, site)
                      - v1.cmp_i(C1, d1, site) * v2.cmp_i(C3, d2, site)) * v3.cmp_i(C2, d3, site)
                     + (v1.cmp_r(C1, d1, site) * v2.cmp_i(C3, d2, site)
                        + v1.cmp_i(C1, d1, site) * v2.cmp_r(C3, d2, site)) * v3.cmp_r(C2, d3, site);
        }
      }
    }
  }

  corr = cmplx(0.0, 0.0);
  for (int s0 = 0; s0 < Nd; ++s0) {
    corr += gm_sink.value(s0) * cmplx(c_r[s0], c_i[s0]);
  }
#endif  // (USE_GROUP_SU3)
}


//====================================================================
//============================================================END=====
#endif
