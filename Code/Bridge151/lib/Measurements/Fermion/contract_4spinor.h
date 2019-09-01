/*!
        @file    contract_4spinor.h

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/

#ifndef CONTRACT_4SPINOR_INCLUDED
#define CONTRACT_4SPINOR_INCLUDED

#include "Field/field_F.h"
#include "Tools/gammaMatrix.h"

#include "bridge_complex.h"

//! Contraction of hadron for 4-spinor fermion.

/*!
  This class BAPI calculates contraction of hadron for 4-spinor fermion.
                                 [15 Mar 2012 H.Matsufuru]
  (Coding history will be recovered from trac.)
  Add momentum of sink.        [21 Mar 2015 Y.Namekawa]
  Add functions contract_at_x(), contract_at_y(), contract_at_z() to measure
  spatial correlation functions.
  Add functions contract_at_?() to measure correlation functions directly by corr_global(t)=sum_{xyz}(v2^*)_alpha(t)(gm_sink)_{alpha,beta}(v1)_beta(t).
  Add functions contract_at_?_cos() to measure Fourier transformed correlation functions in terms of cos by corr=sumcos(px)(v2^*)_alpha(t)(gm_sink)_{alpha,beta}(v1)_beta(t).
  Add functions global_corr_x(), global_corr_y(), global_corr_z()
  which should always be used with contract_at_?().
  A function global_corr_t() is moved here from corr2pt_4spinor.h .
  [25 May 2017 Y.Taniguchi]
 */


//! contraction of two Field_F's at a given time t.
//! corr=sum_{xyz}(v2^*)_alpha(t) (gm_sink)_{alpha,beta} (v1)_beta(t)
void contract_at_t(dcomplex& corr,
                   const GammaMatrix& gm_sink,
                   const Field_F& v1, const Field_F& v2,
                   const int time);

//! contraction of two Field_F's and makes correlator in t.
//! corr_global(t)=sum_{xyz}(v2^*)_alpha(t) (gm_sink)_{alpha,beta} (v1)_beta(t)
void contract_at_t(std::vector<dcomplex>& corr_global,
                   const GammaMatrix& gm_sink,
                   const Field_F& v1, const Field_F& v2);

//! contraction for meson at a given time t with Fourier transformation.
//! corr=sum_{xyz}e^{ip_xx+ip_yy+ip_zz}(v2^*)_alpha(t) (gm_sink)_{alpha,beta} (v1)_beta(t), where (p_x,p_y,p_z) is given by momentum_sink.
void contract_at_t(dcomplex& corr,
                   const std::vector<int>& momentum_sink,
                   const GammaMatrix& gm_sink,
                   const std::vector<int>& source_position,
                   const Field_F& v1, const Field_F& v2,
                   const int time);
void contract_at_t(std::vector<dcomplex>& corr_global,
                   const std::vector<int>& momentum_sink,
                   const GammaMatrix& gm_sink,
                   const std::vector<int>& source_position,
                   const Field_F& v1, const Field_F& v2);

//! contraction for meson at a given time t with Fourier transformation.
//! corr=sum_{xyz}cos(p_xx+p_yy+p_zz)(v2^*)_alpha(t) (gm_sink)_{alpha,beta} (v1)_beta(t), where (p_x,p_y,p_z) is given by momentum_sink. This funciton is intended to be used for a Fourier transformation of the real part Re((v2^*)(gm_sink)(v1)) only when imaginary part is unphysical.
void contract_at_t_cos(dcomplex& corr,
                       const std::vector<int>& momentum_sink,
                       const GammaMatrix& gm_sink,
                       const std::vector<int>& source_position,
                       const Field_F& v1, const Field_F& v2,
                       const int time);
void contract_at_t_cos(std::vector<dcomplex>& corr_global,
                       const std::vector<int>& momentum_sink,
                       const GammaMatrix& gm_sink,
                       const std::vector<int>& source_position,
                       const Field_F& v1, const Field_F& v2);

//! contraction for baryon (Nc=3 case only) at a given time t.
void contract_at_t(dcomplex& corr,
                   const GammaMatrix& gm_sink, const int i_alpha,
                   const Field_F& v1, const Field_F& v2, const Field_F& v3,
                   const int time);

//! contraction for meson at a given x.
//! corr=sum_{yzt}(v2^*)_alpha(x) (gm_sink)_{alpha,beta} (v1)_beta(x)
void contract_at_x(dcomplex& corr,
                   const GammaMatrix& gm_sink,
                   const Field_F& v1, const Field_F& v2, int x);
void contract_at_x(std::vector<dcomplex>& corr_global,
                   const GammaMatrix& gm_sink,
                   const Field_F& v1, const Field_F& v2);

//! contraction for meson at a given x with Fourier transformation, where (p_y,p_z,p_t) is given by momentum_sink.
void contract_at_x(dcomplex& corr,
                   const std::vector<int>& momentum_sink,
                   const GammaMatrix& gm_sink,
                   const std::vector<int>& source_position,
                   const Field_F& v1, const Field_F& v2,
                   const int x);
void contract_at_x(std::vector<dcomplex>& corr_global,
                   const std::vector<int>& momentum_sink,
                   const GammaMatrix& gm_sink,
                   const std::vector<int>& source_position,
                   const Field_F& v1, const Field_F& v2);

//! contraction for meson at a given x with Fourier transformation, where (p_y,p_z,p_t) is given by momentum_sink.
void contract_at_x_cos(dcomplex& corr,
                       const std::vector<int>& momentum_sink,
                       const GammaMatrix& gm_sink,
                       const std::vector<int>& source_position,
                       const Field_F& v1, const Field_F& v2,
                       const int x);
void contract_at_x_cos(std::vector<dcomplex>& corr_global,
                       const std::vector<int>& momentum_sink,
                       const GammaMatrix& gm_sink,
                       const std::vector<int>& source_position,
                       const Field_F& v1, const Field_F& v2);

//! contraction for baryon (Nc=3 case only) at a given x.
void contract_at_x(dcomplex& corr,
                   const GammaMatrix& gm_sink, const int i_alpha,
                   const Field_F& v1, const Field_F& v2, const Field_F& v3,
                   const int x);

//! contraction for meson at a given y.
void contract_at_y(dcomplex& corr,
                   const GammaMatrix& gm_sink,
                   const Field_F& v1, const Field_F& v2,
                   const int y);
void contract_at_y(std::vector<dcomplex>& corr_global,
                   const GammaMatrix& gm_sink,
                   const Field_F& v1, const Field_F& v2);

//! contraction for meson at a given y with Fourier transformation, where (p_x,p_z,p_t) is given by momentum_sink.
void contract_at_y(dcomplex& corr,
                   const std::vector<int>& momentum_sink,
                   const GammaMatrix& gm_sink,
                   const std::vector<int>& source_position,
                   const Field_F& v1, const Field_F& v2,
                   const int y);
void contract_at_y(std::vector<dcomplex>& corr_global,
                   const std::vector<int>& momentum_sink,
                   const GammaMatrix& gm_sink,
                   const std::vector<int>& source_position,
                   const Field_F& v1, const Field_F& v2);

//! contraction for meson at a given y with Fourier transformation, where (p_x,p_z,p_t) is given by momentum_sink.
void contract_at_y_cos(dcomplex& corr,
                       const std::vector<int>& momentum_sink,
                       const GammaMatrix& gm_sink,
                       const std::vector<int>& source_position,
                       const Field_F& v1, const Field_F& v2,
                       const int y);
void contract_at_y_cos(std::vector<dcomplex>& corr_global,
                       const std::vector<int>& momentum_sink,
                       const GammaMatrix& gm_sink,
                       const std::vector<int>& source_position,
                       const Field_F& v1, const Field_F& v2);

//! contraction for meson at a given z.
void contract_at_z(dcomplex& corr,
                   const GammaMatrix& gm_sink,
                   const Field_F& v1, const Field_F& v2,
                   const int z);
void contract_at_z(std::vector<dcomplex>& corr_global,
                   const GammaMatrix& gm_sink,
                   const Field_F& v1, const Field_F& v2);

//! contraction for meson at a given z with Fourier transformation, where (p_x,p_y,p_t) is given by momentum_sink.
void contract_at_z(dcomplex& corr,
                   const std::vector<int>& momentum_sink,
                   const GammaMatrix& gm_sink,
                   const std::vector<int>& source_position,
                   const Field_F& v1, const Field_F& v2,
                   const int z);
void contract_at_z(std::vector<dcomplex>& corr_global,
                   const std::vector<int>& momentum_sink,
                   const GammaMatrix& gm_sink,
                   const std::vector<int>& source_position,
                   const Field_F& v1, const Field_F& v2);

//! contraction for meson at a given z with Fourier transformation, where (p_x,p_y,p_t) is given by momentum_sink.
void contract_at_z_cos(dcomplex& corr,
                       const std::vector<int>& momentum_sink,
                       const GammaMatrix& gm_sink,
                       const std::vector<int>& source_position,
                       const Field_F& v1, const Field_F& v2,
                       const int z);
void contract_at_z_cos(std::vector<dcomplex>& corr_global,
                       const std::vector<int>& momentum_sink,
                       const GammaMatrix& gm_sink,
                       const std::vector<int>& source_position,
                       const Field_F& v1, const Field_F& v2);

//! transform node-local correlator in x to global.
void global_corr_x(std::vector<dcomplex>& corr_global,
                   std::vector<dcomplex>& corr_local);

//! transform node-local correlator in y to global.
void global_corr_y(std::vector<dcomplex>& corr_global,
                   std::vector<dcomplex>& corr_local);

//! transform node-local correlator in z to global.
void global_corr_z(std::vector<dcomplex>& corr_global,
                   std::vector<dcomplex>& corr_local);

//! transform node-local correlator in t to global.
void global_corr_t(std::vector<dcomplex>& corr_global,
                   std::vector<dcomplex>& corr_local);

#endif /* CONTRACT_4SPINOR_INCLUDED */
