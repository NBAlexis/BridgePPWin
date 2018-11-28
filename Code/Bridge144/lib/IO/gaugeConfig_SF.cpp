#include "BridgeLib_Private.h"

/*!
        @file    $Id:: gaugeConfig_SF.cpp #$

        @brief

        @author  Yusuke Taniguchi
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 00:49:46 #$

        @version $LastChangedRevision: 1561 $
*/

#include "gaugeConfig_SF.h"

/*!
  Setup the tree level gauge configuration for the SF boundary condition, which minimize the classical action.
 */

const std::string GaugeConfig_SF::class_name = "GaugeConfig_SF";

//====================================================================
void GaugeConfig_SF::set_cold_SF(Field *U, double *phi, double *phipr)
{
  int Nc   = CommonParameters::Nc();
  int Nvol = CommonParameters::Nvol();
  int Ndim = CommonParameters::Ndim();
  int Nx   = CommonParameters::Nx();
  int Ny   = CommonParameters::Ny();
  int Nz   = CommonParameters::Nz();
  int Nt   = CommonParameters::Nt();
  int Lx   = CommonParameters::Lx();
  int Ly   = CommonParameters::Ly();
  int Lz   = CommonParameters::Lz();
  int Lt   = CommonParameters::Lt();
  int NPEt = CommonParameters::NPEt();

  Index_lex idx;
  int       site;
  Field_G   *up = (Field_G *)U;
  Mat_SU_N  ut(Nc);

  ut.unit();
  Mat_SU_N utree(Nc);
  utree.zero();
  double sr, si;

  for (int t = 0; t < Nt; t++) {
    sr = cos((phi[0] + (phipr[0] - phi[0]) / Lt * (t + (Communicator::ipe(3)) * Nt)) / Lx);
    si = sin((phi[0] + (phipr[0] - phi[0]) / Lt * (t + (Communicator::ipe(3)) * Nt)) / Lx);
    utree.set(0, 0, sr, si);
    sr = cos((phi[1] + (phipr[1] - phi[1]) / Lt * (t + (Communicator::ipe(3)) * Nt)) / Lx);
    si = sin((phi[1] + (phipr[1] - phi[1]) / Lt * (t + (Communicator::ipe(3)) * Nt)) / Lx);
    utree.set(1, 1, sr, si);
    sr = cos((phi[2] + (phipr[2] - phi[2]) / Lt * (t + (Communicator::ipe(3)) * Nt)) / Lx);
    si = sin((phi[2] + (phipr[2] - phi[2]) / Lt * (t + (Communicator::ipe(3)) * Nt)) / Lx);
    utree.set(2, 2, sr, si);
    for (int z = 0; z < Nz; z++) {
      for (int y = 0; y < Ny; y++) {
        for (int x = 0; x < Nx; x++) {
          site = idx.site(x, y, z, t);
          for (int mu = 0; mu < Ndim - 1; ++mu) {
            up->set_mat(site, mu, utree);
          }
        }
      }
    }
  }
  {
    int mu = Ndim - 1;
    for (int site = 0; site < Nvol; ++site) {
      up->set_mat(site, mu, ut);
    }
  }

  vout.detailed(m_vl, "%s: gauge config. was set to the tree level value with SF BC.\n", class_name.c_str());
}


//====================================================================
//============================================================END=====
