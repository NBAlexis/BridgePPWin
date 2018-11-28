#include "BridgeLib_Private.h"

/*!
        @file    $Id:: force_F.cpp #$

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 00:49:46 #$

        @version $LastChangedRevision: 1561 $
*/

#include "force_F.h"

// default templates for core and core1
//====================================================================
void Force::force_core(Field& force_, const Field& eta)
{
  int Nvol = CommonParameters::Nvol();
  int Ndim = CommonParameters::Ndim();

  Field_G force(Nvol, Ndim);

  force_udiv(force, eta);

  mult_generator(force);

  force_ = force;
}


//====================================================================
void Force::force_core1(Field& force_, const Field& zeta, const Field& eta)
{
  int Nvol = CommonParameters::Nvol();
  int Ndim = CommonParameters::Ndim();

  Field_G force(Nvol, Ndim);

  force_udiv1(force, zeta, eta);

  mult_generator(force);

  force_ = force;
}


// utility function
//====================================================================
void Force::mult_generator(Field_G& force)
{
//  int Nvol = CommonParameters::Nvol();
//  int Ndim = CommonParameters::Ndim();
  int Nvol = force.nvol();
  int Ndim = force.nex();

  for (int mu = 0; mu < Ndim; ++mu) {
    for (int isite = 0, nsite = Nvol; isite < nsite; ++isite) {
      Mat_SU_N u = m_U->mat(isite, mu);

      u *= force.mat(isite, mu);
      u.at();
      u *= -2.0;

      force.set_mat(isite, mu, u);
    }
  }
}


//====================================================================
