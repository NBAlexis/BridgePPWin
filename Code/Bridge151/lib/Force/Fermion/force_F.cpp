/*!
        @file    force_F.cpp

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/
#include "BridgeLib_Private.h"
#include "force_F.h"

// default templates for core and core1
//====================================================================
void Force::force_core(Field& force_, const Field& eta)
{
  const int Nvol = CommonParameters::Nvol();
  const int Ndim = CommonParameters::Ndim();

  Field_G force(Nvol, Ndim);

  force_udiv(force, eta);

  mult_generator(force);

  force_ = force;
}


//====================================================================
void Force::force_core1(Field& force_, const Field& zeta, const Field& eta)
{
  const int Nvol = CommonParameters::Nvol();
  const int Ndim = CommonParameters::Ndim();

  Field_G force(Nvol, Ndim);

  force_udiv1(force, zeta, eta);

  mult_generator(force);

  force_ = force;
}


// utility function
//====================================================================
void Force::mult_generator(Field_G& force)
{
  const int Nvol = force.nvol();
  const int Ndim = force.nex();

  for (int mu = 0; mu < Ndim; ++mu) {
    for (int isite = 0; isite < Nvol; ++isite) {
      Mat_SU_N u = m_U->mat(isite, mu);

      u *= force.mat(isite, mu);
      u.at();
      u *= -2.0;

      force.set_mat(isite, mu, u);
    }
  }
}


//====================================================================
