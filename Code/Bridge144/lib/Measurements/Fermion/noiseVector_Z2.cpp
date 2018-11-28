#include "BridgeLib_Private.h"

/*!
        @file    $Id:: noiseVector_Z2.cpp #$

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 18:35:38 #$

        @version $LastChangedRevision: 1571 $
*/

#include "noiseVector_Z2.h"

const std::string NoiseVector_Z2::class_name = "NoiseVector_Z2";

//====================================================================
void NoiseVector_Z2::set(Field& v)
{
  // This implementation assumes the given field v is complex field.

  m_rand->uniform_lex_global(v);

  int Nex  = v.nex();
  int Nvol = v.nvol();
  int Nin  = v.nin();

  double RF2 = 1.0 / sqrt(2.0);

  for (int ex = 0; ex < Nex; ++ex) {
    for (int site = 0; site < Nvol; ++site) {
      for (int in = 0; in < Nin; ++in) {
        double rn1 = v.cmp(in, site, ex);
        double rn2 = floor(2.0 * rn1);
        double rn3 = (2.0 * rn2 - 1.0) * RF2;
        v.set(in, site, ex, rn3);
      }
    }
  }
}


//====================================================================
//============================================================END=====
