#include "BridgeLib_Private.h"

/*!
        @file    $Id:: generatorSet_Mat_SU_N.cpp #$

        @brief

        @author  Hideo Matsufuru  (matsufuru)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 00:49:46 #$

        @version $LastChangedRevision: 1561 $
*/

#include "generatorSet_Mat_SU_N.h"

const std::string GeneratorSet_Mat_SU_N::class_name = "GeneratorSet_Mat_SU_N";

//====================================================================
void GeneratorSet_Mat_SU_N::setup(int Nc)
{
  vout.detailed(m_vl, "%s: SU_N generators setup.\n", class_name.c_str());

  m_Nc  = Nc;
  m_NcA = m_Nc * m_Nc - 1;

  m_Ta.resize(m_NcA);
  for (int ica = 0; ica < m_NcA; ++ica) {
    m_Ta[ica] = new Mat_SU_N(m_Nc);
  }

  int k = 0;
  //  here setup of each generator:
  for (int i = 0; i < m_Nc; ++i) {
    for (int j = i + 1; j < m_Nc; ++j) {
      m_Ta[k]->set(i, j, 1.0, 0.0);
      m_Ta[k]->set(j, i, 1.0, 0.0);

      m_Ta[k + 1]->set(i, j, 0.0, -1.0);
      m_Ta[k + 1]->set(j, i, 0.0, 1.0);

      k += 2;
    }
  }

  for (int i = 1; i < m_Nc; ++i) {
    double elem = sqrt(2.0 / double(i + i * i));
    for (int j = 0; j < i; ++j) {
      m_Ta[k]->set(j, j, elem, 0.0);
    }
    double elem_i = -double(i) * elem;
    m_Ta[k]->set(i, i, elem_i, 0.0);

    ++k;
  }
  assert(k == m_NcA);

  // overall normalization: tr[Ta,Tb] = (1/2) delta_{ab}
  for (int k2 = 0; k2 < m_NcA; ++k2) {
    *m_Ta[k2] *= 0.5;
  }

  //  print();
}


//====================================================================
void GeneratorSet_Mat_SU_N::print()
{
  for (int k = 0; k < m_NcA; ++k) {
    vout.general(m_vl, " generator %d:\n", k);

    for (int i = 0; i < m_Nc; ++i) {
      for (int j = 0; j < m_Nc; ++j) {
        vout.general(m_vl, "   (%7.4f,%7.4f)", m_Ta[k]->r(i, j),
                     m_Ta[k]->i(i, j));
      }
      vout.general(m_vl, "\n");
    }
  }

  vout.general(m_vl, " normalization:\n");
  vout.general(m_vl, "   a   b   tr[T_a T_b]\n");
  Mat_SU_N utmp(m_Nc);
  for (int i = 0; i < m_Nc; ++i) {
    for (int j = 0; j < m_Nc; ++j) {
      utmp.mult_nn(*m_Ta[i], *m_Ta[j]);
      double tr = ReTr(utmp);
      vout.general(m_vl, "   %d   %d    %f\n", i, j, tr);
    }
  }
}


//====================================================================
void GeneratorSet_Mat_SU_N::tidyup()
{
  for (int ica = 0; ica < m_NcA; ++ica) {
    delete m_Ta[ica];
  }

  //  vout.paranoiac(m_vl, "good-bye from GeneratorSet_Mat_SU_N.\n");
}


//====================================================================
//============================================================END=====
