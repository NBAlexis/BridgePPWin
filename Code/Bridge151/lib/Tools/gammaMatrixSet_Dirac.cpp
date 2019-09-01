/*!
        @file    gammaMatrixSet_Dirac.cpp

        @brief

        @author  Hideo Matsufuru  (matsufuru)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/
#include "BridgeLib_Private.h"
#include "gammaMatrixSet_Dirac.h"

#include <cassert>

#include "Parameters/commonParameters.h"
#include "Communicator/communicator.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init = GammaMatrixSet_Dirac::register_factory();
}
#endif

const std::string GammaMatrixSet_Dirac::class_name = "GammaMatrixSet_Dirac";

//====================================================================
void GammaMatrixSet_Dirac::init_GM()
{
  vout.general(m_vl, "Gamma matrix: dirac representation.\n");

  m_gm[UNITY].set(0, 0, icomplex(1, 0));
  m_gm[UNITY].set(1, 1, icomplex(1, 0));
  m_gm[UNITY].set(2, 2, icomplex(1, 0));
  m_gm[UNITY].set(3, 3, icomplex(1, 0));

  m_gm[GAMMA1].set(0, 3, icomplex(0, -1));
  m_gm[GAMMA1].set(1, 2, icomplex(0, -1));
  m_gm[GAMMA1].set(2, 1, icomplex(0, 1));
  m_gm[GAMMA1].set(3, 0, icomplex(0, 1));

  m_gm[GAMMA2].set(0, 3, icomplex(-1, 0));
  m_gm[GAMMA2].set(1, 2, icomplex(1, 0));
  m_gm[GAMMA2].set(2, 1, icomplex(1, 0));
  m_gm[GAMMA2].set(3, 0, icomplex(-1, 0));

  m_gm[GAMMA3].set(0, 2, icomplex(0, -1));
  m_gm[GAMMA3].set(1, 3, icomplex(0, 1));
  m_gm[GAMMA3].set(2, 0, icomplex(0, 1));
  m_gm[GAMMA3].set(3, 1, icomplex(0, -1));

  m_gm[GAMMA4].set(0, 0, icomplex(1, 0));
  m_gm[GAMMA4].set(1, 1, icomplex(1, 0));
  m_gm[GAMMA4].set(2, 2, icomplex(-1, 0));
  m_gm[GAMMA4].set(3, 3, icomplex(-1, 0));

  m_gm[GAMMA5].set(0, 2, icomplex(1, 0));
  m_gm[GAMMA5].set(1, 3, icomplex(1, 0));
  m_gm[GAMMA5].set(2, 0, icomplex(1, 0));
  m_gm[GAMMA5].set(3, 1, icomplex(1, 0));

  m_gm[GAMMA51] = m_gm[GAMMA5].mult(m_gm[GAMMA1]);
  m_gm[GAMMA52] = m_gm[GAMMA5].mult(m_gm[GAMMA2]);
  m_gm[GAMMA53] = m_gm[GAMMA5].mult(m_gm[GAMMA3]);
  m_gm[GAMMA54] = m_gm[GAMMA5].mult(m_gm[GAMMA4]);

  m_gm[GAMMA15] = m_gm[GAMMA1].mult(m_gm[GAMMA5]);
  m_gm[GAMMA25] = m_gm[GAMMA2].mult(m_gm[GAMMA5]);
  m_gm[GAMMA35] = m_gm[GAMMA3].mult(m_gm[GAMMA5]);
  m_gm[GAMMA45] = m_gm[GAMMA4].mult(m_gm[GAMMA5]);

  m_gm[SIGMA12] = m_gm[GAMMA2].mult_i(m_gm[GAMMA1]);
  m_gm[SIGMA23] = m_gm[GAMMA3].mult_i(m_gm[GAMMA2]);
  m_gm[SIGMA31] = m_gm[GAMMA1].mult_i(m_gm[GAMMA3]);

  m_gm[SIGMA41] = m_gm[GAMMA1].mult_i(m_gm[GAMMA4]);
  m_gm[SIGMA42] = m_gm[GAMMA2].mult_i(m_gm[GAMMA4]);
  m_gm[SIGMA43] = m_gm[GAMMA3].mult_i(m_gm[GAMMA4]);

  m_gm[CHARGECONJG] = m_gm[GAMMA4].mult(m_gm[GAMMA2]);

  //print();
}


//====================================================================
void GammaMatrixSet_Dirac::print()
{
  vout.general(m_vl, "\n  unity =\n");
  m_gm[UNITY].print();

  vout.general(m_vl, "\n  gamma1 =\n");
  m_gm[GAMMA1].print();

  vout.general(m_vl, "\n  gamma2 =\n");
  m_gm[GAMMA2].print();

  vout.general(m_vl, "\n  gamma3 =\n");
  m_gm[GAMMA3].print();

  vout.general(m_vl, "\n  gamma4 =\n");
  m_gm[GAMMA4].print();

  vout.general(m_vl, "\n  gamma5 =\n");
  m_gm[GAMMA5].print();

  vout.general(m_vl, "\n  gamma5 * gamma1 =\n");
  m_gm[GAMMA51].print();

  vout.general(m_vl, "\n  gamma5 * gamma2 =\n");
  m_gm[GAMMA52].print();

  vout.general(m_vl, "\n  gamma5 * gamma3 =\n");
  m_gm[GAMMA53].print();

  vout.general(m_vl, "\n  gamma5 * gamma4 =\n");
  m_gm[GAMMA54].print();

  vout.general(m_vl, "\n  sigma12 = -(i/2) [gamma1, gamma2] =\n");
  m_gm[SIGMA12].print();

  vout.general(m_vl, "\n  sigma23 = -(i/2) [gamma2, gamma3] =\n");
  m_gm[SIGMA23].print();

  vout.general(m_vl, "\n  sigma31 = -(i/2) [gamma3, gamma1] =\n");
  m_gm[SIGMA31].print();

  vout.general(m_vl, "\n  sigma41 = -(i/2) [gamma4, gamma1] =\n");
  m_gm[SIGMA41].print();

  vout.general(m_vl, "\n  sigma42 = -(i/2) [gamma4, gamma2] =\n");
  m_gm[SIGMA42].print();

  vout.general(m_vl, "\n  sigma43 = -(i/2) [gamma4, gamma3] =\n");
  m_gm[SIGMA43].print();

  vout.general(m_vl, "\n  charge_conjg = gamma4 * gamma2 =\n");
  m_gm[CHARGECONJG].print();
}


//====================================================================
//============================================================END=====
