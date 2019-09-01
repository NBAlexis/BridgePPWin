/*!
        @file    gammaMatrixSet.cpp

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/
#include "BridgeLib_Private.h"
#include "gammaMatrixSet.h"

#ifdef USE_FACTORY

#ifdef USE_FACTORY_AUTOREGISTER
#else
// setup factories for all subclasses

#include "gammaMatrixSet_Dirac.h"
#include "gammaMatrixSet_Chiral.h"

bool GammaMatrixSet::init_factory()
{
  bool result = true;

  result &= GammaMatrixSet_Dirac::register_factory();
  result &= GammaMatrixSet_Chiral::register_factory();

  return result;
}


#endif /* USE_FACTORY_AUTOREGISTER */

#endif /* USE_FACTORY */
