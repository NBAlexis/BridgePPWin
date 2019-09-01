/*!
        @file    source.cpp

        @brief

        @author  Tatsumi Aoyama  (aoym)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/
#include "BridgeLib_Private.h"
#include "source.h"

#ifdef USE_FACTORY

#ifdef USE_FACTORY_AUTOREGISTER
#else
// setup factories for all subclasses

#include "source_Local.h"
#include "source_Wall.h"
#include "source_Exponential.h"
#include "source_MomentumWall.h"
#include "source_Random.h"

bool Source::init_factory()
{
  bool result = true;

  result &= Source_Local::register_factory();
  result &= Source_Wall::register_factory();
  result &= Source_Exponential::register_factory();
  result &= Source_MomentumWall::register_factory();
  result &= Source_Random::register_factory();

  return result;
}


#endif /* USE_FACTORY_AUTOREGISTER */

#endif /* USE_FACTORY */
