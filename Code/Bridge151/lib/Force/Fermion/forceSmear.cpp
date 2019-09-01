/*!
        @file    forceSmear.cpp

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/
#include "BridgeLib_Private.h"
#include "forceSmear.h"

#ifdef USE_FACTORY

#ifdef USE_FACTORY_AUTOREGISTER
#else
// setup factories for all subclasses

#include "forceSmear_APE.h"
#include "forceSmear_HYP.h"
#include "forceSmear_APE_SF.h"
#include "forceSmear_HYP_SF.h"

bool ForceSmear::init_factory()
{
  bool result = true;

  result &= ForceSmear_APE::register_factory();
  result &= ForceSmear_HYP::register_factory();
  result &= ForceSmear_APE_SF::register_factory();
  result &= ForceSmear_HYP_SF::register_factory();

  return result;
}


#endif /* USE_FACTORY_AUTOREGISTER */

#endif /* USE_FACTORY */
