/*!
        @file    gaugeFixing.cpp

        @brief

        @author  Yusuke Namekawa (namekawa)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/
#include "BridgeLib_Private.h"
#include "gaugeFixing.h"

#ifdef USE_FACTORY

#ifdef USE_FACTORY_AUTOREGISTER
#else
// setup factories for all subclasses

#include "gaugeFixing_None.h"
#include "gaugeFixing_Coulomb.h"
#include "gaugeFixing_Landau.h"

bool GaugeFixing::init_factory()
{
  bool result = true;

  result &= GaugeFixing_None::register_factory();
  result &= GaugeFixing_Coulomb::register_factory();
  result &= GaugeFixing_Landau::register_factory();

  return result;
}


#endif /* USE_FACTORY_AUTOREGISTER */

#endif /* USE_FACTORY */
