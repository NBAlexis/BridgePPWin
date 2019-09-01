/*!
        @file    action.cpp

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/
#include "BridgeLib_Private.h"
#include "action.h"

#ifdef USE_FACTORY

#ifdef USE_FACTORY_AUTOREGISTER
#else
// setup factories for all subclasses

#include "Gauge/action_G_Plaq.h"
#include "Gauge/action_G_Rectangle.h"
#include "Gauge/action_G_Plaq_SF.h"
#include "Gauge/action_G_Rectangle_SF.h"

bool Action::init_factory()
{
  bool result = true;

  result &= Action_G_Plaq::register_factory();
  result &= Action_G_Rectangle::register_factory();
  result &= Action_G_Plaq_SF::register_factory();
  result &= Action_G_Rectangle_SF::register_factory();

  return result;
}


#endif /* USE_FACTORY_AUTOREGISTER */

#endif /* USE_FACTORY */
