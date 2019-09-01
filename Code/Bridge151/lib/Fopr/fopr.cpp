/*!
        @file    fopr.cpp

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/
#include "BridgeLib_Private.h"
#include "fopr.h"

#ifdef USE_FACTORY

#ifdef USE_FACTORY_AUTOREGISTER
#else
// setup factories for all subclasses

#include "fopr_Wilson.h"
#if !USE_IMP
#include "Org/fopr_Wilson_impl.h"
#else
#include "Imp/fopr_Wilson_impl.h"
#endif
#include "fopr_Wilson_eo.h"
#if !USE_IMP
#include "Org/fopr_Wilson_eo_impl.h"
#else
#include "Imp/fopr_Wilson_eo_impl.h"
#endif
#include "fopr_Clover.h"
#include "fopr_Clover_eo.h"
#include "fopr_Rational.h"
#include "fopr_Smeared.h"
#include "fopr_Smeared_eo.h"
#include "fopr_Chebyshev.h"
#include "fopr_WilsonGeneral.h"
#if !USE_IMP
#include "Org/fopr_WilsonGeneral_impl.h"
#else
#include "Imp/fopr_WilsonGeneral_impl.h"
#endif
#include "fopr_CloverGeneral.h"
#include "fopr_Wilson_Isochemical.h"
#include "fopr_Clover_Isochemical.h"
#include "fopr_Wilson_SF.h"
#include "fopr_Clover_SF.h"
#include "fopr_Rational_SF.h"
#include "fopr_CRS.h"


bool Fopr::init_factory()
{
  bool result = true;

  result &= Selector_Fopr_Wilson::register_factory();
#if !USE_IMP
  result &= Org::Fopr_Wilson::register_factory();
#else
  result &= Imp::Fopr_Wilson::register_factory();
#endif
  result &= Selector_Fopr_Wilson_eo::register_factory();
#if !USE_IMP
  result &= Org::Fopr_Wilson_eo::register_factory();
#else
  result &= Imp::Fopr_Wilson_eo::register_factory();
#endif
  result &= Fopr_Clover::register_factory();
  result &= Fopr_Clover_eo::register_factory();
  result &= Fopr_Rational::register_factory();
  result &= Fopr_Smeared::register_factory();
  result &= Fopr_Smeared_eo::register_factory();
  result &= Fopr_Chebyshev::register_factory();
  result &= Selector_Fopr_WilsonGeneral::register_factory();
#if !USE_IMP
  result &= Org::Fopr_WilsonGeneral::register_factory();
#else
  result &= Imp::Fopr_WilsonGeneral::register_factory();
#endif
  result &= Fopr_CloverGeneral::register_factory();
  result &= Fopr_Wilson_Isochemical::register_factory();
  result &= Fopr_Clover_Isochemical::register_factory();
  result &= Fopr_Wilson_SF::register_factory();
  result &= Fopr_Clover_SF::register_factory();
  result &= Fopr_Rational_SF::register_factory();
  result &= Fopr_CRS::register_factory();

  return result;
}


#endif /* USE_FACTORY_AUTOREGISTER */

#endif /* USE_FACTORY */
