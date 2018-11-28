#include "BridgeLib_Private.h"

/*!
        @file    $Id:: fopr_Wilson_eo.cpp #$

        @brief

        @author  UEDA, Satoru  (sueda)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 00:49:46 #$

        @version $LastChangedRevision: 1561 $
*/

#include "fopr_Wilson_eo.h"

#ifdef USE_FACTORY
namespace {
  Fopr *create_object()
  {
    return new Fopr_Wilson_eo();
  }


  Fopr *create_object_with_repr(const std::string& repr)
  {
    return new Fopr_Wilson_eo(repr);
  }


  bool init1 = Fopr::Factory_noarg::Register("Wilson_eo", create_object);
  bool init2 = Fopr::Factory_string::Register("Wilson_eo", create_object_with_repr);
}
#endif
