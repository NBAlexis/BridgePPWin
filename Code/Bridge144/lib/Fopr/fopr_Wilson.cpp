#include "BridgeLib_Private.h"

/*!
        @file    $Id:: fopr_Wilson.cpp #$

        @brief

        @author  Hideo Matsufuru  (matsufuru)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate: 2013-03-21 15:28:34 #$

        @version $LastChangedRevision: 1561 $
*/

#include "fopr_Wilson.h"

#ifdef USE_FACTORY
namespace {
  Fopr *create_object()
  {
    return new Fopr_Wilson();
  }


  Fopr *create_object_with_repr(const std::string& repr)
  {
    return new Fopr_Wilson(repr);
  }


  bool init1 = Fopr::Factory_noarg::Register("Wilson", create_object);
  bool init2 = Fopr::Factory_string::Register("Wilson", create_object_with_repr);
}
#endif
