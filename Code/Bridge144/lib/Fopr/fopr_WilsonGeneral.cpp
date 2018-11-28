#include "BridgeLib_Private.h"

/*!
        @file    $Id:: fopr_WilsonGeneral.cpp #$

        @brief

        @author  Yusuke Namekawa (namekawa)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate: 2013-03-21 15:28:34 #$

        @version $LastChangedRevision: 1561 $
*/

#include "fopr_WilsonGeneral.h"

#ifdef USE_FACTORY
namespace {
  Fopr *create_object()
  {
    return new Fopr_WilsonGeneral();
  }


  Fopr *create_object_with_repr(const std::string& repr)
  {
    return new Fopr_WilsonGeneral(repr);
  }


  bool init1 = Fopr::Factory_noarg::Register("WilsonGeneral", create_object);
  bool init2 = Fopr::Factory_string::Register("WilsonGeneral", create_object_with_repr);
}
#endif
