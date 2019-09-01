/*!
        @file    fopr_WilsonGeneral.cpp

        @brief

        @author  Yusuke Namekawa (namekawa)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate: 2013-03-21 15:28:34 #$

        @version $LastChangedRevision: 1929 $
*/
#include "BridgeLib_Private.h"
#include "fopr_WilsonGeneral.h"

#ifdef USE_FACTORY

namespace Selector_Fopr_WilsonGeneral
{
  namespace {
    Fopr *create_object()
    {
      return new Fopr_WilsonGeneral();
    }


    Fopr *create_object_with_repr(const std::string& repr)
    {
      return new Fopr_WilsonGeneral(repr);
    }
  }

  bool register_factory()
  {
    bool init1 = Fopr::Factory_noarg::Register("WilsonGeneral", create_object);
    bool init2 = Fopr::Factory_string::Register("WilsonGeneral", create_object_with_repr);

    return init1 && init2;
  }


#ifdef USE_FACTORY_AUTOREGISTER
  namespace {
    bool init = register_factory();
  }
#endif
}
#endif
