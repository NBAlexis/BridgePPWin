/*!
        @file    fopr_Wilson_eo.cpp

        @brief

        @author  UEDA, Satoru  (sueda)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/
#include "BridgeLib_Private.h"
#include "fopr_Wilson_eo.h"

#ifdef USE_FACTORY

namespace Selector_Fopr_Wilson_eo
{
  namespace {
    Fopr *create_object()
    {
      return new Fopr_Wilson_eo();
    }


    Fopr *create_object_with_repr(const std::string& repr)
    {
      return new Fopr_Wilson_eo(repr);
    }
  }

  bool register_factory()
  {
    bool init1 = Fopr::Factory_noarg::Register("Wilson_eo", create_object);
    bool init2 = Fopr::Factory_string::Register("Wilson_eo", create_object_with_repr);

    return init1 && init2;
  }


#ifdef USE_FACTORY_AUTOREGISTER
  namespace {
    bool init = register_factory();
  }
#endif
}
#endif
