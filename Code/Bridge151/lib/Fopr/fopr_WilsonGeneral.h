/*!
        @file    fopr_WilsonGeneral.h

        @brief

        @author  Yusuke Namekawa (namekawa)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/

#ifndef FOPR_WILSONGENERAL_INCLUDED
#define FOPR_WILSONGENERAL_INCLUDED

// implementations

// original unoptimised version
#ifdef USE_ORG
#include "Org/fopr_WilsonGeneral_impl.h"
#endif

// improved version
#ifdef USE_IMP
#include "Imp/fopr_WilsonGeneral_impl.h"
#endif

#ifdef LIB_CPP11

#if defined(USE_IMP)
using Fopr_WilsonGeneral = Imp::Fopr_WilsonGeneral;
#else
using Fopr_WilsonGeneral = Org::Fopr_WilsonGeneral;
#endif

#else

#if defined(USE_IMP)
typedef Imp::Fopr_WilsonGeneral   Fopr_WilsonGeneral;
#else
typedef Org::Fopr_WilsonGeneral   Fopr_WilsonGeneral;
#endif
#endif

#ifdef USE_FACTORY
namespace Selector_Fopr_WilsonGeneral
{
  bool register_factory();
}
#endif

#endif /* FOPR_WILSONGENERAL_INCLUDED */
