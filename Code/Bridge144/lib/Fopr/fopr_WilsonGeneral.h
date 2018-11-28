/*!
        @file    $Id:: fopr_WilsonGeneral.h #$

        @brief

        @author  Yusuke Namekawa (namekawa)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-03-05 11:15:58 #$

        @version $LastChangedRevision: 1584 $
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

// improved for BG/Q
#ifdef USE_IMP_BGQ
#include "Imp_BGQ/fopr_WilsonGeneral_impl.h"
#endif

#ifdef LIB_CPP11

#if defined(USE_IMP_BGQ)
using Fopr_WilsonGeneral = Imp_BGQ::Fopr_WilsonGeneral;
#elif defined(USE_IMP)
using Fopr_WilsonGeneral = Imp::Fopr_WilsonGeneral;
#else
using Fopr_WilsonGeneral = Org::Fopr_WilsonGeneral;
#endif

#else

#if defined(USE_IMP_BGQ)
typedef Imp_BGQ::Fopr_WilsonGeneral   Fopr_WilsonGeneral;
#elif defined(USE_IMP)
typedef Imp::Fopr_WilsonGeneral       Fopr_WilsonGeneral;
#else
typedef Org::Fopr_WilsonGeneral       Fopr_WilsonGeneral;
#endif
#endif
#endif
