/*!
        @file    $Id:: fopr_CloverTerm_General.h #$

        @brief

        @author  Yusuke Namekawa (namekawa)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-03-05 11:15:58 #$

        @version $LastChangedRevision: 1584 $
*/

#ifndef FOPR_CLOVERTERM_GENERAL_INCLUDED
#define FOPR_CLOVERTERM_GENERAL_INCLUDED

// implementations

// original unoptimised version
#ifdef USE_ORG
#include "Org/fopr_CloverTerm_General.h"
#endif

// improved version
#ifdef USE_IMP
#include "Imp/fopr_CloverTerm_General.h"
#endif

// improved for BG/Q
#ifdef USE_IMP_BGQ
#include "Imp_BGQ/fopr_CloverTerm_General.h"
#endif

#ifdef LIB_CPP11

#if defined(USE_IMP_BGQ)
using Fopr_CloverTerm_General = Imp_BGQ::Fopr_CloverTerm_General;
#elif defined(USE_IMP)
using Fopr_CloverTerm_General = Imp::Fopr_CloverTerm_General;
#else
using Fopr_CloverTerm_General = Org::Fopr_CloverTerm_General;
#endif

#else

#if defined(USE_IMP_BGQ)
typedef Imp_BGQ::Fopr_CloverTerm_General   Fopr_CloverTerm_General;
#elif defined(USE_IMP)
typedef Imp::Fopr_CloverTerm_General       Fopr_CloverTerm_General;
#else
typedef Org::Fopr_CloverTerm_General       Fopr_CloverTerm_General;
#endif
#endif
#endif
