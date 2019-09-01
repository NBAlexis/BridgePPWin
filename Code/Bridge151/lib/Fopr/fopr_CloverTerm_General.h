/*!
        @file    fopr_CloverTerm_General.h

        @brief

        @author  Yusuke Namekawa (namekawa)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-02-13 20:10:32 #$

        @version $LastChangedRevision: 1942 $
*/

#ifndef FOPR_CLOVERTERM_GENERAL_INCLUDED
#define FOPR_CLOVERTERM_GENERAL_INCLUDED

// implementations

// original unoptimised version
#ifdef USE_ORG
#include "Org/fopr_CloverTerm_General_impl.h"
#endif

// improved version
#ifdef USE_IMP
#include "Imp/fopr_CloverTerm_General_impl.h"
#endif

#ifdef LIB_CPP11

#if defined(USE_IMP)
using Fopr_CloverTerm_General = Imp::Fopr_CloverTerm_General;
#else
using Fopr_CloverTerm_General = Org::Fopr_CloverTerm_General;
#endif

#else

#if defined(USE_IMP)
typedef Imp::Fopr_CloverTerm_General   Fopr_CloverTerm_General;
#else
typedef Org::Fopr_CloverTerm_General   Fopr_CloverTerm_General;
#endif

#endif

#endif
