/*!
        @file    fopr_CloverTerm_eo.h

        @brief

        @author  UEDA, Satoru (sueda)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-02-13 20:10:32 #$

        @version $LastChangedRevision: 1942 $
*/

#ifndef FOPR_CLOVERTERM_EO_INCLUDED
#define FOPR_CLOVERTERM_EO_INCLUDED

// implementations

// original unoptimised version
#ifdef USE_ORG
#include "Org/fopr_CloverTerm_eo_impl.h"
#endif

// improved version
#ifdef USE_IMP
#include "Imp/fopr_CloverTerm_eo_impl.h"
#endif

#ifdef LIB_CPP11

#if defined(USE_IMP)
using Fopr_CloverTerm_eo = Imp::Fopr_CloverTerm_eo;
#else
using Fopr_CloverTerm_eo = Org::Fopr_CloverTerm_eo;
#endif

#else

#if defined(USE_IMP)
typedef Imp::Fopr_CloverTerm_eo   Fopr_CloverTerm_eo;
#else
typedef Org::Fopr_CloverTerm_eo   Fopr_CloverTerm_eo;
#endif

#endif

#endif
