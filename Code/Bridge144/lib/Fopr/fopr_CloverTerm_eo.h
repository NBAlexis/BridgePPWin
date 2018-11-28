/*!
        @file    $Id:: fopr_CloverTerm_eo.h #$

        @brief

        @author  UEDA, Satoru (sueda)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-03-05 11:15:58 #$

        @version $LastChangedRevision: 1584 $
*/

#ifndef FOPR_CLOVERTERM_EO_INCLUDED
#define FOPR_CLOVERTERM_EO_INCLUDED

// implementations

// original unoptimised version
#ifdef USE_ORG
#include "Org/fopr_CloverTerm_eo.h"
#endif

// improved version
#ifdef USE_IMP
#include "Imp/fopr_CloverTerm_eo.h"
#endif

// improved for BG/Q
#ifdef USE_IMP_BGQ
#include "Imp_BGQ/fopr_CloverTerm_eo.h"
#endif

#ifdef LIB_CPP11

#if defined(USE_IMP_BGQ)
using Fopr_CloverTerm_eo = Imp_BGQ::Fopr_CloverTerm_eo;
#elif defined(USE_IMP)
using Fopr_CloverTerm_eo = Imp::Fopr_CloverTerm_eo;
#else
using Fopr_CloverTerm_eo = Org::Fopr_CloverTerm_eo;
#endif

#else

#if defined(USE_IMP_BGQ)
typedef Imp_BGQ::Fopr_CloverTerm_eo   Fopr_CloverTerm_eo;
#elif defined(USE_IMP)
typedef Imp::Fopr_CloverTerm_eo       Fopr_CloverTerm_eo;
#else
typedef Org::Fopr_CloverTerm_eo       Fopr_CloverTerm_eo;
#endif
#endif
#endif
