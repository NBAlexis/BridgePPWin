/*!
        @file    $Id:: fopr_CloverTerm.h #$

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-03-05 11:15:58 #$

        @version $LastChangedRevision: 1584 $
*/

#ifndef FOPR_CLOVERTERM_INCLUDED
#define FOPR_CLOVERTERM_INCLUDED

//! Clover term operator.

/*!
    This class implements the clover term for the clover (improved
    Wilson) fermion operator.
    This part was separated from the Fopr_Clover class.
    The field strength is calculate when the function
    set_config() is called.
    The `mode' for setting fermion operator mode is now only
    defined to the case 'D'.
                [30 Sep 2012 H.Matsufuru,
                 original clover operator: 24 Dec 2011 H.M.]
    (Coding history will be recovered from trac.)
    YAML is implemented.            [14 Nov 2012 Y.Namekawa]
    Selector is implemented.        [03 Mar 2013 Y.Namekawa]
    (Selectors are replaced with factories by Aoyama-san)
    unique_ptr is introduced to avoid memory leaks
                                    [21 Mar 2015 Y.Namekawa]
 */

// implementations

#ifdef USE_ORG
#include "Org/fopr_CloverTerm.h"
#endif

#ifdef USE_IMP
#include "Imp/fopr_CloverTerm.h"
#endif

#ifdef USE_IMP_BGQ
#include "Imp_BGQ/fopr_CloverTerm.h"
#endif

#ifdef LIB_CPP11

#if defined(USE_IMP_BGQ)
using Fopr_CloverTerm = Imp_BGQ::Fopr_CloverTerm;
#elif defined(USE_IMP)
using Fopr_CloverTerm = Imp::Fopr_CloverTerm;
#else
using Fopr_CloverTerm = Org::Fopr_CloverTerm;
#endif

#else

#if defined(USE_IMP_BGQ)
typedef Imp_BGQ::Fopr_CloverTerm   Fopr_CloverTerm;
#elif defined(USE_IMP)
typedef Imp::Fopr_CloverTerm       Fopr_CloverTerm;
#else
typedef Org::Fopr_CloverTerm       Fopr_CloverTerm;
#endif
#endif
#endif
