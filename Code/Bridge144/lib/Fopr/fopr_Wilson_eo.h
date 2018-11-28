/*!
        @file    $Id:: fopr_Wilson_eo.h #$

        @brief

        @author  UEDA, Satoru
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-03-05 11:15:58 #$

        @version $LastChangedRevision: 1584 $
*/


#ifndef FOPR_WILSON_EO_INCLUDED
#define FOPR_WILSON_EO_INCLUDED

//! Even-odd Wilson fermion operator.

/*!
    This class is an even-odd version of Wilson fermion operator.
    At present this is rough implementation, while correctly
    works, and to be updated by supplying complete functionality.
    Only the functions needed for even-odd preconditioned solver
    is ready.
                                     [20 Jun 2012 S.UEDA]
    (Coding history will be recovered from trac.)
    YAML is implemented.             [14 Nov 2012 Y.Namekawa]
    Selector is implemented.         [03 Mar 2013 Y.Namekawa]
    (Selectors are replaced with factories by Aoyama-san)
    Implementation is separated to Fopr_Wilson_eo_impl class.
                                     [06 Jul 2014 H.Matsufuru]
    unique_ptr is introduced to avoid memory leaks
                                     [21 Mar 2015 Y.Namekawa]
 */

// implementations

// original unoptimised version
#ifdef USE_ORG
#include "Org/fopr_Wilson_eo_impl.h"
#endif

// improved version
#ifdef USE_IMP
#include "Imp/fopr_Wilson_eo_impl.h"
#endif

// improved for BG/Q
#ifdef USE_IMP_BGQ
#include "Imp_BGQ/fopr_Wilson_eo_impl.h"
#endif

#ifdef LIB_CPP11

#if defined(USE_IMP_BGQ)
using Fopr_Wilson_eo = Imp_BGQ::Fopr_Wilson_eo;
#elif defined(USE_IMP)
using Fopr_Wilson_eo = Imp::Fopr_Wilson_eo;
#else
using Fopr_Wilson_eo = Org::Fopr_Wilson_eo;
#endif

#else

#if defined(USE_IMP_BGQ)
typedef Imp_BGQ::Fopr_Wilson_eo   Fopr_Wilson_eo;
#elif defined(USE_IMP)
typedef Imp::Fopr_Wilson_eo       Fopr_Wilson_eo;
#else
typedef Org::Fopr_Wilson_eo       Fopr_Wilson_eo;
#endif
#endif
#endif
