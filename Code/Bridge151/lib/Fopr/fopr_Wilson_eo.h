/*!
        @file    fopr_Wilson_eo.h

        @brief

        @author  UEDA, Satoru
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/

#ifndef FOPR_WILSON_EO_INCLUDED
#define FOPR_WILSON_EO_INCLUDED

//! Even-odd Wilson fermion operator.

/*!
    This class BAPI is an even-odd version of Wilson fermion operator.
    At present this is rough implementation, while correctly
    works, and to be updated by supplying complete functionality.
    Only the functions needed for even-odd preconditioned solver
    is ready.
                                     [20 Jun 2012 S.UEDA]
    (Coding history will be recovered from trac.)
    YAML is implemented.             [14 Nov 2012 Y.Namekawa]
    Selector is implemented.         [03 Mar 2013 Y.Namekawa]
    (Selectors are replaced with factories by Aoyama-san)
    Implementation is separated to Fopr_Wilson_eo_impl class BAPI.
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

#ifdef LIB_CPP11

#if defined(USE_IMP)
using Fopr_Wilson_eo = Imp::Fopr_Wilson_eo;
#else
using Fopr_Wilson_eo = Org::Fopr_Wilson_eo;
#endif

#else

#if defined(USE_IMP)
typedef Imp::Fopr_Wilson_eo   Fopr_Wilson_eo;
#else
typedef Org::Fopr_Wilson_eo   Fopr_Wilson_eo;
#endif
#endif

#ifdef USE_FACTORY
namespace Selector_Fopr_Wilson_eo
{
  bool register_factory();
}
#endif

#endif /* FOPR_WILSON_EO_INCLUDED */
