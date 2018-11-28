/*!
        @file    $Id:: fopr_Wilson.h #$

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-03-05 11:15:58 #$

        @version $LastChangedRevision: 1584 $
*/

#ifndef FOPR_WILSON_INCLUDED
#define FOPR_WILSON_INCLUDED

//! Wilson fermion operator.

/*!
    This fermion operator defines the standard Wilson fermion.
    The gamma matrix representation is given as control string
    "Dirac"(default) or "Chiral" at the construction, which is
    used to construct the Fopr_Wilson instance.
    The `mode', which of D, Ddag, H, DdagD are multiplied, is
    controlled by setting the pointers to these functions,
    m_mult and m_mult_dag.
    At the beginning, they are set to point mult_undef() which
    just represent the mode has not been set.
    set_mode(string) must be called before mult() is called.
                                    [24 Dec 2011 H,Matsufuru]
    (Coding history will be recovered from trac.)
    YAML is implemented.            [14 Nov 2012 Y.Namekawa]
    unique_ptr is introduced to avoid memory leaks
                                    [21 Mar 2015 Y.Namekawa]
    Use of pimpl pattern removed.
                                    [29 July 2016 T.Aoyama]
 */

// implementations

// original unoptimised version
#ifdef USE_ORG
#include "Org/fopr_Wilson_impl.h"
#endif

// improved version
#ifdef USE_IMP
#include "Imp/fopr_Wilson_impl.h"
#endif

// improved for BG/Q
#ifdef USE_IMP_BGQ
#include "Imp_BGQ/fopr_Wilson_impl.h"
#endif

#ifdef LIB_CPP11

#if defined(USE_IMP_BGQ)
using Fopr_Wilson = Imp_BGQ::Fopr_Wilson;
#elif defined(USE_IMP)
using Fopr_Wilson = Imp::Fopr_Wilson;
#else
using Fopr_Wilson = Org::Fopr_Wilson;
#endif

#else

#if defined(USE_IMP_BGQ)
typedef Imp_BGQ::Fopr_Wilson   Fopr_Wilson;
#elif defined(USE_IMP)
typedef Imp::Fopr_Wilson       Fopr_Wilson;
#else
typedef Org::Fopr_Wilson       Fopr_Wilson;
#endif
#endif
#endif
