/*!
        @file    configure.h

        @brief

        @author  Tatsumi Aoyama  (aoym)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-08-09 10:37:08 #$

        @version $LastChangedRevision: 1987 $
*/

#ifndef CONFIGURE_INCLUDED
#define CONFIGURE_INCLUDED

#define BRIDGE_VERSION    "1.5.1"

// some #define options are set as compiler options from Makefile

//#define USE_MPI
#define ENABLE_MULTI_INSTANCE

// smart pointer support

#ifdef cpp11_available
#undef cpp11_available
#endif

#if defined(__INTEL_COMPILER)
// Intel compiler does not set __cplusplus macro properly
#if (__INTEL_COMPILER >= 1400)
#define cpp11_available
#endif

#else

//#if defined (__cplusplus) && __cplusplus >= 201103L
//In MSVC, 199711L is cpp11
#if defined (__cplusplus) && __cplusplus >= 199711L
#define cpp11_available
#endif
#endif

#if defined(LIB_CPP11) && defined(cpp11_available)
// use c++11 unique_ptr
#include <memory>
using std::unique_ptr;
#else
// use bridge alternative
#include "Tools/unique_pointer.h"
using Bridge::unique_ptr;
#endif

#undef cpp11_available
#endif /* CONFIGURE_INCLUDED */
