/*!
        @file    $Id:: configure.h #$

        @brief

        @author  Tatsumi Aoyama  (aoym)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 00:49:46 #$

        @version $LastChangedRevision: 1561 $
*/

#ifndef CONFIGURE_INCLUDED
#define CONFIGURE_INCLUDED

// #define BRIDGE_VERSION "1.2rc"
#define BRIDGE_VERSION    "1.4.0"

// some #define options are set as compiler options from Makefile

//#define USE_MPI
#define ENABLE_MULTI_INSTANCE

// smart pointer support

#ifdef cpp11_available
#undef cpp11_available
#endif

#if defined(__INTEL_COMPILER)
// intel compiler does not set __cplusplus macro properly
#if (__INTEL_COMPILER >= 1400)
#define cpp11_available
#endif

#else

#if defined (__cplusplus) && __cplusplus >= 201103L
#define cpp11_available
#endif
#endif

//see:
//https://stackoverflow.com/questions/47043869/is-c11-available-in-visual-studio-2017
//cpp11 is for vs2017
#define LIB_CPP11
#define cpp11_available

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
