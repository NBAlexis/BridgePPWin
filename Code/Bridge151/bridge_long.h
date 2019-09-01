/*!
        @file    bridge_long.h

        @brief

        @author  Yusuke Namekawa  (namekawa)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/

#ifndef BRIDGE_LONG_INCLUDED
#define BRIDGE_LONG_INCLUDED

//! definition of long for Bridge++

/*!
    This header defines long for Bridge++ on 32/64 bit machine.
                                   [15 Mar 2018 Y.Namekawa]
    Add system bit selection by WORDSIZE, thanks to Aoyama-san
                                   [18 May 2018 Y.Namekawa]
 */

//-- choose long element_type
//- NB. size_t is not recommended due to its environment dependence

#ifdef LIB_CPP11

#include <cstdint>

#if (__WORDSIZE == 64)
typedef int64_t   long_t;
#elif (__WORDSIZE == 32)
#warning 32bit environment may overflow for int variables
typedef int32_t   long_t;
#else
#error unknown environment, not 32bit nor 64bit ?
#endif

#else

// NB. some variables such as Lvol may overflow
//   ex. Lvol = 256^3 x 128 = 2147483648 (max_int = 2147483647)
typedef long long_t;

#endif


#endif // BRIDGE_LONG_INCLUDED
