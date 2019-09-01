/*!
        @file    bridge_defs.h

        @brief

        @author  Tatsumi Aoyama  (aoym)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/

#ifndef DEFS_INCLUDED
#define DEFS_INCLUDED

#include <string>

// for debug
//#define LOG printf(">>> %s\n", __PRETTY_FUNCTION__)
#define LOG

// direction label
enum BAPI Direction
{
  XDIR = 0,
  YDIR = 1,
  ZDIR = 2,
  TDIR = 3,
  WDIR = 4
};

enum BAPI ForwardBackward
{
  Forward  =  1,  // +mu
  Backward = -1   // -mu
};

namespace Element_type
{
  enum BAPI element_type
  {
    REAL = 1, COMPLEX = 2
  };
}

#endif /* DEFS_INCLUDED */
