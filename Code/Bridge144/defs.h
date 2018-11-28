/*!
        @file    $Id:: defs.h #$

        @brief

        @author  Tatsumi Aoyama  (aoym)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 18:35:38 #$

        @version $LastChangedRevision: 1571 $
*/

#ifndef DEFS_INCLUDED
#define DEFS_INCLUDED

#include <string>

// for debug
#define LOG printf(">>> %s\n", __PRETTY_FUNCTION__)
//#define LOG

// direction label
enum Direction
{
  XDIR = 0,
  YDIR = 1,
  ZDIR = 2,
  TDIR = 3,
  WDIR = 4,
};

enum ForwardBackward
{
  Forward  =  1,  // +mu
  Backward = -1,  // -mu
};
#endif /* DEFS_INCLUDED */
