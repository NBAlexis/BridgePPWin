#include "BridgeLib_Private.h"

/*!
        @file    $Id:: file_utils.cpp #$

        @brief

        @author  Tatsumi Aoyama  (aoym)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 00:49:46 #$

        @version $LastChangedRevision: 1561 $
*/

#include "file_utils.h"

//====================================================================
std::string FileUtils::generate_filename(const char *fmt, ...)
{
  static const int buf_size = FILENAME_MAX;
  static char      buf[buf_size]; // not thread-safe.

  va_list arg;

  va_start(arg, fmt);
  vsnprintf(buf, buf_size, fmt, arg);
  va_end(arg);

  return std::string(buf);
}


//==========================================================
//==================================================END=====
