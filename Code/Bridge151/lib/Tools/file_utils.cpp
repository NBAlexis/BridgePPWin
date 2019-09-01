/*!
        @file    file_utils.cpp

        @brief

        @author  Tatsumi Aoyama  (aoym)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/
#include "BridgeLib_Private.h"
//#include "file_utils.h"

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
