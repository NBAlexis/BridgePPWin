/*!
        @file    $Id:: file_utils.h #$

        @brief

        @author  Tatsumi Aoyama  (aoym)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 00:49:46 #$

        @version $LastChangedRevision: 1561 $
*/

#ifndef FILE_UTILS_INCLUDED
#define FILE_UTILS_INCLUDED

#include <cstdio>
#include <cstdlib>
#include <cstdarg>

#include "Parameters/commonParameters.h"

//! File utility.

/*!
    FileUtils provides file utilities, made by Aoyama-san.
                                   [28 May 2013 Y.Namekawa]
 */

namespace FileUtils
{
  std::string generate_filename(const char *fmt, ...);
}
#endif
