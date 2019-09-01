/*!
        @file    io_format.cpp

        @brief

        @author  Tatsumi Aoyama (aoym)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/
#include "BridgeLib_Private.h"
#include "io_format.h"
#include "io_format_gauge.h"

/**
   Several Data layouts are predefined in namespace IO_Format.
 */

namespace IO_Format {
  namespace {
    const Trivial_Format      trivial_format_ = Trivial_Format();
    const Gauge::ILDG_Format  ildg_format_    = Gauge::ILDG_Format();
    const Gauge::JLQCD_Format jlqcd_format_   = Gauge::JLQCD_Format();
  }

  const Format *Trivial      = &trivial_format_;
  const Format *Gauge::ILDG  = &ildg_format_;
  const Format *Gauge::JLQCD = &jlqcd_format_;
}
