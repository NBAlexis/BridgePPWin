#include "BridgeLib_Private.h"

/*!
        @file    $Id:: io_format.cpp #$

        @brief

        @author  Tatsumi Aoyama (aoym)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 00:49:46 #$

        @version $LastChangedRevision: 1561 $
*/
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
