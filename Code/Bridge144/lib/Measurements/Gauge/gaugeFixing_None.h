/*!
        @file    $Id:: gaugeFixing_None.h #$

        @brief

        @author  Yusuke Namekawa (namekawa)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 00:49:46 #$

        @version $LastChangedRevision: 1561 $
*/

#ifndef GAUGEFIXING_NONE_INCLUDED
#define GAUGEFIXING_NONE_INCLUDED

#include "gaugeFixing.h"

#include "IO/bridgeIO.h"
using Bridge::vout;


//! None for gauge fixing.

/*
    This class deals with no gauge fixing to manage spectroscopy
    w/o gauge fixing in a single code, proposed by Aoyama-san.
                                        [30 Jun 2016 Y.Namekawa]
*/

class GaugeFixing_None : public GaugeFixing
{
 public:
  static const std::string class_name;

 public:
  GaugeFixing_None() : GaugeFixing() {}

  ~GaugeFixing_None() {}

  void set_parameters(const Parameters& params);
  void set_parameters(const int Niter, const int Nnaive,
                      const int Nmeas, const int Nreset,
                      const double Enorm, const double wp);

  void fix(Field_G& Ufix, const Field_G& Uorg);
};
#endif
