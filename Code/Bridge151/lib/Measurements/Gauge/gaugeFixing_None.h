/*!
        @file    gaugeFixing_None.h

        @brief

        @author  Yusuke Namekawa (namekawa)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/

#ifndef GAUGEFIXING_NONE_INCLUDED
#define GAUGEFIXING_NONE_INCLUDED

#include "gaugeFixing.h"

#include "IO/bridgeIO.h"
using Bridge::vout;


//! None for gauge fixing.

/*
    This class BAPI deals with no gauge fixing to manage spectroscopy
    w/o gauge fixing in a single code, proposed by Aoyama-san.
                                        [30 Jun 2016 Y.Namekawa]
*/

class BAPI GaugeFixing_None : public GaugeFixing
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

#ifdef USE_FACTORY
 private:
  static GaugeFixing *create_object()
  {
    return new GaugeFixing_None();
  }

 public:
  static bool register_factory()
  {
    return GaugeFixing::Factory::Register("None", create_object);
  }
#endif
};
#endif
