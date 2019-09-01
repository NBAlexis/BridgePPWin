/*!
        @file    gaugeFixing_None.cpp

        @brief

        @author  Yusuke Namekawa (namekawa)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/
#include "BridgeLib_Private.h"
#include "gaugeFixing_None.h"

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init = GaugeFixing_None::register_factory();
}
#endif

const std::string GaugeFixing_None::class_name = "GaugeFixing_None";

//====================================================================
void GaugeFixing_None::set_parameters(const Parameters& params)
{
  //- No parameters are set.
}


//====================================================================
void GaugeFixing_None::set_parameters(const int Niter, const int Nnaive,
                                      const int Nmeas, const int Nreset,
                                      const double Enorm, const double wp)
{
  //- No parameters are set.
}


//====================================================================
void GaugeFixing_None::fix(Field_G& Ufix, const Field_G& Uorg)
{
  copy(Ufix, Uorg);  // do nothing for gauge fixing, just copy input.
}


//====================================================================
//============================================================END=====
