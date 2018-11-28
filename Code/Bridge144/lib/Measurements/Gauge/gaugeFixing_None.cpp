#include "BridgeLib_Private.h"

/*!
        @file    $Id:: gaugeFixing_None.cpp #$

        @brief

        @author  Yusuke Namekawa (namekawa)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2017-11-13 18:52:56 #$

        @version $LastChangedRevision: 1677 $
*/

#include "gaugeFixing_None.h"

#ifdef USE_FACTORY
namespace {
  GaugeFixing *create_object()
  {
    return new GaugeFixing_None();
  }

  bool init = GaugeFixing::Factory::Register("None", create_object);
}
#endif

const std::string GaugeFixing_None::class_name = "GaugeFixing_None";

//====================================================================
void GaugeFixing_None::set_parameters(const Parameters& )
{
  //- No parameters are set.
}


//====================================================================
void GaugeFixing_None::set_parameters(const int , const int ,
                                      const int , const int ,
                                      const double , const double )
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
