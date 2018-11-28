#include "BridgeLib_Private.h"

/*!
        @file    $Id:: checker.cpp #$

        @brief

        @author  Tatsumi Aoyama  (aoym)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 00:49:46 #$

        @version $LastChangedRevision: 1561 $
*/

#include "checker.h"

namespace ParameterCheck
{
  Bridge::VerboseLevel vl;


  int non_negative(const int v)
  {
    if (v < 0) {
      vout.crucial(vl, "ParameterCheck: range check error, negative int.\n");
      return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
  }


  int non_zero(const double v)
  {
    if (fabs(v) < CommonParameters::epsilon_criterion()) {
      vout.crucial(vl, "ParameterCheck: range check error, zero double.\n");
      return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
  }


  int square_non_zero(const double v)
  {
    if (fabs(v) < CommonParameters::epsilon_criterion2()) {
      vout.crucial(vl, "ParameterCheck: range check error, square_zero double.\n");
      return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
  }


  int non_zero(const int v)
  {
    if (v == 0) {
      vout.crucial(vl, "ParameterCheck: range check error, zero int.\n");
      return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
  }


  int non_NULL(const std::string v)
  {
    if (v == "NULL") {
      vout.crucial(vl, "ParameterCheck: range check error, NULL string.\n");
      return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
  }
}
