/*!
        @file    parameterCheck.cpp

        @brief

        @author  Tatsumi Aoyama  (aoym)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/
#include "BridgeLib_Private.h"
//#include "parameterCheck.h"

namespace ParameterCheck
{
  Bridge::VerboseLevel vl;

  //====================================================================
  int non_negative(const int v)
  {
    if (v < 0) {
      vout.crucial(vl, "ParameterCheck: range check error, negative int.\n");
      return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
  }

  int non_negativef(const double v)
  {
      if (v < 0) {
          vout.crucial(vl, "ParameterCheck: range check error, negative int.\n");
          return EXIT_FAILURE;
      }
      return EXIT_SUCCESS;
  }

  //====================================================================
  int non_zero(const double v)
  {
    if (fabs(v) < CommonParameters::epsilon_criterion()) {
      vout.crucial(vl, "ParameterCheck: range check error, zero double.\n");
      return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
  }


  //====================================================================
  int square_non_zero(const double v)
  {
    if (fabs(v) < CommonParameters::epsilon_criterion2()) {
      vout.crucial(vl, "ParameterCheck: range check error, square_zero double.\n");
      return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
  }


  //====================================================================
  int non_zero(const int v)
  {
    if (v == 0) {
      vout.crucial(vl, "ParameterCheck: range check error, zero int.\n");
      return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
  }


  //====================================================================
  int non_NULL(const std::string v)
  {
    if (v == "NULL") {
      vout.crucial(vl, "ParameterCheck: range check error, NULL string.\n");
      return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
  }


  //====================================================================
  int is_satisfied(const bool cond)
  {
    return cond ? EXIT_SUCCESS : EXIT_FAILURE;
  }
}
