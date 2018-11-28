/*!
        @file    $Id:: checker.h #$

        @brief

        @author  Tatsumi Aoyama  (aoym)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 00:49:46 #$

        @version $LastChangedRevision: 1561 $
*/

//! checker of input parameters.

/*!
    This checker examines if a parameter satisfies the specified
    condition such as non-zero. The checker returns EXIT_SUCCESS
    or EXIT_FAILURE, instead of bool in Aoyama-san's code.
                                      [16 Jun 2013 Y.Namekawa]
*/

#ifndef CHECKER_INCLUDED
#define CHECKER_INCLUDED

#include "commonParameters.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

namespace ParameterCheck
{
  typedef bool (*valid_double)(const double);
  typedef bool (*valid_int)(const int);
  typedef bool (*valid_double_vector)(const std::vector<double>&);
  typedef bool (*valid_int_vector)(const std::vector<int>&);
  typedef bool (*valid_string)(const std::string&);

  int non_negative(const double v);
  int non_negative(const int v);
  int non_zero(const int v);
  int non_zero(const double v);
  int square_non_zero(const double v);
  int non_NULL(const std::string v);
}
#endif
