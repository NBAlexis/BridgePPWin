/*!
        @file    $Id: test.h #$

        @brief

        @author  Tatsumi Aoyama  (aoym)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate: 2013-03-21 16:21:38 #$

        @version $LastChangedRevision: 1561 $
*/

#ifndef TEST_INCLUDED
#define TEST_INCLUDED

#define EXIT_SKIP    -1

namespace Test {
  static const int default_precision        = 11;
  static const int default_output_precision = 15;

  // utility routines for verifying test result.
  //   returns 0 if result and expected agree within specified criterion
  //   where the criterion is 10^-precision.

  int verify(const double result, const double expected, double eps = 0.0);
}
#endif
