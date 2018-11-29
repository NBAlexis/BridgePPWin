/*!
@file    $Id: test.cpp #$

@brief

@author  Tatsumi Aoyama  (aoym)
$LastChangedBy: aoym $

@date    $LastChangedDate: 2013-01-22 13:51:53 #$

@version $LastChangedRevision: 1561 $
*/
#include "BppSmallTest.h"

// utility routines for verifying test result.

namespace Test {
    int verify(const double result, const double expected, double eps)
    {
        Bridge::VerboseLevel vl = CommonParameters::Vlevel();

        int verify_result = 0;  // 0 for successful, 1 for failure.

        int precision = default_precision;

#if USE_TESTMANAGER
        // obtain current precision specified through menu.
        precision = TestManager::Instance().getCheckPrecision();
#endif

        if (eps == 0.0) {  // if eps is unspecified
            eps = pow(10.0, -precision);
        }

        vout.general(vl, "\n");
        vout.general(vl, "check test data\n");
        // vout.general(vl, "%s \n", testname.c_str());

        //vout.general(vl, "check precision = %d\n", precision);
        vout.general(vl, "check precision = %2.1e\n", eps);

        vout.general(vl, "result: %15.14le  expected: %15.14le\n", result, expected);

        double diff;
        if (fabs(expected) > eps) {
            diff = (result - expected) / expected;
        }
        else {
            diff = result - expected;
        }

        if (fabs(diff) <= eps) {
            vout.general(vl, "test result: correct! (^_^)\n");
            verify_result = 0;
        }
        else {
            vout.general(vl, "test result: incorrect! (T_T)\n");
            verify_result = 1;
        }
        vout.general(vl, "diff: expected - result = %15.14le\n", diff);
        //vout.general(vl, "diff1: expected - result = %15.14le\n", diff);
        //vout.general(vl, "diff2: exp(-diff) = %15.14le\n", exp(-diff));
        vout.general(vl, "\n");

        return verify_result;
    }
}
