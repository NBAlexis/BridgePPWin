/*!
        @file    $Id: test_Rational_Approx.cpp #$

        @brief

        @author  Hideo Matsufuru  (matsufuru)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 00:49:46 #$

        @version $LastChangedRevision: 1561 $
*/

#include "Tests/test.h"

#include "Tools/math_Rational.h"

//====================================================================
//! Test of rational approximation of fermion operators.

/*!
                                [28 Dec 2011 H.Matsufuru]
    (Coding history will be recovered from trac.)
    YAML is implemented.        [14 Nov 2012 Y.Namekawa]
    unique_ptr is introduced to avoid memory leaks
                                [21 Mar 2015 Y.Namekawa]
 */

namespace Test_Rational {
  const std::string test_name = "Rational.Approx";

  //- test-private parameters
  namespace {
    const std::string filename_input = "test_Rational_Approx.yaml";
  }

  //- prototype declaration
  int approx(void);

#ifdef USE_TESTMANAGER_AUTOREGISTER
  namespace {
#if defined(USE_GROUP_SU2)
    // Nc=2 is not available.
#else
    static const bool is_registered = TestManager::RegisterTest(
      test_name,
      approx
      );
#endif
  }
#endif

  //====================================================================
  int approx(void)
  {
    // ####  parameter setup  ####

    Parameters params_all = ParameterManager::read(filename_input);

    Parameters params_test     = params_all.lookup("Test_Rational");
    Parameters params_rational = params_all.lookup("Math_Rational");

    const string str_vlevel = params_test.get_string("verbose_level");

    const bool   do_check        = params_test.is_set("expected_result");
    const double expected_result = do_check ? params_test.get_double("expected_result") : 0.0;

    const int    n_exp  = params_rational.get_int("exponent_numerator");
    const int    d_exp  = params_rational.get_int("exponent_denominator");
    const double x_min  = params_rational.get_double("lower_bound");
    const double x_max  = params_rational.get_double("upper_bound");
    const int    Npoint = params_rational.get_int("number_of_partitions");

    Bridge::VerboseLevel vl = vout.set_verbose_level(str_vlevel);

    //- print input parameters
    vout.general(vl, "  vlevel = %s\n", str_vlevel.c_str());
    vout.general(vl, "  n_exp  = %d\n", n_exp);
    vout.general(vl, "  d_exp  = %d\n", d_exp);
    vout.general(vl, "  x_min  = %12.8f\n", x_min);
    vout.general(vl, "  x_max  = %12.8f\n", x_max);
    vout.general(vl, "  Npoint = %d\n", Npoint);
    vout.general(vl, "\n");


    // ####  object setup  ####
    Math_Rational rational;
    rational.set_parameters(params_rational);

    unique_ptr<Timer> timer(new Timer(test_name));


    // ####  Execution main part  ####
    timer->start();

    double x     = x_min;
    double dx    = (x_max - x_min) / Npoint;
    double p_exp = ((double)n_exp) / ((double)d_exp);

    for (int j = 0; j < Npoint + 1; ++j) {
      double y1 = rational.func(x);
      double y2 = pow(x, p_exp);

      vout.general(vl, "x,rational,x^(n/d) = %12.8f  %12.8f  %12.8f \n", x, y1, y2);

      x += dx;
    }

    double result = x;

    timer->report();


    if (do_check) {
      return Test::verify(result, expected_result);
    } else {
      vout.detailed(vl, "check skipped: expected_result not set.\n\n");
      return EXIT_SKIP;
    }
  }
} // namespace Test_Rational
