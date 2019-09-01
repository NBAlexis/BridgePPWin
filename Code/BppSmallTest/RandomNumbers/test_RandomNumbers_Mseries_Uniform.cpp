/*!
        @file    test_RandomNumbers_Mseries_Uniform.cpp

        @brief

        @author  Hideo Matsufuru  (matsufuru)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/
#include "BppSmallTest.h"
#include "test.h"

#include "Tools/randomNumbers_Mseries.h"

//====================================================================
//! Test of random number generator.

/*!
                                [28 Dec 2011 H.Matsufuru]
    (Coding history will be recovered from trac.)
    YAML is implemented.        [14 Nov 2012 Y.Namekawa]
    unique_ptr is introduced to avoid memory leaks
                                [21 Mar 2015 Y.Namekawa]
 */

namespace Test_RandomNumbers_Mseries {
  const std::string test_name = "RandomNumbers.Mseries.Uniform";

  //- test-private parameters
  namespace {
    const std::string filename_input = "test_RandomNumbers_Mseries_Uniform.yaml";
  }

  //- prototype declaration
  int uniform_calc_pi(void);

#ifdef USE_TESTMANAGER_AUTOREGISTER
  namespace {
    static const bool is_registered = TestManager::RegisterTest(
      test_name,
      uniform_calc_pi
      );
  }
#endif

  //====================================================================
  int uniform_calc_pi(void)
  {
    // ####  parameter setup  ####

    const Parameters params_all = ParameterManager::read(filename_input);

    const Parameters params_test = params_all.lookup("Test_RandomNumbers");

    const int    Nseed      = params_test.get_int("number_of_seeds");
    const int    iseed_base = params_test.get_int("int_seed_base");
    const int    Nrand      = params_test.get_int("number_of_samples");
    const string str_vlevel = params_test.get_string("verbose_level");

    const bool   do_check        = params_test.is_set("expected_result");
    const double expected_result = do_check ? params_test.get_double("expected_result") : 0.0;

    const Bridge::VerboseLevel vl = vout.set_verbose_level(str_vlevel);

    //- print input parameters
    vout.general(vl, "  Nseed      = %d\n", Nseed);
    vout.general(vl, "  iseed_base = %d\n", iseed_base);
    vout.general(vl, "  Nrand      = %d\n", Nrand);
    vout.general(vl, "  vlevel     = %s\n", str_vlevel.c_str());
    vout.general(vl, "\n");


    // #### object setup #####
    const unique_ptr<Timer> timer(new Timer(test_name));


    // ####  Execution main part  ####
    timer->start();

    vout.general(vl, "\n");
    vout.general(vl, "Monte Carlo estimate of pi:\n");
    vout.general(vl, "  number of samples = %10d\n", Nrand);
    vout.general(vl, "        seed    estimate of pi\n");

    double t1 = 0.0;
    double t2 = 0.0;
    for (int iseed = 0; iseed < Nseed; ++iseed) {
      int iseed2 = iseed_base + iseed;

      unique_ptr<RandomNumbers> rand(new RandomNumbers_Mseries(iseed2));

      int Npi = 0;
      for (int i = 0; i < Nrand; ++i) {
        double rand1 = rand->get();
        double rand2 = rand->get();
        double r     = rand1 * rand1 + rand2 * rand2;

        if (r < 1.0) { ++Npi; }
        //  vout.general(vl, "  %10.8f  %10.8f\n",rand1,rand2);
      }

      double pi_exp = (4.0 * Npi) / Nrand;

      t1 += pi_exp;
      t2 += pi_exp * pi_exp;

      //vout.general(vl, "  estimate of pi    = %10.8f\n",pi_exp);
      vout.general(vl, "  %10d    %14.10f\n", iseed2, pi_exp);
    }

    const double api = t1 / (double)Nseed;
    const double vpi = t2 / (double)Nseed - api * api;
    const double dpi = sqrt(vpi);
    const double epi = dpi / sqrt((double)Nseed - 1);

    const double pi = 3.141592653589793;
    vout.general(vl, "  true value = %10.8f\n", pi);
    vout.general(vl, "  average    = %10.8f\n", api);
    vout.general(vl, "  variance   = %10.8f\n", vpi);
    vout.general(vl, "  deviation  = %10.8f\n", dpi);
    vout.general(vl, "  error      = %10.8f\n", epi);

    const double result = api;
    // changed to check the obtained value of pi itseft. [25 May 2014 H.M.]

    timer->report();


    if (do_check) {
      return Test::verify(result, expected_result);
    } else {
      vout.detailed(vl, "check skipped: expected_result not set.\n\n");
      return EXIT_SKIP;
    }
  }
} // namespace Test_RandomNumbers
