/*!
        @file    $Id:: test_RandomNumbers_Mseries_Gaussian.cpp #$

        @brief

        @author  Hideo Matsufuru  (matsufuru)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 18:35:38 #$

        @version $LastChangedRevision: 1571 $
*/

#include "Tests/test.h"

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
  const std::string test_name = "RandomNumbers.Mseries.Gaussian";

  //- test-private parameters
  namespace {
    const std::string filename_input = "test_RandomNumbers_Mseries_Gaussian.yaml";
  }

  //- prototype declaration
  int gaussian(void);

#ifdef USE_TESTMANAGER_AUTOREGISTER
  namespace {
    static const bool is_registered = TestManager::RegisterTest(
      test_name,
      gaussian
      );
  }
#endif

  //====================================================================
  int gaussian(void)
  {
    // ####  parameter setup  ####

    Parameters params_all = ParameterManager::read(filename_input);

    Parameters params_test = params_all.lookup("Test_RandomNumbers");

    int          iseed      = params_test.get_int("int_seed");
    int          Nrand      = params_test.get_int("number_of_samples");
    const string str_vlevel = params_test.get_string("verbose_level");

    const bool   do_check        = params_test.is_set("expected_result");
    const double expected_result = do_check ? params_test.get_double("expected_result") : 0.0;

    Bridge::VerboseLevel vl = vout.set_verbose_level(str_vlevel);

    //- print input parameters
    vout.general(vl, "  iseed  = %d\n", iseed);
    vout.general(vl, "  Nrand  = %d\n", Nrand);
    vout.general(vl, "  vlevel = %s\n", str_vlevel.c_str());
    vout.general(vl, "\n");


    // ####  object setup  ####
    unique_ptr<RandomNumbers> rand(new RandomNumbers_Mseries(iseed));
    unique_ptr<Timer>         timer(new Timer(test_name));


    // ####  Execution main part  ####
    timer->start();

    double av = 0.0;
    double vr = 0.0;

    double rand1, rand2;
    for (int i = 0; i < Nrand; ++i) {
      rand->gauss(rand1, rand2);
      av += rand1 + rand2;
      vr += rand1 * rand1 + rand2 * rand2;
      // vout.general(vl, "  %10.8f  %10.8f\n",rand1,rand2);
    }
    av = av / (2.0 * Nrand);
    vr = vr / (2.0 * Nrand) - av * av;
    vr = sqrt(vr);

    vout.general(vl, "\n");
    vout.general(vl, "Gaussian distribution:\n");
    vout.general(vl, "  number of samples = %10d\n", Nrand);
    vout.general(vl, "  average           = %10.8f\n", av);
    vout.general(vl, "  variance          = %10.8f\n", vr);
    vout.general(vl, "  variance(expect)  = %10.8f\n", 1.0 / sqrt(2.0));

    double result = vr;

    timer->report();


    if (do_check) {
      return Test::verify(result, expected_result);
    } else {
      vout.detailed(vl, "check skipped: expected_result not set.\n\n");
      return EXIT_SKIP;
    }
  }
} // namespace Test_RandomNumbers
