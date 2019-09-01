/*!
        @file    $Id:: test_RandomNumbers_MT19937_GaussianField.cpp #$

        @brief

        @author  Hideo Matsufuru  (matsufuru)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 18:35:38 #$

        @version $LastChangedRevision: 1571 $
*/

#include "Tests/test.h"

#include "Tools/randomNumbers_MT19937.h"

//====================================================================
//! Test of random number generator.

/*!
                                [28 Dec 2011 H.Matsufuru]
    (Coding history will be recovered from trac.)
    YAML is implemented.        [14 Nov 2012 Y.Namekawa]
    unique_ptr is introduced to avoid memory leaks
                                [21 Mar 2015 Y.Namekawa]
 */

namespace Test_RandomNumbers_MT19937 {
  const std::string test_name = "RandomNumbers.MT19937.GaussianField";

  //- test-private parameters
  namespace {
    const std::string filename_input = "test_RandomNumbers_MT19937_GaussianField.yaml";
  }

  //- prototype declaration
  int gaussian_field(void);

#ifdef USE_TESTMANAGER_AUTOREGISTER
  namespace {
    static const bool is_registered = TestManager::RegisterTest(
      test_name,
      gaussian_field
      );
  }
#endif

  //====================================================================
  int gaussian_field(void)
  {
    // ####  parameter setup  ####

    Parameters params_all = ParameterManager::read(filename_input);

    Parameters params_test = params_all.lookup("Test_RandomNumbers");

    int          iseed      = params_test.get_int("int_seed");
    const string str_vlevel = params_test.get_string("verbose_level");

    const bool   do_check        = params_test.is_set("expected_result");
    const double expected_result = do_check ? params_test.get_double("expected_result") : 0.0;

    Bridge::VerboseLevel vl = vout.set_verbose_level(str_vlevel);

    //- print input parameters
    vout.general(vl, "  iseed  = %d\n", iseed);
    vout.general(vl, "  vlevel = %s\n", str_vlevel.c_str());
    vout.general(vl, "\n");


    // ####  object setup  ####
    unique_ptr<RandomNumbers> rand(new RandomNumbers_MT19937(iseed));
    unique_ptr<Timer>         timer(new Timer(test_name));


    // ####  Execution main part  ####
    timer->start();

    int   Nin  = 29;
    int   Nvol = CommonParameters::Nvol();
    int   Nex  = 33;
    Field v(Nin, Nvol, Nex);

    double av = 0.0;
    double vr = 0.0;

    rand->gauss_lex_global(v);

    int size = v.size();
    for (int i = 0; i < size; ++i) {
      av += v.cmp(i);
      vr += v.cmp(i) * v.cmp(i);
      // vout.general(vl, "  %10.8f\n",v.cmp(i));
    }

    double av_all = Communicator::reduce_sum(av);
    double vr_all = Communicator::reduce_sum(vr);

    int global_size = CommonParameters::Lvol() * Nin * Nex;

    av = av_all / global_size;
    vr = vr_all / global_size - av * av;
    vr = sqrt(vr);

    vout.general(vl, "\n");
    vout.general(vl, "Gaussian distribution (Field):\n");
    vout.general(vl, "  number of samples = %10d\n", size);
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
