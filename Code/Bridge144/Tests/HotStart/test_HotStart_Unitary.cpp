/*!
        @file    $Id: test_HotStart_Unitary.cpp #$

        @brief

        @author  Hideo Matsufuru  (matsufuru)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate: 2013-01-22 13:51:53 #$

        @version $LastChangedRevision: 1561 $
*/

#include "Tests/test.h"

#include "IO/gaugeConfig.h"

#include "Tools/randomNumbers_Mseries.h"

//====================================================================
//! Test of hot start.

/*!
                               [19 Nov 2013 S.Ueda]
    unique_ptr is introduced to avoid memory leaks
                               [21 Mar 2015 Y.Namekawa]
 */

namespace Test_HotStart {
  const std::string test_name = "HotStart.Unitary";

  //- test-private parameters
  namespace {
    const std::string filename_input = "test_HotStart_Unitary.yaml";
  }

  //- prototype declaration
  int unitary(void);

#ifdef USE_TESTMANAGER_AUTOREGISTER
  namespace {
#if defined(USE_GROUP_SU2)
    // Nc=2 is not available.
#else
    static const bool is_registered = TestManager::RegisterTest(
      test_name,
      unitary
      );
#endif
  }
#endif

  //====================================================================
  int unitary(void)
  {
    int Nc   = CommonParameters::Nc();
    int Nvol = CommonParameters::Nvol();
    int Ndim = CommonParameters::Ndim();

    Parameters params_all = ParameterManager::read(filename_input);

    Parameters params_test = params_all.lookup("Test_HotStart");

    const string str_vlevel = params_test.get_string("verbose_level");

    const bool   do_check        = params_test.is_set("expected_result");
    const double expected_result = do_check ? params_test.get_double("expected_result") : 0.0;

    Bridge::VerboseLevel vl = vout.set_verbose_level(str_vlevel);


    // ####  Set up a gauge configuration  ####
    unique_ptr<Field_G> U(new Field_G(Nvol, Ndim));

    int i_seed_noise = 1234567;
    unique_ptr<RandomNumbers> rand(new RandomNumbers_Mseries(i_seed_noise));
    U->set_random(rand);


    // #### object setup #####
    unique_ptr<Timer> timer(new Timer(test_name));


    // ####  Execution main part  ####
    timer->start();

    Mat_SU_N link(Nc);
    Mat_SU_N unity(Nc);

    double av = 0;
    for (int site = 0; site < Nvol; ++site) {
      for (int mu = 0; mu < Ndim; ++mu) {
        U->mat(link, site, mu);
        unity.unit();
        unity *= -1.0;

        // |UU^\dag - I|
        unity.multadd_nd(link, link);
        av += sqrt(unity.norm2());
      }
    }

    double av_all = Communicator::reduce_sum(av);
    int    Nlink  = CommonParameters::Lvol() * Ndim;

    vout.general(vl, "\n");
    vout.general(vl, "Random SU(%d):\n", Nc);
    vout.general(vl, "  number of matrix  = %10d\n", Nlink);
    vout.general(vl, "  ave |UU^dag - I| = %23.16e\n", av_all / Nlink);

    double result = av_all / Nlink;

    timer->report();


    if (do_check) {
      return Test::verify(result, expected_result);
    } else {
      vout.detailed(vl, "check skipped: expected_result not set.\n\n");
      return EXIT_SKIP;
    }
  }
} // namespace Test_HotStart
