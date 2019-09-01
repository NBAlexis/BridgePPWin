/*!
        @file    test_HotStart_Eigenvalue.cpp

        @brief

        @author  Hideo Matsufuru  (matsufuru)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate: 2013-01-22 13:51:53 #$

        @version $LastChangedRevision: 1929 $
*/
#include "BppSmallTest.h"
#include "test.h"

#include "IO/gaugeConfig.h"

#include "Tools/eigen_QR_Cmplx.h"
#include "Tools/randomNumbers_Mseries.h"

//====================================================================
//! Test of hot start.

/*!
                               [19 Nov 2013 S.Ueda]
    unique_ptr is introduced to avoid memory leaks
                               [21 Mar 2015 Y.Namekawa]
 */

namespace Test_HotStart {
  const std::string test_name = "HotStart.Eigenvalue";

  //- test-private parameters
  namespace {
    const std::string filename_input = "test_HotStart_Eigenvalue.yaml";
  }

  //- prototype declaration
  int eigenvalue(void);

#ifdef USE_TESTMANAGER_AUTOREGISTER
  namespace {
#if defined(USE_GROUP_SU2)
    // Nc=2 is not available.
#else
    static const bool is_registered = TestManager::RegisterTest(
      test_name,
      eigenvalue
      );
#endif
  }
#endif

  //====================================================================
  int eigenvalue(void)
  {
    const int Nc   = CommonParameters::Nc();
    const int Ndim = CommonParameters::Ndim();
    const int Nvol = CommonParameters::Nvol();
    const int NPE  = CommonParameters::NPE();

    const double pi = 4.0 * atan(1.0);

    const Parameters params_all = ParameterManager::read(filename_input);

    const Parameters params_test = params_all.lookup("Test_HotStart");

    const string str_vlevel = params_test.get_string("verbose_level");

    const bool   do_check        = params_test.is_set("expected_result");
    const double expected_result = do_check ? params_test.get_double("expected_result") : 0.0;

    const Bridge::VerboseLevel vl = vout.set_verbose_level(str_vlevel);


    // ####  Set up a gauge configuration  ####
    unique_ptr<Field_G> U(new Field_G(Nvol, Ndim));

    const int                 i_seed_noise = 1234567;
    unique_ptr<RandomNumbers> rand(new RandomNumbers_Mseries(i_seed_noise));
    U->set_random(rand);


    // #### object setup #####
    const unique_ptr<Timer> timer(new Timer(test_name));


    // ####  Execution main part  ####
    timer->start();

    double                av = 0.0;
    double                va = 0.0;
    Eigen_QR_Cmplx        eigen_qr(Nc);
    std::valarray<double> vec(2 * Nc);

    for (int site = 0; site < Nvol; ++site) {
      for (int mu = 0; mu < Ndim; ++mu) {
        vec = eigen_qr.solve(U->ptr(0, site, mu));
        for (int icolor = 0; icolor < Nc; ++icolor) {
          double arg = atan2(vec[2 * icolor + 1], vec[2 * icolor]);

          vout.detailed(vl, "%10.8f\n", arg);

          av += arg;
          va += arg * arg;
        }
      }
    }

    const double av_all = Communicator::reduce_sum(av);
    const double va_all = Communicator::reduce_sum(va);

    av = av_all / Nvol / NPE / Ndim / Nc;
    va = va_all / Nvol / NPE / Ndim / Nc - av * av;

    va = sqrt(va);

    vout.general(vl, "\n");
    vout.general(vl, "Random SU(%d):\n", Nc);
    // vout.general(vl, "  # of eigenvalue  = %10d\n", Nvol * NPE * Ndim * Nc);
    vout.general(vl, "  ave              = %23.16e\n", av);
    vout.general(vl, "  var              = %23.16e\n", va);
    vout.general(vl, "  var(expected)    = %23.16e\n", pi / sqrt(3.0));

    const double result = av;

    timer->report();


    if (do_check) {
      return Test::verify(result, expected_result);
    } else {
      vout.detailed(vl, "check skipped: expected_result not set.\n\n");
      return EXIT_SKIP;
    }
  }
} // namespace Test_HotStart
