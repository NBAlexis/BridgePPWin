/*!
        @file    test_HotStart_Determinant.cpp

        @brief

        @author  Hideo Matsufuru  (matsufuru)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate: 2013-01-22 13:51:53 #$

        @version $LastChangedRevision: 1929 $
*/
#include "BppSmallTest.h"
#include "test.h"

#include "IO/gaugeConfig.h"

#include "Tools/decompose_LUP_Cmplx.h"
#include "Tools/randomNumbers_Mseries.h"

//====================================================================
//! Test of hot start.

/*!
                               [19 Nov 2013 S.Ueda]
    unique_ptr is introduced to avoid memory leaks
                               [21 Mar 2015 Y.Namekawa]
 */

namespace Test_HotStart {
  const std::string test_name = "HotStart.Determinant";

  //- test-private parameters
  namespace {
    const std::string filename_input = "test_HotStart_Determinant.yaml";
  }

  //- prototype declaration
  int determinant(void);

#ifdef USE_TESTMANAGER_AUTOREGISTER
  namespace {
#if defined(USE_GROUP_SU2)
    // Nc=2 is not available.
#else
    static const bool is_registered = TestManager::RegisterTest(
      test_name,
      determinant
      );
#endif
  }
#endif

  //====================================================================
  int determinant(void)
  {
    // #####  parameter setup  #####
    const int Nc   = CommonParameters::Nc();
    const int Ndim = CommonParameters::Ndim();
    const int Nvol = CommonParameters::Nvol();
    const int NPE  = CommonParameters::NPE();

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

    double re_av = 0.0;
    double re_va = 0.0;
    double im_av = 0.0;
    double im_va = 0.0;

    Decompose_LUP_Cmplx lup(Nc);
    for (int site = 0; site < Nvol; ++site) {
      for (int mu = 0; mu < Ndim; ++mu) {
        lup.set_matrix(U->ptr(0, site, mu));

        dcomplex det = lup.determinant();

        re_av += real(det);
        re_va += real(det) * real(det);
        im_av += imag(det);
        im_va += imag(det) * imag(det);
      }
    }

    const double re_av_all = Communicator::reduce_sum(re_av);
    const double re_va_all = Communicator::reduce_sum(re_va);
    const double im_av_all = Communicator::reduce_sum(im_av);
    const double im_va_all = Communicator::reduce_sum(im_va);

    re_av = re_av_all / Nvol / NPE / Ndim;
    im_av = im_av_all / Nvol / NPE / Ndim;
    re_va = re_va_all / Nvol / NPE / Ndim - re_av * re_av;
    im_va = im_va_all / Nvol / NPE / Ndim - im_av * im_av;

    re_va = sqrt(re_va);
    im_va = sqrt(im_va);

    vout.general(vl, "\n");
    vout.general(vl, "Random SU(%d):\n", Nc);
    // vout.general(vl, "  number of matrix  = %10d\n", Nvol * NPE * Ndim);
    vout.general(vl, "  ave Re(det)       = %23.16e\n", re_av);
    vout.general(vl, "  var Re(det)       = %23.16e\n", re_va);
    vout.general(vl, "  ave Im(det)       = %23.16e\n", im_av);
    vout.general(vl, "  var Im(det)       = %23.16e\n", im_va);

    const double result = re_av;

    timer->report();


    if (do_check) {
      return Test::verify(result, expected_result);
    } else {
      vout.detailed(vl, "check skipped: expected_result not set.\n\n");
      return EXIT_SKIP;
    }
  }
} // namespace Test_HotStart
