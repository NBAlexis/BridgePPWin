/*!
        @file    $Id:: test_Mult.cpp #$

        @brief

        @author  Yusuke Namekawa  (namekawa)
                 $LastChangedBy: namekawa $

        @date    $LastChangedDate:: 2017-02-25 16:51:21 #$

        @version $LastChangedRevision: 1573 $
*/

#include "Tests/test.h"

#include "Fopr/fopr.h"

#include "IO/gaugeConfig.h"

#include "Tools/gammaMatrixSet.h"
#include "Tools/randomNumberManager.h"

//- profiler of fx10
//#include "fj_tool/fapp.h"

//====================================================================
//! Test of Mult.

/*!
                               [06 Jun 2014 Y.Namekawa]
    unique_ptr is introduced to avoid memory leaks
                               [21 Mar 2015 Y.Namekawa]
 */

namespace Test_Mult {
  const std::string test_name = "Mult";

  //- test-private parameters
  namespace {
    // const std::string filename_input  = "test_Mult.yaml";
  }

  //- prototype declaration
  int mult(const std::string&);

  //- mult for various fermions
  int mult_Clover()
  {
    return mult("test_Mult_Clover.yaml");
  }


  int mult_CloverGeneral()
  {
    return mult("test_Mult_CloverGeneral.yaml");
  }


  int mult_Clover_Isochemical()
  {
    return mult("test_Mult_Clover_Isochemical.yaml");
  }


  //- NB. Fopr_Clover_SF is implemented only in Chiral rep.
  // int mult_Clover_SF()
  // {
  //   return mult("test_Mult_Clover_SF.yaml");
  // }


  //- NB. test_Mult_Wilson is implemented separately for beginners
  // int mult_Wilson()
  // {
  //   return mult("test_Mult_Wilson.yaml");
  // }


  int mult_WilsonGeneral()
  {
    return mult("test_Mult_WilsonGeneral.yaml");
  }


  int mult_Wilson_Isochemical()
  {
    return mult("test_Mult_Wilson_Isochemical.yaml");
  }


  //- NB. Fopr_Wilson_SF is implemented only in Chiral rep.
  // int mult_Wilson_SF()
  // {
  //   return mult("test_Mult_Wilson_SF.yaml");
  // }


#ifdef USE_TESTMANAGER_AUTOREGISTER
  namespace {
#if defined(USE_GROUP_SU2)
    // Nc=2 is not available.
#else
    static const bool is_registered_Clover = TestManager::RegisterTest(
      "Mult.Clover",
      mult_Clover
      );

#if defined(USE_IMP_BGQ)
    // Imp_BGQ has not been available for CloverGeneral yet.
#else
    static const bool is_registered_CloverGeneral = TestManager::RegisterTest(
      "Mult.CloverGeneral",
      mult_CloverGeneral
      );
#endif

    static const bool is_registered_Clover_Isochemical = TestManager::RegisterTest(
      "Mult.Clover_Isochemical",
      mult_Clover_Isochemical
      );

    //- NB. Fopr_Clover_SF is implemented only in Chiral rep.
    // static const bool is_registered_Clover_SF = TestManager::RegisterTest(
    //   "Mult.Clover_SF",
    //   mult_Clover_SF
    //   );

    //- NB. test_Mult_Wilson is implemented separately for beginners
    // static const bool is_registered_Wilson = TestManager::RegisterTest(
    //   "Mult.Wilson",
    //   mult_Wilson
    // );

#if defined(USE_IMP_BGQ)
    // Imp_BGQ has not been available for WilsonGeneral yet.
#else
    static const bool is_registered_WilsonGeneral = TestManager::RegisterTest(
      "Mult.WilsonGeneral",
      mult_WilsonGeneral
      );
#endif

    static const bool is_registered_Wilson_Isochemical = TestManager::RegisterTest(
      "Mult.Wilson_Isochemical",
      mult_Wilson_Isochemical
      );

    //- NB. Fopr_Wilson_SF is implemented only in Chiral rep.
    // static const bool is_registered_Wilson_SF = TestManager::RegisterTest(
    //   "Mult.Wilson_SF",
    //   mult_Wilson_SF
    //   );
#endif
  }
#endif

  //====================================================================
  int mult(const std::string& filename_input)
  {
    // ####  parameter setup  ####
    int Nc   = CommonParameters::Nc();
    int Nd   = CommonParameters::Nd();
    int Ndim = CommonParameters::Ndim();
    int Nvol = CommonParameters::Nvol();

    int Lvol    = CommonParameters::Lvol();
    int NPE     = CommonParameters::NPE();
    int Nthread = ThreadManager_OpenMP::get_num_threads_available();

    Parameters params_all = ParameterManager::read(filename_input);

    Parameters params_test = params_all.lookup("Test_Mult");
    Parameters params_fopr = params_all.lookup("Fopr");

    const string        str_gconf_status = params_test.get_string("gauge_config_status");
    const string        str_gconf_read   = params_test.get_string("gauge_config_type_input");
    const string        readfile         = params_test.get_string("config_filename_input");
    const string        str_rand_type    = params_test.get_string("random_number_type");
    const unsigned long seed             = params_test.get_unsigned_long("seed_for_random_number");
    const int           Nmult            = params_test.get_int("number_of_mult");
    const string        str_vlevel       = params_test.get_string("verbose_level");

    const bool   do_check        = params_test.is_set("expected_result");
    const double expected_result = do_check ? params_test.get_double("expected_result") : 0.0;
    const double tolerance       = params_test.get_double("tolerance");

    const string str_fopr_type  = params_fopr.get_string("fermion_type");
    const string str_gmset_type = params_fopr.get_string("gamma_matrix_type");

    Bridge::VerboseLevel vl = vout.set_verbose_level(str_vlevel);

    //- print input parameters
    vout.general(vl, "  gconf_status = %s\n", str_gconf_status.c_str());
    vout.general(vl, "  gconf_read   = %s\n", str_gconf_read.c_str());
    vout.general(vl, "  readfile     = %s\n", readfile.c_str());
    vout.general(vl, "  rand_type    = %s\n", str_rand_type.c_str());
    vout.general(vl, "  seed         = %lu\n", seed);
    vout.general(vl, "  Nmult        = %d\n", Nmult);
    vout.general(vl, "  vlevel       = %s\n", str_vlevel.c_str());
    vout.general(vl, "  gmset_type   = %s\n", str_gmset_type.c_str());

    //- input parameter check
    int err = 0;
    err += ParameterCheck::non_NULL(str_gconf_status);

    if (err) {
      vout.crucial(vl, "Error at %s: input parameters have not been set\n", test_name.c_str());
      exit(EXIT_FAILURE);
    }


    RandomNumberManager::initialize(str_rand_type, seed);


    // ####  Set up a gauge configuration  ####
    unique_ptr<Field_G> U(new Field_G(Nvol, Ndim));

    if (str_gconf_status == "Continue") {
      GaugeConfig(str_gconf_read).read(U, readfile);
    } else if (str_gconf_status == "Cold_start") {
      GaugeConfig("Unit").read(U);
    } else if (str_gconf_status == "Hot_start") {
      GaugeConfig("Random").read(U);
    } else {
      vout.crucial(vl, "Error at %s: unsupported gconf status \"%s\"\n", test_name.c_str(), str_gconf_status.c_str());
      exit(EXIT_FAILURE);
    }


    // ####  object setup  #####
    unique_ptr<GammaMatrixSet> gmset(GammaMatrixSet::New(str_gmset_type));

    unique_ptr<Fopr> fopr(Fopr::New(str_fopr_type, str_gmset_type));
    fopr->set_parameters(params_fopr);
    fopr->set_config(U);
    fopr->set_mode("D");

    unique_ptr<Timer> timer(new Timer(test_name));

    vout.general(vl, "\n");


    // ####  Execution main part  ####
    Field_F b, y;
    b.set(1.0);

    double result = 0.0;


    timer->start();
    //- profiler of fx10 starts
    // fapp_start("Mult",1,1);

#pragma omp parallel
    {
      for (int i = 0; i < Nmult; ++i) {
        fopr->mult(y, b);
      }
      result = y.norm();
    }

    //- profiler of fx10 ends
    // fapp_stop("Mult",1,1);
    timer->stop();
    double elapse_sec = timer->elapsed_sec();


    //- additional verify for Nthread > 1
    int thread_test = 0;

    if (Nthread > 1) {
      double result_single = 0.0;

      fopr->mult(y, b);
      result_single = y.norm();

      if (Test::verify(result, result_single, tolerance)) {
        vout.crucial("%s: result(multithread) not equal to result(single thread)\n", test_name.c_str());
        vout.crucial("%s:   result(multithread)   = %22.14e\n", test_name.c_str(), result);
        vout.crucial("%s:   result(single thread) = %22.14e\n", test_name.c_str(), result_single);

        thread_test = 1;  // test failed.
      }
    }


    //- Flops counting
    double gflo_mult   = fopr->flop_count() / 1.0e+9;
    double gflops_mult = gflo_mult * Nmult / (elapse_sec * NPE * Nthread);

    vout.general(vl, "%s: %12.4e GFlops / core\n", test_name.c_str(), gflops_mult);
    vout.general(vl, "\n");


    RandomNumberManager::finalize();


    if (do_check) {
      int verify_test = Test::verify(result, expected_result, tolerance);

      if ((thread_test == 0) && (verify_test == 0)) {
        return EXIT_SUCCESS;
      } else {
        return EXIT_FAILURE;
      }
    } else {
      vout.detailed(vl, "check skipped: expected_result not set.\n\n");
      return EXIT_SKIP;
    }
  }
} // namespace Test_Mult