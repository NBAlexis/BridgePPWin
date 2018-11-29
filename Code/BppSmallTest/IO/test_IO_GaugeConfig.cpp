/*!
        @file    $Id: test_IO_GaugeConfig.cpp #$

        @brief

        @author  Tatsumi Aoyama  (aoym)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate: 2013-03-21 15:28:34 #$

        @version $LastChangedRevision: 1593 $
*/

#include "BppSmallTest.h"

//====================================================================
//! Test of I/O.

/*!
    (Implemented by Aoyama-san.)
    (Coding history will be recovered from trac.)
    YAML is implemented.      [14 Nov 2012 Y.Namekawa]
    unique_ptr is introduced to avoid memory leaks
                              [21 Mar 2015 Y.Namekawa]
*/

namespace Test_IO_GaugeConfig {
  const std::string test_name = "IO.GaugeConfig";

  //- test-private parameters
  namespace {
    // const std::string filename_input  = "test_IO_GaugeConfig_Text.yaml";
  }

  //- prototype declaration
  int test_io_gconf_main(const std::string&);
  int check_conf(const unique_ptr<Field_G>& f, const unique_ptr<Field_G>& g);

  //- tests for various file formats
  int test_io_gconf_text()
  {
    return test_io_gconf_main("test_IO_GaugeConfig_Text.yaml");
  }


  int test_io_gconf_binary()
  {
    return test_io_gconf_main("test_IO_GaugeConfig_Binary.yaml");
  }


  int test_io_gconf_binary_parallel()
  {
    return test_io_gconf_main("test_IO_GaugeConfig_BinaryParallel.yaml");
  }


  int test_io_gconf_binary_distributed()
  {
    return test_io_gconf_main("test_IO_GaugeConfig_BinaryDistributed.yaml");
  }


  int test_io_gconf_fortran()
  {
    return test_io_gconf_main("test_IO_GaugeConfig_Fortran.yaml");
  }


  int test_io_gconf_ILDG()
  {
    return test_io_gconf_main("test_IO_GaugeConfig_ILDG.yaml");
  }


  int test_io_gconf_ILDG_parallel()
  {
    return test_io_gconf_main("test_IO_GaugeConfig_ILDG_Parallel.yaml");
  }


#ifdef USE_TESTMANAGER_AUTOREGISTER
  namespace {
#if defined(USE_GROUP_SU2)
    // Nc=2 is not available.
#else
    static const bool is_registered_text = TestManager::RegisterTest(
      "IO.GaugeConfig.Text",
      test_io_gconf_text
      );
    static const bool is_registered_binary = TestManager::RegisterTest(
      "IO.GaugeConfig.Binary",
      test_io_gconf_binary
      );
#ifdef USE_MPI
    static const bool is_registered_binary_parallel = TestManager::RegisterTest(
      "IO.GaugeConfig.BinaryParallel",
      test_io_gconf_binary_parallel
      );
    static const bool is_registered_binary_distributed = TestManager::RegisterTest(
      "IO.GaugeConfig.BinaryDistributed",
      test_io_gconf_binary_distributed
      );
#endif
    static const bool is_registered_fortran = TestManager::RegisterTest(
      "IO.GaugeConfig.Fortran",
      test_io_gconf_fortran
      );

#ifdef USE_LIMELIB
    static const bool is_registered_ILDG = TestManager::RegisterTest(
      "IO.GaugeConfig.ILDG",
      test_io_gconf_ILDG
      );
#ifdef USE_MPI
    static const bool is_registered_ILDG_parallel = TestManager::RegisterTest(
      "IO.GaugeConfig.ILDG_Parallel",
      test_io_gconf_ILDG_parallel
      );
#endif
#endif
#endif
  }
#endif

  //====================================================================
  int test_io_gconf_main(const std::string& filename_input)
  {
    // ####  parameter setup  ####
    int Ndim = CommonParameters::Ndim();
    int Nvol = CommonParameters::Nvol();

    Parameters params_all = ParameterManager::read(filename_input);

    Parameters params_test = params_all.lookup("Test_IO_GaugeConfig");

    const string        str_gconf_status = params_test.get_string("gauge_config_status");
    const string        str_gconf_read   = params_test.get_string("gauge_config_type_input");
    const string        str_gconf_write  = params_test.get_string("gauge_config_type_output");
    const string        readfile         = params_test.get_string("config_filename_input");
    const string        testfile         = params_test.get_string("config_filename_output");
    const string        str_rand_type    = params_test.get_string("random_number_type");
    const unsigned long seed             = params_test.get_unsigned_long("seed_for_random_number");
    const string        str_vlevel       = params_test.get_string("verbose_level");

    const bool   do_check        = params_test.is_set("expected_result");
    const double expected_result = do_check ? params_test.get_double("expected_result") : 0.0;

    Bridge::VerboseLevel vl = vout.set_verbose_level(str_vlevel);

    //- print input parameters
    vout.general(vl, "  gconf_status = %s\n", str_gconf_status.c_str());
    vout.general(vl, "  gconf_read   = %s\n", str_gconf_read.c_str());
    vout.general(vl, "  gconf_write  = %s\n", str_gconf_write.c_str());
    vout.general(vl, "  readfile     = %s\n", readfile.c_str());
    vout.general(vl, "  testfile     = %s\n", testfile.c_str());
    vout.general(vl, "  rand_type    = %s\n", str_rand_type.c_str());
    vout.general(vl, "  seed         = %lu\n", seed);
    vout.general(vl, "  vlevel       = %s\n", str_vlevel.c_str());

    //- input parameter check
    int err = 0;
    err += ParameterCheck::non_NULL(str_gconf_status);
    err += ParameterCheck::non_NULL(str_gconf_write);
    err += ParameterCheck::non_NULL(testfile);

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

    unique_ptr<Field_G>     Utest(new Field_G(Nvol, Ndim));
    unique_ptr<GaugeConfig> gconf_test(new GaugeConfig(str_gconf_write));


    // #### object setup #####
    unique_ptr<Timer> timer(new Timer(test_name));


    // ####  Execution main part  ####
    timer->start();

    err = 0;

    //- first write to file.
    gconf_test->write_file(U, testfile);

    //- then read back from file.
    gconf_test->read_file(Utest, testfile);

    //- check with reference config.
    // bool is_equal = (*Utest == *U);
    bool is_equal = (check_conf(Utest, U) == 0);
    if (!is_equal) ++err;

    vout.general(vl, "%s: \t%s\n", test_name.c_str(), is_equal ? "ok" : "failed");


    int result = err;

    timer->report();

    RandomNumberManager::finalize();


    if (do_check) {
      return Test::verify(result, expected_result);
    } else {
      vout.detailed(vl, "check skipped: expected_result not set.\n\n");
      return EXIT_SKIP;
    }
  }


  //====================================================================
  static inline bool is_equal(const double x, const double y)
  {
    const double eps = CommonParameters::epsilon_criterion();

    if ((x == 0) && (y == 0)) return true;

    if (x == 0) return fabs(y) < eps;

    if (y == 0) return fabs(x) < eps;

    return fabs((x - y) / y) < eps;
  }


  //====================================================================
  int check_conf(const unique_ptr<Field_G>& f, const unique_ptr<Field_G>& g)
  {
    Bridge::VerboseLevel vl = CommonParameters::Vlevel();

    int err = 0;

    int Ndim = CommonParameters::Ndim();
    int Nvol = CommonParameters::Nvol();
    int Nc   = CommonParameters::Nc();

    for (int idir = 0; idir < Ndim; ++idir) {
      for (int isite = 0; isite < Nvol; ++isite) {
        for (int i = 0; i < Nc * Nc; ++i) {
          double v1r = f->cmp_r(i, isite, idir);
          double v1i = f->cmp_i(i, isite, idir);

          double v2r = g->cmp_r(i, isite, idir);
          double v2i = g->cmp_i(i, isite, idir);

//          if (!((v1r == v2r) && (v1i == v2i))) ++err;
          bool is_ok = (is_equal(v1r, v2r) && is_equal(v1i, v2i));
          if (!is_ok) ++err;

          if (!is_ok) {
            vout.general(vl, "%6d : %4d: %2d: %19.15f %19.15f\n                   %19.15f %19.15f : %s\n",
                         isite, i, idir,
                         v1r, v1i, v2r, v2i,
                         is_ok ? "ok" : "fail");
          } else {
            vout.paranoiac(vl, "%6d : %4d: %2d: %19.15f %19.15f\n                   %19.15f %19.15f : %s\n",
                           isite, i, idir,
                           v1r, v1i, v2r, v2i,
                           is_ok ? "ok" : "fail");
          }
        }
      }
    }

    vout.general(vl, "%s: error=%d\n", __func__, err);

    return err;
  }
} // namespace Test_IO_GaugeConfig
