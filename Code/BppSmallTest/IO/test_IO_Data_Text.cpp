/*!
        @file    test_IO_Data_Text.cpp

        @brief

        @author  Tatsumi Aoyama  (aoym)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate: 2013-03-21 15:28:34 #$

        @version $LastChangedRevision: 1929 $
*/
#include "BppSmallTest.h"
#include "test.h"

#include "IO/dataIO_Text.h"
#include "IO/gaugeConfig.h"

#include "Tools/randomNumberManager.h"

//====================================================================
//! Test of I/O.

/*!
    (Implemented by Aoyama-san.)
    (Coding history will be recovered from trac.)
    YAML is implemented.      [14 Nov 2012 Y.Namekawa]
    unique_ptr is introduced to avoid memory leaks
                              [21 Mar 2015 Y.Namekawa]
*/

namespace Test_IO_Data {
  const std::string test_name = "IO.Data.Text";

  //- test-private parameters
  namespace {
    const std::string filename_input = "test_IO_Data_Text.yaml";
  }

  //- prototype declaration
  int test_io_data_text(void);

#ifdef USE_TESTMANAGER_AUTOREGISTER
  namespace {
#if defined(USE_GROUP_SU2)
    // Nc=2 is not available.
#else
    static const bool is_registered = TestManager::RegisterTest(
      test_name,
      test_io_data_text
      );
#endif
  }
#endif

  //====================================================================
  int test_io_data_text(void)
  {
    // ####  parameter setup  ####
    const int    Ndim = CommonParameters::Ndim();
    const int    Nvol = CommonParameters::Nvol();
    const double tiny = CommonParameters::epsilon_criterion();

    const Parameters params_all = ParameterManager::read(filename_input);

    const Parameters params_test = params_all.lookup("Test_IO_Data");

    const string        str_gconf_status = params_test.get_string("gauge_config_status");
    const string        str_gconf_read   = params_test.get_string("gauge_config_type_input");
    const string        readfile         = params_test.get_string("config_filename_input");
    const string        testfile         = params_test.get_string("config_filename_output");
    const string        str_rand_type    = params_test.get_string("random_number_type");
    const unsigned long seed             = params_test.get_unsigned_long("seed_for_random_number");
    const int           data_size        = params_test.get_int("data_size");
    const string        str_vlevel       = params_test.get_string("verbose_level");

    const bool   do_check        = params_test.is_set("expected_result");
    const double expected_result = do_check ? params_test.get_double("expected_result") : 0.0;

    const Bridge::VerboseLevel vl = vout.set_verbose_level(str_vlevel);

    //- print input parameters
    vout.general(vl, "  gconf_status = %s\n", str_gconf_status.c_str());
    vout.general(vl, "  gconf_read   = %s\n", str_gconf_read.c_str());
    vout.general(vl, "  readfile     = %s\n", readfile.c_str());
    vout.general(vl, "  testfile     = %s\n", testfile.c_str());
    vout.general(vl, "  rand_type    = %s\n", str_rand_type.c_str());
    vout.general(vl, "  seed         = %lu\n", seed);
    vout.general(vl, "  data_size    = %d\n", data_size);
    vout.general(vl, "  vlevel       = %s\n", str_vlevel.c_str());
    vout.general(vl, "\n");

    //- input parameter check
    int err = 0;
    err += ParameterCheck::non_NULL(str_gconf_status);
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


    // #####  object setup  #####
    const unique_ptr<DataIO> data_io(new DataIO_Text);
    const unique_ptr<Timer>  timer(new Timer(test_name));


    // ####  Execution main part  ####
    timer->start();

    err = 0;

    //- first write to file.
    std::vector<double> array(data_size);
    double              *p = U->ptr(0);

    for (size_t i = 0; i < data_size; ++i) {
      array[i] = *p++;
    }

    data_io->write_file(array, testfile, false);


    //- write again with complex data.
    std::vector<dcomplex> arrayc(data_size);
    p = U->ptr(0);

    for (size_t i = 0; i < data_size; ++i) {
      double x = *p++;
      double y = *p++;

      arrayc[i] = cmplx(x, y);
    }

    data_io->write_file(arrayc, testfile, true);


    //- then read back from file.
    std::vector<double> array2(data_size);

    data_io->read_file(array2, testfile);


    //- check with reference config.
    for (size_t i = 0; i < data_size; ++i) {
      if (abs(array[i] - array2[i]) > tiny) ++err;
    }

    vout.general(vl, "%s: \t%s\n", test_name.c_str(), (err == 0) ? "ok" : "failed");


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
} // namespace Test_IO_Data
