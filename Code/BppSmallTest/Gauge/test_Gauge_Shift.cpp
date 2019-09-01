/*!
        @file    test_Gauge_Shift.cpp

        @brief

        @author  Hideo Matsufuru  (matsufuru)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate: 2013-01-22 13:51:53 #$

        @version $LastChangedRevision: 1929 $
*/
#include "BppSmallTest.h"
#include "test.h"

#include "IO/gaugeConfig.h"

#include "Measurements/Gauge/staple_lex.h"

#include "Tools/randomNumberManager.h"

//====================================================================
//! Test of gauge quantities.

/*!
                          [28 Dec 2011 H.Matsufuru]
    (Coding history will be recovered from trac.)
    YAML is implemented.  [14 Nov 2012 Y.Namekawa]
    unique_ptr is introduced to avoid memory leaks
                          [21 Mar 2015 Y.Namekawa]
 */

namespace Test_Gauge {
  const std::string test_name = "Gauge.Shift";

  //- test-private parameters
  namespace {
    const std::string filename_input = "test_Gauge_lex.yaml";
  }

  //- prototype declaration
  int shift(void);

#ifdef USE_TESTMANAGER_AUTOREGISTER
  namespace {
#if defined(USE_GROUP_SU2)
    // Nc=2 is not available.
#else
    static const bool is_registered = TestManager::RegisterTest(
      test_name,
      shift
      );
#endif
  }
#endif

  //====================================================================
  int shift(void)
  {
    // #####  parameter setup  #####
    const int Nvol = CommonParameters::Nvol();
    const int Ndim = CommonParameters::Ndim();

    const Parameters params_all = ParameterManager::read(filename_input);

    const Parameters params_test = params_all.lookup("Test_Gauge");

    const string        str_gconf_status = params_test.get_string("gauge_config_status");
    const string        str_gconf_read   = params_test.get_string("gauge_config_type_input");
    const string        readfile         = params_test.get_string("config_filename_input");
    const string        str_rand_type    = params_test.get_string("random_number_type");
    const unsigned long seed             = params_test.get_unsigned_long("seed_for_random_number");
    const string        str_vlevel       = params_test.get_string("verbose_level");

    const bool   do_check        = params_test.is_set("expected_result");
    const double expected_result = do_check ? params_test.get_double("expected_result") : 0.0;

    const Bridge::VerboseLevel vl = vout.set_verbose_level(str_vlevel);

    //- print input parameters
    vout.general(vl, "  gconf_status = %s\n", str_gconf_status.c_str());
    vout.general(vl, "  gconf_read   = %s\n", str_gconf_read.c_str());
    vout.general(vl, "  readfile     = %s\n", readfile.c_str());
    vout.general(vl, "  rand_type    = %s\n", str_rand_type.c_str());
    vout.general(vl, "  seed         = %lu\n", seed);
    vout.general(vl, "  vlevel       = %s\n", str_vlevel.c_str());
    vout.general(vl, "\n");

    //- input parameter check
    int err = 0;
    err += ParameterCheck::non_NULL(str_gconf_status);

    if (err) {
      vout.crucial(vl, "Error at %s: input parameters have not been set\n", test_name.c_str());
      exit(EXIT_FAILURE);
    }


    RandomNumberManager::initialize(str_rand_type, seed);


    // #### object setup #####
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

    const unique_ptr<Index_lex>  index(new Index_lex);
    const unique_ptr<Staple_lex> staple(new Staple_lex);

    ShiftField_lex      shift;
    unique_ptr<Field_G> U2(new Field_G(Nvol, Ndim));

    const unique_ptr<Timer> timer(new Timer(test_name));


    // ####  Execution main part  ####
    timer->start();

    shift.backward(*U2, *U, 3);
    vout.general(vl, "U  = %.8f\n", U->cmp(17, index->site(0, 0, 0, 0), 2));
    vout.general(vl, "U2 = %.8f\n", U2->cmp(17, index->site(0, 0, 0, 7), 2));
    vout.general(vl, "\n");

    shift.forward(*U2, *U, 0);
    vout.general(vl, "U  = %.8f\n", U->cmp(17, index->site(0, 0, 0, 0), 3));
    vout.general(vl, "U2 = %.8f\n", U2->cmp(17, index->site(1, 0, 0, 0), 3));
    vout.general(vl, "\n");

    double    result = 0.0;
    const int nodeid = Communicator::nodeid();

    result = staple->plaquette(*U);
    vout.general(vl, "plaq (original field) = %d  %.8f\n", nodeid, result);

    result = staple->plaquette(*U2);
    vout.general(vl, "plaq (shifted field)  = %d  %.8f\n", nodeid, result);


    timer->report();

    RandomNumberManager::finalize();


    if (do_check) {
      return Test::verify(result, expected_result);
    } else {
      vout.detailed(vl, "check skipped: expected_result not set.\n\n");
      return EXIT_SKIP;
    }
  }
} // namespace Test_Gauge
