/*!
        @file    $Id: test_WilsonLoop.cpp #$

        @brief

        @author  Hideo Matsufuru  (matsufuru)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate: 2013-01-22 13:51:53 #$

        @version $LastChangedRevision: 1571 $
*/

#include "BppSmallTest.h"

//====================================================================
//! Test of Wilson loop measurement.

/*!
                          [27 May 2012 H.Matsufuru]
    (Coding history will be recovered from trac.)
    YAML is implemented.  [14 Nov 2012 Y.Namekawa]
    unique_ptr is introduced to avoid memory leaks
                          [21 Mar 2015 Y.Namekawa]
 */

namespace Test_WilsonLoop {
  const std::string test_name = "WilsonLoop";

  //- test-private parameters
  namespace {
    const std::string filename_input = "test_WilsonLoop.yaml";
  }

  //- prototype declaration
  int wilsonloop(void);

#ifdef USE_TESTMANAGER_AUTOREGISTER
  namespace {
#if defined(USE_GROUP_SU2)
    // Nc=2 is not available.
#else
    static const bool is_registered = TestManager::RegisterTest(
      test_name,
      wilsonloop
      );
#endif
  }
#endif

  //====================================================================
  int wilsonloop(void)
  {
    // ####  parameter setup  ####
    int Nvol = CommonParameters::Nvol();
    int Ndim = CommonParameters::Ndim();

    Parameters params_all = ParameterManager::read(filename_input);

    Parameters params_test       = params_all.lookup("Test_WilsonLoop");
    Parameters params_wilsonloop = params_all.lookup("WilsonLoop");
    Parameters params_proj       = params_all.lookup("Projection");
    Parameters params_smear      = params_all.lookup("Smear_APE_spatial");

    const string        str_gconf_status = params_test.get_string("gauge_config_status");
    const string        str_gconf_read   = params_test.get_string("gauge_config_type_input");
    const string        readfile         = params_test.get_string("config_filename_input");
    const string        str_rand_type    = params_test.get_string("random_number_type");
    const unsigned long seed             = params_test.get_unsigned_long("seed_for_random_number");
    const int           Nsmear           = params_test.get_int("number_of_max_smearing");
    const int           Nmeas            = params_test.get_int("number_of_measurement_step");
    const string        str_vlevel       = params_test.get_string("verbose_level");

    const bool   do_check        = params_test.is_set("expected_result");
    const double expected_result = do_check ? params_test.get_double("expected_result") : 0.0;

    const string str_proj_type = params_proj.get_string("projection_type");

    Bridge::VerboseLevel vl = vout.set_verbose_level(str_vlevel);

    //- print input parameters
    vout.general(vl, "  gconf_status = %s\n", str_gconf_status.c_str());
    vout.general(vl, "  gconf_read   = %s\n", str_gconf_read.c_str());
    vout.general(vl, "  readfile     = %s\n", readfile.c_str());
    vout.general(vl, "  rand_type    = %s\n", str_rand_type.c_str());
    vout.general(vl, "  seed         = %lu\n", seed);
    vout.general(vl, "  Nsmear       = %d\n", Nsmear);
    vout.general(vl, "  Nmeas        = %d\n", Nmeas);
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


    // #### Set up a gauge configuration  ####
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


    // #### object setup #####
    WilsonLoop wilsonloop;
    wilsonloop.set_parameters(params_wilsonloop);

    unique_ptr<Projection> proj(Projection::New(str_proj_type));
    proj->set_parameters(params_proj);

    unique_ptr<Smear> smear(Smear::New("APE_spatial", proj));
    smear->set_parameters(params_smear);

    unique_ptr<Timer> timer(new Timer(test_name));


    // ####  Execution main part  ####
    timer->start();

    Field_G Uorg(Nvol, Ndim), Usmear(Nvol, Ndim);
    Uorg   = *U;
    Usmear = *U;

    double result = 0.0;
    for (int i_smear = 0; i_smear <= Nsmear; ++i_smear) {
      vout.general(vl, "i_smear = %d\n", i_smear);

      if (i_smear > 0) smear->smear(Usmear, Uorg);

      if ((i_smear % Nmeas) == 0) {
        result = wilsonloop.measure(Usmear);
      }

      Uorg = Usmear;
    }

    timer->report();

    RandomNumberManager::finalize();


    if (do_check) {
      return Test::verify(result, expected_result);
    } else {
      vout.detailed(vl, "check skipped: expected_result not set.\n\n");
      return EXIT_SKIP;
    }
  }
} // namespace Test_WilsonLoop
