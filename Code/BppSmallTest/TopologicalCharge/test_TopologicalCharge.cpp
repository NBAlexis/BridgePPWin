/*!
        @file    $Id: test_TopologicalCharge.cpp #$

        @brief

        @author  Yusuke Namekawa  (namekawa)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate: 2013-01-22 13:51:53 #$

        @version $LastChangedRevision: 1571 $
*/

#include "BppSmallTest.h"

//====================================================================
//! Test of Topological Charge measurement.

/*!
    This class measures the topological charge
    defined by the clover leaf on the lattice.
                          [01 Jan 2014 Y.Namekawa]
    unique_ptr is introduced to avoid memory leaks
                          [21 Mar 2015 Y.Namekawa]
 */

namespace Test_TopologicalCharge {
  const std::string test_name = "TopologicalCharge";

  //- test-private parameters
  namespace {
    const std::string filename_input = "test_TopologicalCharge.yaml";
  }

  //- prototype declaration
  int topological_charge(void);

#ifdef USE_TESTMANAGER_AUTOREGISTER
  namespace {
#if defined(USE_GROUP_SU2)
    // Nc=2 is not available.
#else
    static const bool is_registered = TestManager::RegisterTest(
      test_name,
      topological_charge
      );
#endif
  }
#endif

  //====================================================================
  int topological_charge(void)
  {
    // ####  parameter setup  ####
    int Nvol = CommonParameters::Nvol();
    int Ndim = CommonParameters::Ndim();

    Parameters params_all = ParameterManager::read(filename_input);

    Parameters params_test  = params_all.lookup("Test_TopologicalCharge");
    Parameters params_topo  = params_all.lookup("TopologicalCharge");
    Parameters params_proj  = params_all.lookup("Projection");
    Parameters params_smear = params_all.lookup("Smear");

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

    const string str_proj_type  = params_proj.get_string("projection_type");
    const string str_smear_type = params_smear.get_string("smear_type");

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
    vout.general(vl, "  proj_type    = %s\n", str_proj_type.c_str());
    vout.general(vl, "  smear_type   = %s\n", str_smear_type.c_str());
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
    unique_ptr<TopologicalCharge> topological_charge(new TopologicalCharge());
    topological_charge->set_parameters(params_topo);

    unique_ptr<Projection> proj(Projection::New(str_proj_type));
    proj->set_parameters(params_proj);

    unique_ptr<Smear> smear(Smear::New(str_smear_type, proj));
    smear->set_parameters(params_smear);

    unique_ptr<Timer> timer(new Timer(test_name));


    // ####  Execution main part  ####
    timer->start();

    Field_G Uorg(Nvol, Ndim), Usmear(Nvol, Ndim);
    Uorg   = *U;
    Usmear = *U;

    double result = 0.0;
    for (int i_smear = 0; i_smear <= Nsmear; ++i_smear) {
      if (i_smear > 0) smear->smear(Usmear, Uorg);

      if ((i_smear % Nmeas) == 0) {
        result = topological_charge->measure(Usmear);

        vout.general(vl, "i_smear,Q = %d %20.16e\n", i_smear, result);
        vout.general(vl, "\n");
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
} // namespace Test_TopologicalCharge
