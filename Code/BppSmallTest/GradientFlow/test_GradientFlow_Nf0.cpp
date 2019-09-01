/*!
        @file    test_GradientFlow_Nf0.cpp

        @brief

        @author  Sinya Aoki (saoki)

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/
#include "BppSmallTest.h"
#include "test.h"

#include "IO/gaugeConfig.h"

#include "Measurements/Gauge/energyDensity.h"
#include "Measurements/Gauge/gradientFlow.h"

#include "Tools/randomNumberManager.h"

//====================================================================
//! Test of gradientFlow.

/*!
                               [08 Jul 2012 S.Aoki]
    (Coding history will be recovered from trac.)
    YAML is implemented.       [14 Nov 2012 Y.Namekawa]
    4th order Runge-Kutta in commutator-free method
    formulated by E.Celledoni et al. FGCS 19, 341 (2003),
    as well as 1st and 2nd order Runge-Kutta are implemented.
                               [10 Oct 2014 Y.Namekawa]
    unique_ptr is introduced to avoid memory leaks.
                               [21 Mar 2015 Y.Namekawa]
    Adaptive stepsize control is implemented.
                               [01 May 2015 Y.Namekawa]
 */

namespace Test_GradientFlow {
  const std::string test_name = "GradientFlow.Nf0";

  const std::string test_name_RK1         = "GradientFlow.Nf0.RK1";
  const std::string test_name_RK2         = "GradientFlow.Nf0.RK2";
  const std::string test_name_RK3         = "GradientFlow.Nf0.RK3";
  const std::string test_name_RK4         = "GradientFlow.Nf0.RK4";
  const std::string test_name_RK_adaptive = "GradientFlow.Nf0.RK_adaptive";

  //- test-private parameters
  namespace {
    const std::string filename_input_RK1         = "test_GradientFlow_Nf0_RK1.yaml";
    const std::string filename_input_RK2         = "test_GradientFlow_Nf0_RK2.yaml";
    const std::string filename_input_RK3         = "test_GradientFlow_Nf0_RK3.yaml";
    const std::string filename_input_RK4         = "test_GradientFlow_Nf0_RK4.yaml";
    const std::string filename_input_RK_adaptive = "test_GradientFlow_Nf0_RK_adaptive.yaml";
  }

  //- prototype declaration
  int update(const std::string& filename_input);

  int run_test_RK1()
  { return update(filename_input_RK1); }

  int run_test_RK2()
  { return update(filename_input_RK2); }

  int run_test_RK3()
  { return update(filename_input_RK3); }

  int run_test_RK4()
  { return update(filename_input_RK4); }

  int run_test_RK_adaptive()
  { return update(filename_input_RK_adaptive); }

#ifdef USE_TESTMANAGER_AUTOREGISTER
  namespace {
#if defined(USE_GROUP_SU2)
    // Nc=2 is not available.
#else
    static const bool is_registered_RK1 = TestManager::RegisterTest(
      test_name_RK1,
      run_test_RK1
      );
    static const bool is_registered_RK2 = TestManager::RegisterTest(
      test_name_RK2,
      run_test_RK2
      );
    static const bool is_registered_RK3 = TestManager::RegisterTest(
      test_name_RK3,
      run_test_RK3
      );
    static const bool is_registered_RK4 = TestManager::RegisterTest(
      test_name_RK4,
      run_test_RK4
      );
    static const bool is_registered_RK_adaptive = TestManager::RegisterTest(
      test_name_RK_adaptive,
      run_test_RK_adaptive
      );
#endif
  }
#endif

  //====================================================================
  int update(const std::string& filename_input)
  {
    // ####  parameter setup  ####
    const int Nvol = CommonParameters::Nvol();
    const int Ndim = CommonParameters::Ndim();
    const int Nc   = CommonParameters::Nc();

    const Parameters params_all = ParameterManager::read(filename_input);

    const Parameters params_test = params_all.lookup("Test_GradientFlow");
    //- NB. additional parameter for action_G is set in the following object setup
    Parameters       params_action_G       = params_all.lookup("Action_G");
    const Parameters params_g_flow         = params_all.lookup("GradientFlow");
    const Parameters params_energy_density = params_all.lookup("EnergyDensity");

    const string        str_gconf_status = params_test.get_string("gauge_config_status");
    const string        str_gconf_read   = params_test.get_string("gauge_config_type_input");
    const string        readfile         = params_test.get_string("config_filename_input");
    const string        str_gconf_write  = params_test.get_string("gauge_config_type_output");
    const string        writefile        = params_test.get_string("config_filename_output");
    const string        str_rand_type    = params_test.get_string("random_number_type");
    const unsigned long seed             = params_test.get_unsigned_long("seed_for_random_number");
    const int           Nstep            = params_test.get_int("number_of_steps");
    const double        t_flow_max       = params_test.get_double("max_flow_time");
    const string        str_vlevel       = params_test.get_string("verbose_level");

    const bool   do_check        = params_test.is_set("expected_result");
    const double expected_result = do_check ? params_test.get_double("expected_result") : 0.0;

    const string str_action_G_type = params_action_G.get_string("action_type");

    const Bridge::VerboseLevel vl = vout.set_verbose_level(str_vlevel);

    //- print input parameters
    vout.general(vl, "  gconf_status = %s\n", str_gconf_status.c_str());
    vout.general(vl, "  gconf_read   = %s\n", str_gconf_read.c_str());
    vout.general(vl, "  readfile     = %s\n", readfile.c_str());
    vout.general(vl, "  gconf_write  = %s\n", str_gconf_write.c_str());
    vout.general(vl, "  writefile    = %s\n", writefile.c_str());
    vout.general(vl, "  rand_type    = %s\n", str_rand_type.c_str());
    vout.general(vl, "  seed         = %lu\n", seed);
    vout.general(vl, "  Nstep        = %d\n", Nstep);
    vout.general(vl, "  t_flow_max   = %10.6f\n", t_flow_max);
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

    unique_ptr<Action> action_G(Action::New(str_action_G_type));
    //- beta is overwritten to be Nc in GradientFlow
    params_action_G.set_double("beta", static_cast<double>(Nc));
    action_G->set_parameters(params_action_G);

    const unique_ptr<GradientFlow> g_flow(new GradientFlow(action_G));
    g_flow->set_parameters(params_g_flow);

    const unique_ptr<EnergyDensity> energy_density(new EnergyDensity);
    energy_density->set_parameters(params_energy_density);

    const unique_ptr<Timer> timer(new Timer(test_name));


    // ####  Execution main part  ####
    timer->start();

    double t_flow = 0.0;
    double result = 0.0;

    for (int i = 0; i < Nstep; ++i) {
      result = g_flow->evolve(t_flow, *U);

      double t2       = t_flow * t_flow;
      double E_plaq   = energy_density->E_plaq(*U);
      double E_clover = energy_density->E_clover(*U);

      vout.general(vl, "  (t, t^2 E_plaq, t^2 E_clover) = %.8f %.16f %.16f\n",
                   t_flow, t2 * E_plaq, t2 * E_clover);


      if (t_flow > t_flow_max) break;
    }

    timer->report();

    RandomNumberManager::finalize();


    if (do_check) {
      //- NB. verification criterion is loosed for RK_adaptive,
      //      due to its severe sensitivity to rounding errors
      if (filename_input == "test_GradientFlow_Nf0_RK_adaptive.yaml") {
        return Test::verify(result, expected_result, 1.0e-10);
      } else {
        return Test::verify(result, expected_result);
      }
    } else {
      vout.detailed(vl, "check skipped: expected_result not set.\n\n");
      return EXIT_SKIP;
    }
  }
} // namespace Test_GradientFlow
