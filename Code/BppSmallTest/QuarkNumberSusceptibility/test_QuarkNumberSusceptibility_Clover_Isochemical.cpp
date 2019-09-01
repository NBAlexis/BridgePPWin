/*!
        @file    test_QuarkNumberSusceptibility_Clover_Isochemical.cpp

        @brief

        @author  Hideo Matsufuru  (matsufuru)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate: 2013-01-22 13:51:53 #$

        @version $LastChangedRevision: 1929 $
*/
#include "BppSmallTest.h"
#include "test.h"

#include "Fopr/fopr_Clover_Isochemical.h"
#include "Fopr/fopr_Smeared.h"

#include "IO/gaugeConfig.h"

#include "Measurements/Fermion/fprop_Standard_lex.h"
#include "Measurements/Fermion/noiseVector_Z2.h"
#include "Measurements/Fermion/quarkNumberSusceptibility_Wilson.h"

#include "Tools/randomNumberManager.h"
#include "Tools/randomNumbers_Mseries.h"

//====================================================================
//! Quark number susceptibility for the Wilson-element_type fermion.

/*!
    This class tests measurements of the traces which is used to determine
    the quark number susceptibility.
                                [02 Sep 2012 H.Matsufuru]
    (Coding history will be recovered from trac.)
    YAML is implemented.        [14 Nov 2012 Y.Namekawa]
    Selectors are implemented.  [03 Mar 2013 Y.Namekawa]
    (Selectors are replaced with factories by Aoyama-san)
    unique_ptr is introduced to avoid memory leaks
                                [21 Mar 2015 Y.Namekawa]
 */

namespace Test_QuarkNumSuscept {
  const std::string test_name = "QuarkNumberSusceptibility.Clover_Isochemical";

  // test-private parameters
  namespace {
    const std::string filename_input = "test_QuarkNumberSusceptibility_Clover_Isochemical.yaml";
  }

  // prototype declaration
  int quark_num_suscept(void);

#ifdef USE_TESTMANAGER_AUTOREGISTER
  namespace {
#if defined(USE_GROUP_SU2)
    // Nc=2 is not available.
#else
    static const bool is_registered = TestManager::RegisterTest(
      test_name,
      quark_num_suscept
      );
#endif
  }
#endif

  //====================================================================
  int quark_num_suscept(void)
  {
    // #####  parameter setup  #####
    const int Nc   = CommonParameters::Nc();
    const int Nvol = CommonParameters::Nvol();
    const int Ndim = CommonParameters::Ndim();

    const Parameters params_all = ParameterManager::read(filename_input);

    const Parameters params_test     = params_all.lookup("Test_QuarkNumSuscept_Clover_Isochemical");
    const Parameters params_clover   = params_all.lookup("Fopr_Clover_Isochemical");
    const Parameters params_proj     = params_all.lookup("Projection");
    const Parameters params_smear    = params_all.lookup("Smear");
    const Parameters params_dr_smear = params_all.lookup("Director_Smear");
    const Parameters params_solver   = params_all.lookup("Solver");

    const string        str_gconf_status = params_test.get_string("gauge_config_status");
    const string        str_gconf_read   = params_test.get_string("gauge_config_type_input");
    const string        readfile         = params_test.get_string("config_filename_input");
    const string        str_rand_type    = params_test.get_string("random_number_type");
    const unsigned long seed             = params_test.get_unsigned_long("seed_for_random_number");
    const int           i_seed_noise     = params_test.get_int("int_seed_for_noise");
    const string        str_vlevel       = params_test.get_string("verbose_level");

    const bool   do_check        = params_test.is_set("expected_result");
    const double expected_result = do_check ? params_test.get_double("expected_result") : 0.0;

    const string str_gmset_type  = params_clover.get_string("gamma_matrix_type");
    const string str_proj_type   = params_proj.get_string("projection_type");
    const string str_smear_type  = params_smear.get_string("smear_type");
    const string str_solver_type = params_solver.get_string("solver_type");

    const Bridge::VerboseLevel vl = vout.set_verbose_level(str_vlevel);

    //- print input parameters
    vout.general(vl, "  gconf_status = %s\n", str_gconf_status.c_str());
    vout.general(vl, "  gconf_read   = %s\n", str_gconf_read.c_str());
    vout.general(vl, "  readfile     = %s\n", readfile.c_str());
    vout.general(vl, "  rand_type    = %s\n", str_rand_type.c_str());
    vout.general(vl, "  seed         = %lu\n", seed);
    vout.general(vl, "  i_seed_noise = %d\n", i_seed_noise);
    vout.general(vl, "  vlevel       = %s\n", str_vlevel.c_str());
    vout.general(vl, "  gmset_type   = %s\n", str_gmset_type.c_str());
    vout.general(vl, "  proj_type    = %s\n", str_proj_type.c_str());
    vout.general(vl, "  smear_type   = %s\n", str_smear_type.c_str());
    vout.general(vl, "  solver_type  = %s\n", str_solver_type.c_str());

    //- input parameter check
    int err = 0;
    err += ParameterCheck::non_NULL(str_gconf_status);
    err += ParameterCheck::non_zero(i_seed_noise);

    if (err) {
      vout.crucial(vl, "Error at %s: input parameters have not been set\n", test_name.c_str());
      exit(EXIT_FAILURE);
    }


    RandomNumberManager::initialize(str_rand_type, seed);


    // #####  Set up a gauge configuration  ####
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
    unique_ptr<Fopr> fopr_c(Fopr::New("Clover_Isochemical", str_gmset_type));
    fopr_c->set_parameters(params_clover);

    unique_ptr<Projection> proj(Projection::New(str_proj_type));
    proj->set_parameters(params_proj);

    unique_ptr<Smear> smear(Smear::New(str_smear_type, proj));
    smear->set_parameters(params_smear);

    unique_ptr<Director> dr_smear(new Director_Smear(smear));
    dr_smear->set_parameters(params_dr_smear);

    unique_ptr<Fopr> fopr_smear(Fopr::New("Smeared", fopr_c, dr_smear));
    fopr_smear->set_config(U);


    //- Random number is initialized with a parameter specified by iseed
    unique_ptr<RandomNumbers> rand(new RandomNumbers_Mseries(i_seed_noise));
    unique_ptr<NoiseVector>   nv(new NoiseVector_Z2(rand));

    unique_ptr<Solver> solver(Solver::New(str_solver_type, fopr_smear));
    solver->set_parameters(params_solver);
    unique_ptr<Fprop> fprop_lex(new Fprop_Standard_lex(solver));

    const unique_ptr<QuarkNumberSusceptibility_Wilson> quark_suscept(new QuarkNumberSusceptibility_Wilson(fopr_smear, fprop_lex, nv));
    quark_suscept->set_parameters(params_test);

    const unique_ptr<Timer> timer(new Timer(test_name));


    // ####  Execution main part  ####
    timer->start();

    const double result = quark_suscept->measure();

    timer->report();


    RandomNumberManager::finalize();


    if (do_check) {
      return Test::verify(result, expected_result);
    } else {
      vout.detailed(vl, "check skipped: expected_result not set.\n\n");
      return EXIT_SKIP;
    }
  }
} // namespace Test_QuarkNumSuscept
