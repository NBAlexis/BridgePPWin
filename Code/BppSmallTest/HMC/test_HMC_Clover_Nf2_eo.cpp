/*!
        @file    test_HMC_Clover_Nf2_eo.cpp

        @brief

        @author  Yusuke Namekawa (namekawa)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate: 2013-01-22 13:51:53 #$

        @version $LastChangedRevision: 1929 $
*/
#include "BppSmallTest.h"
#include "test.h"

#include "Action/Fermion/action_F_Standard_eo.h"

#include "Fopr/fopr_Clover_eo.h"
#include "Fopr/fopr_Smeared_eo.h"
#include "Fopr/fopr_Smeared.h"

#include "Force/Fermion/force_F_Clover_Nf2.h"
#include "Force/Fermion/force_F_Smeared.h"

#include "HMC/hmc_General.h"
#include "HMC/builder_Integrator.h"

#include "IO/gaugeConfig.h"

#include "Measurements/Fermion/fprop_Standard_eo.h"

#include "Tools/file_utils.h"
#include "Tools/randomNumberManager.h"
#include "Tools/randomNumbers_Mseries.h"

//====================================================================
//! Test of spectroscopy with clover fermion.

/*!
    This test class obtains the quark propagator for clover fermion
    and calculates typical hadron correlators.
    The quantum numbers of hadrons are specified with gamma matrices
    (GammaMatrix class instance) whose set is defined in a subclass
    of GammaMatrixSet class.
                                           [12 Apr 2012 H.Matsufuru]
    (Coding history will be recovered from trac.)
    YAML is implemented.                   [14 Nov 2012 Y.Namekawa]
    Selectors are implemented.             [02 Feb 2013 Y.Namekawa]
    (Selectors are replaced with factories by Aoyama-san)
    Even-odd and smearing are implemented. [03 Mar 2013 Y.Namekawa]
    unique_ptr is introduced to avoid memory leaks
                                           [21 Mar 2015 Y.Namekawa]
 */

namespace Test_HMC_Clover {
  const std::string test_name = "HMC.Clover.Nf2_eo";

  //- test-private parameters
  namespace {
    const std::string filename_input = "test_HMC_Clover_Nf2.yaml";
  }

  //- prototype declaration
  int update_Nf2_eo(void);

#ifdef USE_TESTMANAGER_AUTOREGISTER
  namespace {
#if defined(USE_GROUP_SU2)
    // Nc=2 is not available.
#else
    static const bool is_registered = TestManager::RegisterTest(
      test_name,
      update_Nf2_eo
      );
#endif
  }
#endif

  //====================================================================
  int update_Nf2_eo(void)
  {
    // #####  parameter setup  #####
    const int Nc   = CommonParameters::Nc();
    const int Nvol = CommonParameters::Nvol();
    const int Ndim = CommonParameters::Ndim();

    const Parameters params_all = ParameterManager::read(filename_input);

    const Parameters params_test       = params_all.lookup("Test_HMC_Clover");
    const Parameters params_action_G   = params_all.lookup("Action_G");
    const Parameters params_fopr       = params_all.lookup("Fopr");
    const Parameters params_proj       = params_all.lookup("Projection");
    const Parameters params_smear      = params_all.lookup("Smear");
    const Parameters params_dr_smear   = params_all.lookup("Director_Smear");
    const Parameters params_solver_MD  = params_all.lookup("Solver_MD");
    const Parameters params_solver_H   = params_all.lookup("Solver_H");
    const Parameters params_integrator = params_all.lookup("Builder_Integrator");
    const Parameters params_hmc        = params_all.lookup("HMC_General");

    const string        str_gconf_status = params_test.get_string("gauge_config_status");
    const string        str_gconf_read   = params_test.get_string("gauge_config_type_input");
    const string        readfile         = params_test.get_string("config_filename_input");
    const string        str_gconf_write  = params_test.get_string("gauge_config_type_output");
    const string        writefile        = params_test.get_string("config_filename_output");
    const string        str_rand_type    = params_test.get_string("random_number_type");
    const unsigned long seed             = params_test.get_unsigned_long("seed_for_random_number");
    int                 i_conf           = params_test.get_int("trajectory_number");
    const int           Ntraj            = params_test.get_int("trajectory_number_step");
    const int           i_save_conf      = params_test.get_int("save_config_interval");
    const string        str_vlevel       = params_test.get_string("verbose_level");

    const bool   do_check        = params_test.is_set("expected_result");
    const double expected_result = do_check ? params_test.get_double("expected_result") : 0.0;

    const string           str_action_G_type  = params_action_G.get_string("action_type");
    const string           str_fopr_type      = params_fopr.get_string("fermion_type");
    const string           str_gmset_type     = params_fopr.get_string("gamma_matrix_type");
    const string           str_proj_type      = params_proj.get_string("projection_type");
    const string           str_smear_type     = params_smear.get_string("smear_type");
    const string           str_solver_MD_type = params_solver_MD.get_string("solver_type");
    const string           str_solver_H_type  = params_solver_H.get_string("solver_type");
    const int              Nlevels            = params_integrator.get_int("number_of_levels");
    const std::vector<int> level_action       = params_integrator.get_int_vector("level_of_actions");

    const Bridge::VerboseLevel vl = vout.set_verbose_level(str_vlevel);

    //- print input parameters
    vout.general(vl, "  gconf_status = %s\n", str_gconf_status.c_str());
    vout.general(vl, "  gconf_read     = %s\n", str_gconf_read.c_str());
    vout.general(vl, "  readfile       = %s\n", readfile.c_str());
    vout.general(vl, "  gconf_write    = %s\n", str_gconf_write.c_str());
    vout.general(vl, "  writefile      = %s\n", writefile.c_str());
    vout.general(vl, "  rand_type      = %s\n", str_rand_type.c_str());
    vout.general(vl, "  seed           = %lu\n", seed);
    vout.general(vl, "  i_conf         = %d\n", i_conf);
    vout.general(vl, "  Ntraj          = %d\n", Ntraj);
    vout.general(vl, "  i_save_conf    = %d\n", i_save_conf);
    vout.general(vl, "  vlevel         = %s\n", str_vlevel.c_str());
    vout.general(vl, "  gmset_type     = %s\n", str_gmset_type.c_str());
    vout.general(vl, "  proj_type      = %s\n", str_proj_type.c_str());
    vout.general(vl, "  smear_type     = %s\n", str_smear_type.c_str());
    vout.general(vl, "  solver_MD_type = %s\n", str_solver_MD_type.c_str());
    vout.general(vl, "  solver_H_type  = %s\n", str_solver_H_type.c_str());
    vout.general(vl, "\n");

    //- input parameter check
    int err = 0;
    err += ParameterCheck::non_NULL(str_gconf_status);
    err += ParameterCheck::non_negative(i_conf);
    err += ParameterCheck::non_negative(Ntraj);
    err += ParameterCheck::non_negative(i_save_conf);

    if (err) {
      vout.crucial(vl, "Error at %s: input parameters have not been set\n", test_name.c_str());
      exit(EXIT_FAILURE);
    }

    if ((str_solver_MD_type == "CG") || (str_solver_H_type == "CG")) {
      vout.crucial(vl, "Error at %s: CG can not be adopted. Use CGNE,CGNR, instead.\n", test_name.c_str());
      exit(EXIT_FAILURE);
    }

    RandomNumberManager::initialize(str_rand_type, seed);


    // #####  object setup  #####
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

    const unique_ptr<GaugeConfig> gconf_write(new GaugeConfig(str_gconf_write));


    unique_ptr<Action> action_G(Action::New(str_action_G_type));
    action_G->set_parameters(params_action_G);

    //-- N_f=2 part
    unique_ptr<Fopr> fopr(Fopr::New(str_fopr_type, str_gmset_type));
    fopr->set_parameters(params_fopr);

    unique_ptr<Force> force_fopr(new Force_F_Clover_Nf2(str_gmset_type));
    // define fermion force (SA)
    force_fopr->set_parameters(params_fopr);

    // define smearing method (SA)
    unique_ptr<Projection> proj(Projection::New(str_proj_type));
    proj->set_parameters(params_proj);

    unique_ptr<Smear> smear(Smear::New(str_smear_type, proj));
    smear->set_parameters(params_smear);

    // define force smearing method (SA)
    unique_ptr<ForceSmear> force_smear(ForceSmear::New(str_smear_type, proj));
    force_smear->set_parameters(params_smear);

    unique_ptr<Director> dr_smear(new Director_Smear(smear));
    dr_smear->set_parameters(params_dr_smear);

    unique_ptr<Fopr> fopr_smear(Fopr::New("Smeared", fopr, dr_smear));
    // define smeared fermion operator (SA)
    unique_ptr<Force> force_fopr_smear(
      new Force_F_Smeared(force_fopr, force_smear, dr_smear));
    // define smeared fermion force (SA)

    //- NB1 fopr_eo->set_config is performed in Action_F_Standard_eo.
    //  NB2 fopr_eo->set_mode   is performed in fprop.
    unique_ptr<Fopr> fopr_eo(Fopr::New(str_fopr_type + "_eo", str_gmset_type));
    fopr_eo->set_parameters(params_fopr);

    unique_ptr<Fopr> fopr_smear_eo(Fopr::New("Smeared_eo", fopr_eo, dr_smear));

    unique_ptr<Solver> solver_eo_MD(Solver::New(str_solver_MD_type, fopr_smear_eo));
    solver_eo_MD->set_parameters(params_solver_MD);
    unique_ptr<Fprop> fprop_eo_MD(new Fprop_Standard_eo(solver_eo_MD));

    unique_ptr<Solver> solver_eo_H(Solver::New(str_solver_H_type, fopr_smear_eo));
    solver_eo_H->set_parameters(params_solver_H);
    unique_ptr<Fprop> fprop_eo_H(new Fprop_Standard_eo(solver_eo_H));

    unique_ptr<Action> action_F(
      new Action_F_Standard_eo(fopr_smear, force_fopr_smear,
                               fprop_eo_MD, fprop_eo_H));
    // define fermion action (SA)


    ActionList actions(Nlevels);
    actions.append(level_action[0], action_F);
    actions.append(level_action[1], action_G);

    std::vector<Director *> directors(1);
    directors[0] = (Director *)dr_smear.get(); // register director[0] (SA)

    const unique_ptr<Builder_Integrator> builder(new Builder_Integrator(actions, directors));
    builder->set_parameters(params_integrator);
    Integrator *integrator = builder->build();

    //- Random number is initialized with a parameter specified by i_conf
    unique_ptr<RandomNumbers> rand(new RandomNumbers_Mseries(i_conf));

    HMC_General hmc(actions, directors, integrator, rand); // define hmc_leapfrog (SA)
    hmc.set_parameters(params_hmc);

    const unique_ptr<Timer> timer(new Timer(test_name));


    // ####  Execution main part  ####
    timer->start();

    vout.general(vl, "HMC: Ntraj = %d\n", Ntraj);  // a number of trajectory (SA)

    double result = 0.0;
    for (int traj = 0; traj < Ntraj; ++traj) {
      vout.general(vl, "\n");
      vout.general(vl, "traj = %d\n", traj);

      result = hmc.update(*U); // hmc update (SA)

      if ((i_conf + traj + 1) % i_save_conf == 0) {
        std::string filename = FileUtils::generate_filename("%s-%06d", writefile.c_str(), (i_conf + traj + 1));
        gconf_write->write_file(U, filename);
      }
    }

    gconf_write->write_file(U, writefile);

    timer->report();

    RandomNumberManager::finalize();


    if (do_check) {
      return Test::verify(result, expected_result);
    } else {
      vout.detailed(vl, "check skipped: expected_result not set.\n\n");
      return EXIT_SKIP;
    }
  }
} // namespace Test_HMC_Clover
