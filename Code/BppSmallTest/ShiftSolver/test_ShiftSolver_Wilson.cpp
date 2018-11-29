/*!
        @file    $Id: test_ShiftSolver_Wilson.cpp #$

        @brief

        @author  Hideo Matsufuru  (matsufuru)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate: 2013-04-27 12:28:50 #$

        @version $LastChangedRevision: 1571 $
*/

#include "BppSmallTest.h"

//====================================================================
//! Test of multishift solver.

/*!
                               [28 Dec 2011 H.Matsufuru]
    (Coding history will be recovered from trac.)
    YAML is implemented.       [14 Nov 2012 Y.Namekawa]
    Selectors are implemented. [03 Mar 2013 Y.Namekawa]
    (Selectors are replaced with factories by Aoyama-san)
    unique_ptr is introduced to avoid memory leaks
                               [21 Mar 2015 Y.Namekawa]
*/

namespace Test_ShiftSolver {
  const std::string test_name = "ShiftSolver.Wilson";

  //- test-private parameters
  namespace {
    const std::string filename_input = "test_ShiftSolver_Wilson.yaml";
  }

  //- prototype declaration
  int solve(void);

#ifdef USE_TESTMANAGER_AUTOREGISTER
  namespace {
#if defined(USE_GROUP_SU2)
    // Nc=2 is not available.
#else
    static const bool is_registered = TestManager::RegisterTest(
      test_name,
      solve
      );
#endif
  }
#endif

  //====================================================================
  int solve(void)
  {
    // ####  parameter setup  ####
    int Nc   = CommonParameters::Nc();
    int Nd   = CommonParameters::Nd();
    int Ndim = CommonParameters::Ndim();
    int Nvol = CommonParameters::Nvol();

    Parameters params_all = ParameterManager::read(filename_input);

    Parameters params_test   = params_all.lookup("Test_ShiftSolver");
    Parameters params_wilson = params_all.lookup("Fopr_Wilson");
    Parameters params_solver = params_all.lookup("Fprop_Shift");
    Parameters params_source = params_all.lookup("Source");

    const string        str_gconf_status = params_test.get_string("gauge_config_status");
    const string        str_gconf_read   = params_test.get_string("gauge_config_type_input");
    const string        readfile         = params_test.get_string("config_filename_input");
    const string        str_rand_type    = params_test.get_string("random_number_type");
    const unsigned long seed             = params_test.get_unsigned_long("seed_for_random_number");
    const string        str_vlevel       = params_test.get_string("verbose_level");

    const bool   do_check        = params_test.is_set("expected_result");
    const double expected_result = do_check ? params_test.get_double("expected_result") : 0.0;

    const string str_gmset_type  = params_wilson.get_string("gamma_matrix_type");
    const int    Nshift          = params_solver.get_int("number_of_shifts");
    const string str_source_type = params_source.get_string("source_type");

    Bridge::VerboseLevel vl = vout.set_verbose_level(str_vlevel);

    //- print input parameters
    vout.general(vl, "  gconf_status = %s\n", str_gconf_status.c_str());
    vout.general(vl, "  gconf_read   = %s\n", str_gconf_read.c_str());
    vout.general(vl, "  readfile     = %s\n", readfile.c_str());
    vout.general(vl, "  rand_type    = %s\n", str_rand_type.c_str());
    vout.general(vl, "  seed         = %lu\n", seed);
    vout.general(vl, "  vlevel       = %s\n", str_vlevel.c_str());
    vout.general(vl, "  gmset_type   = %s\n", str_gmset_type.c_str());
    vout.general(vl, "  Nshift       = %d\n", Nshift);
    vout.general(vl, "  source_type  = %s\n", str_source_type.c_str());
    vout.general(vl, "\n");

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


    // ####  object setup  ####
    unique_ptr<Fopr> fopr_w(Fopr::New("Wilson", str_gmset_type));
    fopr_w->set_parameters(params_wilson);
    fopr_w->set_config(U);

    unique_ptr<Source> source(Source::New(str_source_type));
    source->set_parameters(params_source);

    unique_ptr<Index_lex> index(new Index_lex);

    unique_ptr<Fprop_Wilson_Shift> fprop_shift(new Fprop_Wilson_Shift(fopr_w, index));
    fprop_shift->set_parameters(params_solver);

    unique_ptr<Timer> timer(new Timer(test_name));


    // ####  Execution main part  ####
    timer->start();

    Field_F              b;
    std::vector<Field_F> xq_shift(Nshift);

    double result;
    {
      int ispin = 0;
      {
        int icolor = 0;
        //  for(int ispin = 0; ispin < Nd; ++ispin){
        //    for(int icolor = 0; icolor < Nc; ++icolor){

        int idx = icolor + Nc * ispin;
        source->set(b, idx);
        result = fprop_shift->invert_D(&xq_shift, b);
      }
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
} // namespace Test_ShiftSolver