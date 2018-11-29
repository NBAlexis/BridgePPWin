/*!
        @file    $Id:: test_Spectrum_Wilson_2ptFunction.cpp #$

        @brief

        @author  UEDA, Satoru  (sueda)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 18:35:38 #$

        @version $LastChangedRevision: 1571 $
*/

#include "BppSmallTest.h"

//====================================================================
//! Test of Spectrum with Wilson fermion, prepared for beginners.

/*!
                                          [28 Dec 2011 H.Matsufuru]
    (Coding history will be recovered from trac.)
    YAML is implemented.                  [14 Nov 2012 Y.Namekawa]
    Selectors are implemented.            [03 Mar 2013 Y.Namekawa]
    (Selectors are replaced with factories by Aoyama-san)
    unique_ptr is introduced to avoid memory leaks
                                          [21 Mar 2015 Y.Namekawa]
 */

namespace Test_Spectrum_Wilson {
  const std::string test_name = "Spectrum.Wilson.Hadron2ptFunction";

  //- test-private parameters
  namespace {
    const std::string filename_input = "test_Spectrum_Wilson_Hadron2ptFunction.yaml";
  }

  //- prototype declaration
  int hadron_2ptFunction(void);

#ifdef USE_TESTMANAGER_AUTOREGISTER
  namespace {
#if defined(USE_GROUP_SU2)
    // Nc=2 is not available.
#else
    static const bool is_registered = TestManager::RegisterTest(
      test_name,
      hadron_2ptFunction
      );
#endif
  }
#endif

  //====================================================================
  int hadron_2ptFunction(void)
  {
    // ####  parameter setup  ####
    int Nc   = CommonParameters::Nc();
    int Nd   = CommonParameters::Nd();
    int Ndim = CommonParameters::Ndim();
    int Nvol = CommonParameters::Nvol();

    Parameters params_all = ParameterManager::read(filename_input);

    Parameters params_test   = params_all.lookup("Test_Spectrum");
    Parameters params_gfix   = params_all.lookup("GaugeFixing");
    Parameters params_fopr   = params_all.lookup("Fopr");
    Parameters params_solver = params_all.lookup("Solver");
    Parameters params_source = params_all.lookup("Source");

    const string        str_gconf_status = params_test.get_string("gauge_config_status");
    const string        str_gconf_read   = params_test.get_string("gauge_config_type_input");
    const string        readfile         = params_test.get_string("config_filename_input");
    const string        str_rand_type    = params_test.get_string("random_number_type");
    const unsigned long seed             = params_test.get_unsigned_long("seed_for_random_number");
    const string        str_vlevel       = params_test.get_string("verbose_level");

    const bool   do_check        = params_test.is_set("expected_result");
    const double expected_result = do_check ? params_test.get_double("expected_result") : 0.0;

    const string str_gfix_type   = params_gfix.get_string("gauge_fixing_type");
    const string str_fopr_type   = params_fopr.get_string("fermion_type");
    const string str_gmset_type  = params_fopr.get_string("gamma_matrix_type");
    const string str_solver_type = params_solver.get_string("solver_type");
    const string str_source_type = params_source.get_string("source_type");

    Bridge::VerboseLevel vl = vout.set_verbose_level(str_vlevel);

    //- print input parameters
    vout.general(vl, "  gconf_status = %s\n", str_gconf_status.c_str());
    vout.general(vl, "  gconf_read   = %s\n", str_gconf_read.c_str());
    vout.general(vl, "  readfile     = %s\n", readfile.c_str());
    vout.general(vl, "  rand_type    = %s\n", str_rand_type.c_str());
    vout.general(vl, "  seed         = %lu\n", seed);
    vout.general(vl, "  vlevel       = %s\n", str_vlevel.c_str());
    vout.general(vl, "  gfix_type    = %s\n", str_gfix_type.c_str());
    vout.general(vl, "  gmset_type   = %s\n", str_gmset_type.c_str());
    vout.general(vl, "  solver_type  = %s\n", str_solver_type.c_str());
    vout.general(vl, "  source_type  = %s\n", str_source_type.c_str());

    //- input parameter check
    int err = 0;
    err += ParameterCheck::non_NULL(str_gconf_status);

    if (err) {
      vout.crucial(vl, "Error at %s: input parameters have not been set\n", test_name.c_str());
      exit(EXIT_FAILURE);
    }

    if (str_solver_type == "CG") {
      vout.crucial(vl, "Error at %s: CG can not be adopted. Use CGNE,CGNR, instead.\n", test_name.c_str());
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


    // ####  Gauge fixing  ####
    unique_ptr<Field_G>     Ufix(new Field_G(Nvol, Ndim));
    unique_ptr<GaugeFixing> gfix(GaugeFixing::New(str_gfix_type));
    gfix->set_parameters(params_gfix);

    gfix->fix(*Ufix, *U);

    copy(*U, *Ufix);


    // ####  object setup  #####
    unique_ptr<GammaMatrixSet> gmset(GammaMatrixSet::New(str_gmset_type));

    unique_ptr<Fopr> fopr(Fopr::New(str_fopr_type, str_gmset_type));
    fopr->set_parameters(params_fopr);
    fopr->set_config(U);

    unique_ptr<Solver> solver(Solver::New(str_solver_type, fopr));
    solver->set_parameters(params_solver);

    unique_ptr<Fprop> fprop_lex(new Fprop_Standard_lex(solver));

    unique_ptr<Source> source(Source::New(str_source_type));
    source->set_parameters(params_source);

    unique_ptr<Timer> timer(new Timer(test_name));

    Corr2pt_4spinor corr(gmset);
    corr.set_parameters(params_all.lookup("Corr2pt_4spinor"));


    // ####  Execution main part  ####
    timer->start();

    std::vector<Field_F> sq(Nc * Nd);
    for (int i = 0; i < Nc * Nd; ++i) {
      sq[i].set(0.0);
    }

    Field_F b;
    b.set(0.0);

    int    Nconv;
    double diff;

    vout.general(vl, "\n");
    vout.general(vl, "Solving quark propagator:\n");
    vout.general(vl, "  color spin   Nconv      diff           diff2\n");

    for (int ispin = 0; ispin < Nd; ++ispin) {
      for (int icolor = 0; icolor < Nc; ++icolor) {
        int idx = icolor + Nc * ispin;
        source->set(b, idx);

        fprop_lex->invert_D(sq[idx], b, Nconv, diff);

        Field_F y(b);
        fopr->set_mode("D");
        fopr->mult(y, sq[idx]);  // y  = fopr->mult(sq[idx]);
        axpy(y, -1.0, b);        // y -= b;
        double diff2 = y.norm2() / b.norm2();

        vout.general(vl, "   %2d   %2d   %6d   %12.4e   %12.4e\n",
                     icolor, ispin, Nconv, diff, diff2);
      }
    }

    //- meson correlators
    vout.general(vl, "\n");
    vout.general(vl, "2-point correlator:\n");

    double result = corr.meson_all(sq, sq);

    timer->report();

    RandomNumberManager::finalize();


    if (do_check) {
      return Test::verify(result, expected_result);
    } else {
      vout.detailed(vl, "check skipped: expected_result not set.\n\n");
      return EXIT_SKIP;
    }
  }
} // namespace Test_Spectrum_Wilson
