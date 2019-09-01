/*!
        @file    test_Solver_Wilson.cpp

        @brief

        @author  Yusuke Namekawa (namekawa)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/
#include "BppSmallTest.h"
#include "test.h"

#include "IO/gaugeConfig.h"

#include "Measurements/Fermion/fprop_Standard_lex.h"
#include "Measurements/Fermion/source.h"
#include "Measurements/Gauge/gaugeFixing.h"

#include "Tools/gammaMatrixSet.h"
#include "Tools/randomNumberManager.h"

#ifdef LIB_CPP11
using std::round;
#endif

//====================================================================
//! Test of Solver with Wilson fermion

/*!
    This class executes solvers, mainly as a performance monitor.
                          [16 May 2018 Y.Namekawa]
 */

namespace Test_Solver_Wilson {
  const std::string test_name = "Solver.Wilson";

  //- test-private parameters
  namespace {
    // const std::string filename_input = "test_Solver_Wilson.yaml";
  }

  //- prototype declaration
  int solver(const std::string&);

  //- solver for various algorithms
  int solver_BiCGStab_Cmplx()
  {
    return solver("test_Solver_Wilson_BiCGStab_Cmplx.yaml");
  }


  int solver_BiCGStab_DS_L_Cmplx()
  {
    return solver("test_Solver_Wilson_BiCGStab_DS_L_Cmplx.yaml");
  }


  int solver_BiCGStab_IDS_L_Cmplx()
  {
    return solver("test_Solver_Wilson_BiCGStab_IDS_L_Cmplx.yaml");
  }


  int solver_BiCGStab_L_Cmplx()
  {
    return solver("test_Solver_Wilson_BiCGStab_L_Cmplx.yaml");
  }


  int solver_CGNE()
  {
    return solver("test_Solver_Wilson_CGNE.yaml");
  }


  int solver_CGNR()
  {
    return solver("test_Solver_Wilson_CGNR.yaml");
  }


  int solver_GMRES_m_Cmplx()
  {
    return solver("test_Solver_Wilson_GMRES_m_Cmplx.yaml");
  }


#ifdef USE_TESTMANAGER_AUTOREGISTER
  namespace {
#if defined(USE_GROUP_SU2)
    // Nc=2 is not available.
#else
    static const bool is_registered_BiCGStab_Cmplx = TestManager::RegisterTest(
      "Solver.Wilson.BiCGStab_Cmplx",
      solver_BiCGStab_Cmplx
      );

    static const bool is_registered_BiCGStab_DS_L_Cmplx = TestManager::RegisterTest(
      "Solver.Wilson.BiCGStab_DS_L_Cmplx",
      solver_BiCGStab_DS_L_Cmplx
      );

    static const bool is_registered_BiCGStab_IDS_L_Cmplx = TestManager::RegisterTest(
      "Solver.Wilson.BiCGStab_IDS_L_Cmplx",
      solver_BiCGStab_IDS_L_Cmplx
      );

    static const bool is_registered_BiCGStab_L_Cmplx = TestManager::RegisterTest(
      "Solver.Wilson.BiCGStab_L_Cmplx",
      solver_BiCGStab_L_Cmplx
      );

    static const bool is_registered_CGNE = TestManager::RegisterTest(
      "Solver.Wilson.CGNE",
      solver_CGNE
      );

    static const bool is_registered_CGNR = TestManager::RegisterTest(
      "Solver.Wilson.CGNR",
      solver_CGNR
      );

    static const bool is_registered_GMRES_m_Cmplx = TestManager::RegisterTest(
      "Solver.Wilson.GMRES_m_Cmplx",
      solver_GMRES_m_Cmplx
      );
#endif
  }
#endif

  //====================================================================
  int solver(const std::string& filename_input)
  {
    // ####  parameter setup  ####
    const int Nc   = CommonParameters::Nc();
    const int Nd   = CommonParameters::Nd();
    const int Ndim = CommonParameters::Ndim();
    const int Nvol = CommonParameters::Nvol();

    const Parameters params_all = ParameterManager::read(filename_input);

    const Parameters params_test   = params_all.lookup("Test_Solver");
    const Parameters params_gfix   = params_all.lookup("GaugeFixing");
    const Parameters params_fopr   = params_all.lookup("Fopr");
    const Parameters params_solver = params_all.lookup("Solver");
    const Parameters params_source = params_all.lookup("Source");

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

    const Bridge::VerboseLevel vl = vout.set_verbose_level(str_vlevel);

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
    {
      unique_ptr<Field_G>           Ufix(new Field_G(Nvol, Ndim));
      const unique_ptr<GaugeFixing> gfix(GaugeFixing::New(str_gfix_type));
      gfix->set_parameters(params_gfix);

      gfix->fix(*Ufix, *U);

      copy(*U, *Ufix);
    }


    // ####  object setup  #####
    unique_ptr<GammaMatrixSet> gmset(GammaMatrixSet::New(str_gmset_type));

    unique_ptr<Fopr> fopr(Fopr::New(str_fopr_type, str_gmset_type));
    fopr->set_parameters(params_fopr);
    fopr->set_config(U);

    unique_ptr<Solver> solver(Solver::New(str_solver_type, fopr));
    solver->set_parameters(params_solver);

    const unique_ptr<Fprop> fprop_lex(new Fprop_Standard_lex(solver));

    const unique_ptr<Source> source(Source::New(str_source_type));
    source->set_parameters(params_source);

    const unique_ptr<Timer> timer(new Timer(test_name));


    // ####  Execution main part  ####
    timer->start();

    std::vector<Field_F> sq(Nc * Nd);
    for (int i_cd = 0; i_cd < Nc * Nd; ++i_cd) {
      sq[i_cd].set(0.0);
    }

    vout.general(vl, "\n");
    vout.general(vl, "Solving quark propagator:\n");
    vout.general(vl, "  color spin   Nconv      diff           diff2\n");

    int ispin = 0;
    // for (int ispin = 0; ispin < Nd; ++ispin) {

    int icolor = 0;
    // for (int icolor = 0; icolor < Nc; ++icolor) {
    int i_cd = icolor + Nc * ispin;

    Field_F b;      // b.set(0.0);
    source->set(b, i_cd);

    int    Nconv;
    double diff;
    fprop_lex->invert_D(sq[i_cd], b, Nconv, diff);

    Field_F y(b);
    fopr->set_mode("D");
    fopr->mult(y, sq[i_cd]);     // y  = fopr->mult(sq[i_cd]);
    axpy(y, -1.0, b);            // y -= b;
    double diff2 = y.norm2() / b.norm2();

    vout.general(vl, "   %2d   %2d   %6d   %12.4e   %12.4e\n",
                 icolor, ispin, Nconv, diff, diff2);
    //   }
    // }

    const double result = round(Nconv / 10.0);

    timer->report();

    RandomNumberManager::finalize();


    if (do_check) {
      return Test::verify(result, expected_result);
    } else {
      vout.detailed(vl, "check skipped: expected_result not set.\n\n");
      return EXIT_SKIP;
    }
  }
} // namespace Test_Solver_Wilson