/*!
        @file    $Id: test_Rational_Smeared.cpp #$

        @brief

        @author  Hideo Matsufuru  (matsufuru)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 18:35:38 #$

        @version $LastChangedRevision: 1571 $
*/

#include "BppSmallTest.h"

//====================================================================
//! Test of rational approximation of fermion operators.

/*!
                                [28 Dec 2011 H.Matsufuru]
    (Coding history will be recovered from trac.)
    YAML is implemented.        [14 Nov 2012 Y.Namekawa]
    Selectors are implemented.  [03 Mar 2013 Y.Namekawa]
    (Selectors are replaced with factories by Aoyama-san)
    unique_ptr is introduced to avoid memory leaks
                                [21 Mar 2015 Y.Namekawa]
 */

namespace Test_Rational {
  const std::string test_name = "Rational.Smeared_Rational";

  //- test-private parameters
  namespace {
    const std::string filename_input = "test_Rational_Smeared.yaml";
  }

  //- prototype declaration
  int smeared_rational(void);

#ifdef USE_TESTMANAGER_AUTOREGISTER
  namespace {
#if defined(USE_GROUP_SU2)
    // Nc=2 is not available.
#else
    static const bool is_registered = TestManager::RegisterTest(
      test_name,
      smeared_rational
      );
#endif
  }
#endif

  //====================================================================
  int smeared_rational(void)
  {
    // #####  parameter setup  #####
    int Ndim = CommonParameters::Ndim();
    int Nc   = CommonParameters::Nc();
    int Nd   = CommonParameters::Nd();
    int Nvol = CommonParameters::Nvol();
    int NinF = 2 * Nc * Nd;

    Parameters params_all = ParameterManager::read(filename_input);

    Parameters params_test     = params_all.lookup("Test_Rational");
    Parameters params_clover   = params_all.lookup("Fopr_Clover");
    Parameters params_proj     = params_all.lookup("Projection");
    Parameters params_smear    = params_all.lookup("Smear");
    Parameters params_dr_smear = params_all.lookup("Director_Smear");
    Parameters params_rational = params_all.lookup("Fopr_Rational");
    Parameters params_source   = params_all.lookup("Source");

    const string        str_gconf_status = params_test.get_string("gauge_config_status");
    const string        str_gconf_read   = params_test.get_string("gauge_config_type_input");
    const string        readfile         = params_test.get_string("config_filename_input");
    const string        str_rand_type    = params_test.get_string("random_number_type");
    const unsigned long seed             = params_test.get_unsigned_long("seed_for_random_number");
    const string        str_vlevel       = params_test.get_string("verbose_level");

    const bool   do_check        = params_test.is_set("expected_result");
    const double expected_result = do_check ? params_test.get_double("expected_result") : 0.0;

    const string str_gmset_type  = params_clover.get_string("gamma_matrix_type");
    const string str_proj_type   = params_proj.get_string("projection_type");
    const string str_smear_type  = params_smear.get_string("smear_type");
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
    vout.general(vl, "  proj_type    = %s\n", str_proj_type.c_str());
    vout.general(vl, "  smear_type   = %s\n", str_smear_type.c_str());
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
    unique_ptr<Fopr> fopr_c(Fopr::New("Clover", str_gmset_type));
    fopr_c->set_parameters(params_clover);

    unique_ptr<Projection> proj(Projection::New(str_proj_type));
    proj->set_parameters(params_proj);

    unique_ptr<Smear> smear(Smear::New(str_smear_type, proj));
    smear->set_parameters(params_smear);

    unique_ptr<Director> dr_smear(new Director_Smear(smear));
    dr_smear->set_parameters(params_dr_smear);

    unique_ptr<Fopr> fopr_r(Fopr::New("Rational", fopr_c));
    fopr_r->set_parameters(params_rational);

    unique_ptr<Source> source(Source::New(str_source_type));
    source->set_parameters(params_source);

    //- Smeared-rational
    Fopr_Smeared fopr_smear_1(fopr_c, dr_smear);
    Fopr_Smeared fopr_smear_2(fopr_r, dr_smear);
    fopr_smear_1.set_config(U);
    fopr_smear_2.set_config(U);

    //- Rational-smeared
    Fopr_Rational fopr_r_smear_1((Fopr *)&fopr_smear_1);
    fopr_r_smear_1.set_parameters(params_rational);
    fopr_r_smear_1.set_config(U);

    unique_ptr<Timer> timer(new Timer(test_name));


    // ####  Execution main part  ####
    timer->start();

    Field   xq(NinF, Nvol, 1), b(NinF, Nvol, 1);
    Field   v(NinF, Nvol, 1), w(NinF, Nvol, 1);
    Field_F b2;

    vout.general(vl, "Consistency check of smeared-rational and rational-smeared.\n");

    vout.general(vl, "Smeared-rational:\n");

    {
      int ispin = 0;
      {
        int icolor = 0;
        //for(int ispin = 0; ispin < Nd; ++ispin){
        //    for(int icolor = 0; icolor < Nc; ++icolor){

        int idx = icolor + Nc * ispin;
        source->set(b2, idx);

        b = (Field)b2;

        fopr_smear_2.mult(v, b);

        fopr_smear_1.set_mode("DdagD");
        fopr_smear_1.mult(w, v);

        fopr_smear_2.mult(v, w);

        axpy(v, -1.0, b);  // v -= b;
        double vv = v.norm2();

        vout.general(vl, "  standard norm2 = %.8e\n", vv);
      }
    }

    vout.general(vl, "Rational-smeared:\n");

    double vv;
    {
      int ispin = 0;
      {
        int icolor = 0;
        //for(int ispin = 0; ispin < Nd; ++ispin){
        //    for(int icolor = 0; icolor < Nc; ++icolor){

        int idx = icolor + Nc * ispin;
        source->set(b2, idx);

        b = (Field)b2;

        fopr_r_smear_1.mult(v, b);

        fopr_smear_1.set_mode("DdagD");
        fopr_smear_1.mult(w, v);

        fopr_r_smear_1.mult(v, w);

        axpy(v, -1.0, b);  // v -= b;
        vv = v.norm2();

        vout.general(vl, "  standard norm2 = %.8e\n", vv);
      }
    }

    double result = vv;

    timer->report();

    RandomNumberManager::finalize();


    if (do_check) {
      return Test::verify(result, expected_result);
    } else {
      vout.detailed(vl, "check skipped: expected_result not set.\n\n");
      return EXIT_SKIP;
    }
  }
} // namespace Test_Rational
