/*!
        @file    $Id: test_Rational_Inverse.cpp #$

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
  const std::string test_name = "Rational.Inverse";

  //- test-private parameters
  namespace {
    const std::string filename_input = "test_Rational_Inverse.yaml";
  }

  //- prototype declaration
  int inverse(void);

#ifdef USE_TESTMANAGER_AUTOREGISTER
  namespace {
#if defined(USE_GROUP_SU2)
    // Nc=2 is not available.
#else
    static const bool is_registered = TestManager::RegisterTest(
      test_name,
      inverse
      );
#endif
  }
#endif

  //====================================================================
  int inverse(void)
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
    unique_ptr<Fopr> fopr_c(Fopr::New("Clover", str_gmset_type));
    fopr_c->set_parameters(params_clover);
    fopr_c->set_config(U);

    unique_ptr<Fopr_Rational> fopr_r(new Fopr_Rational(fopr_c));
    fopr_r->set_parameters(params_rational);

    unique_ptr<Source> source(Source::New(str_source_type));
    source->set_parameters(params_source);

    unique_ptr<Timer> timer(new Timer(test_name));


    // ####  Execution main part  ####
    timer->start();

    Field   xq(NinF, Nvol, 1), b(NinF, Nvol, 1);
    Field   v(NinF, Nvol, 1), w(NinF, Nvol, 1);
    Field_F b2;

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

        fopr_r->mult(v, b);

        fopr_c->set_mode("DdagD");
        fopr_c->mult(w, v);

        fopr_r->mult(v, w);
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
