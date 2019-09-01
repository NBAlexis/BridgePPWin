/*!
        @file    test_RandomNumbers_Field.cpp

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/
#include "BppSmallTest.h"
#include "test.h"

#include "Tools/randomNumbers.h"

//====================================================================
//! Test of random number generator.

/*!
                                [28 Dec 2011 H.Matsufuru]
    (Coding history will be recovered from trac.)
    YAML is implemented.        [14 Nov 2012 Y.Namekawa]
    unique_ptr is introduced to avoid memory leaks.
                                [21 Mar 2015 Y.Namekawa]
    Enhanced for several types of RandomField.
                                [ 2 Feb 2017 Y.Namekawa]
 */

namespace Test_RandomNumbers {
  const std::string test_name = "RandomNumbers.Field";

  //- test-private parameters
  namespace {
    // const std::string filename_input = "test_RandomNumbers_MT19937_Field_Gaussian.yaml";
  }

  //- prototype declaration
  int rand_field(const std::string&);

  //- rand_field for various kinds of RandomNumbers
  int rand_field_MT19937_Gaussian()
  {
    return rand_field("test_RandomNumbers_MT19937_Field_Gaussian.yaml");
  }


  int rand_field_MT19937_U1()
  {
    return rand_field("test_RandomNumbers_MT19937_Field_U1.yaml");
  }


  int rand_field_MT19937_Z2()
  {
    return rand_field("test_RandomNumbers_MT19937_Field_Z2.yaml");
  }


  int rand_field_Mseries_Gaussian()
  {
    return rand_field("test_RandomNumbers_Mseries_Field_Gaussian.yaml");
  }


  int rand_field_Mseries_U1()
  {
    return rand_field("test_RandomNumbers_Mseries_Field_U1.yaml");
  }


  int rand_field_Mseries_Z2()
  {
    return rand_field("test_RandomNumbers_Mseries_Field_Z2.yaml");
  }


#ifdef USE_SFMTLIB
  int rand_field_SFMT_Gaussian()
  {
    return rand_field("test_RandomNumbers_SFMT_Field_Gaussian.yaml");
  }


  int rand_field_SFMT_U1()
  {
    return rand_field("test_RandomNumbers_SFMT_Field_U1.yaml");
  }


  int rand_field_SFMT_Z2()
  {
    return rand_field("test_RandomNumbers_SFMT_Field_Z2.yaml");
  }
#endif


#ifdef USE_TESTMANAGER_AUTOREGISTER
  namespace {
    static const bool is_registered_MT19937_Gaussian = TestManager::RegisterTest(
      "RandomNumbers.MT19937.Field.Gaussian",
      rand_field_MT19937_Gaussian
      );

    static const bool is_registered_MT19937_U1 = TestManager::RegisterTest(
      "RandomNumbers.MT19937.Field.U1",
      rand_field_MT19937_U1
      );

    static const bool is_registered_MT19937_Z2 = TestManager::RegisterTest(
      "RandomNumbers.MT19937.Field.Z2",
      rand_field_MT19937_Z2
      );


    static const bool is_registered_Mseries_Gaussian = TestManager::RegisterTest(
      "RandomNumbers.Mseries.Field.Gaussian",
      rand_field_Mseries_Gaussian
      );

    static const bool is_registered_Mseries_U1 = TestManager::RegisterTest(
      "RandomNumbers.Mseries.Field.U1",
      rand_field_Mseries_U1
      );

    static const bool is_registered_Mseries_Z2 = TestManager::RegisterTest(
      "RandomNumbers.Mseries.Field.Z2",
      rand_field_Mseries_Z2
      );

#ifdef USE_SFMTLIB
    static const bool is_registered_SFMT_Gaussian = TestManager::RegisterTest(
      "RandomNumbers.SFMT.Field.Gaussian",
      rand_field_SFMT_Gaussian
      );

    static const bool is_registered_SFMT_U1 = TestManager::RegisterTest(
      "RandomNumbers.SFMT.Field.U1",
      rand_field_SFMT_U1
      );

    static const bool is_registered_SFMT_Z2 = TestManager::RegisterTest(
      "RandomNumbers.SFMT.Field.Z2",
      rand_field_SFMT_Z2
      );
#endif
  }
#endif

  //====================================================================
  int rand_field(const std::string& filename_input)
  {
    // ####  parameter setup  ####
    const Parameters params_all  = ParameterManager::read(filename_input);
    const Parameters params_test = params_all.lookup("Test_RandomNumbers");

    const string str_rand_type  = params_test.get_string("random_number_type");
    const string str_noise_type = params_test.get_string("noise_type");
    const int    iseed          = params_test.get_int("int_seed");
    const string str_vlevel     = params_test.get_string("verbose_level");

    const bool   do_check        = params_test.is_set("expected_result");
    const double expected_result = do_check ? params_test.get_double("expected_result") : 0.0;

    const Bridge::VerboseLevel vl = vout.set_verbose_level(str_vlevel);

    //- print input parameters
    vout.general(vl, "  rand_type  = %s\n", str_rand_type.c_str());
    vout.general(vl, "  noise_type = %s\n", str_noise_type.c_str());
    vout.general(vl, "  iseed      = %d\n", iseed);
    vout.general(vl, "  vlevel     = %s\n", str_vlevel.c_str());
    vout.general(vl, "\n");


    // ####  object setup  ####
    const unique_ptr<RandomNumbers> rand(RandomNumbers::New(str_rand_type, iseed));
    const unique_ptr<Timer>         timer(new Timer(test_name));


    // ####  Execution main part  ####
    timer->start();

    //- NB. Nin must be even, due to complex of U1,Z2 fields
    const int Nin  = 30;
    const int Nvol = CommonParameters::Nvol();
    const int NPE  = CommonParameters::NPE();
    const int Nex  = 33;
    Field     v(Nin, Nvol, Nex);

    double av = 0.0;
    double vr = 0.0;

    rand->lex_global(str_noise_type, v);

    const int size = v.size();
    for (int i = 0; i < size; ++i) {
      av += v.cmp(i);
      vr += v.cmp(i) * v.cmp(i);
      // vout.general(vl, "  %10.8f\n",v.cmp(i));
    }

    double av_all = Communicator::reduce_sum(av);
    double vr_all = Communicator::reduce_sum(vr);

    av = av_all / Nvol / NPE / Nin / Nex;
    vr = vr_all / Nvol / NPE / Nin / Nex - av * av;
    vr = sqrt(vr);

    vout.general(vl, "\n");
    vout.general(vl, "Distribution (Field):\n");
    vout.general(vl, "  number of samples = %10d\n", size);
    vout.general(vl, "  average           = %10.8f\n", av);
    vout.general(vl, "  variance          = %10.8f\n", vr);
    vout.general(vl, "  variance(expect)  = %10.8f\n", 1.0 / sqrt(2.0));

    double result = vr;

    timer->report();


    if (do_check) {
      return Test::verify(result, expected_result);
    } else {
      vout.detailed(vl, "check skipped: expected_result not set.\n\n");
      return EXIT_SKIP;
    }
  }
} // namespace Test_RandomNumbers
