/*!
        @file    test_FFT.cpp

        @brief

        @author  Yusuke Namekawa (namekawa)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate: 2013-01-22 13:51:53 #$

        @version $LastChangedRevision: 1929 $
*/

#include "test.h"

#include "Tools/fft.h"

#include "Field/field_F.h"
#include "Measurements/Fermion/source.h"

//====================================================================
//! Test of Fast Fourier Transformation.

/*!
                               [06 Jun 2015 Y.Namekawa]
 */

#ifdef USE_FFTWLIB
namespace Test_FFT {
  const std::string test_name = "FFT.fft";

  //- test-private parameters
  namespace {
    const std::string filename_input = "test_FFT.yaml";
  }

  //- prototype declaration
  int fft(void);

#ifdef USE_TESTMANAGER_AUTOREGISTER
  namespace {
    static const bool is_registered = TestManager::RegisterTest(
      test_name,
      fft
      );
  }
#endif

  //====================================================================
  int fft(void)
  {
    // ####  parameter setup  ####
    //- global lattice size
    const int Lt = CommonParameters::Lt();

    //- local size
    const int Nx   = CommonParameters::Nx();
    const int Ny   = CommonParameters::Ny();
    const int Nz   = CommonParameters::Nz();
    const int Nvol = CommonParameters::Nvol();

    const int Nc = CommonParameters::Nc();

    const Parameters params_all = ParameterManager::read(filename_input);

    const Parameters params_test   = params_all.lookup("Test_FFT");
    const Parameters params_fft    = params_all.lookup("FFT");
    const Parameters params_source = params_all.lookup("Source");

    const string str_vlevel = params_test.get_string("verbose_level");

    const bool   do_check        = params_test.is_set("expected_result");
    const double expected_result = do_check ? params_test.get_double("expected_result") : 0.0;

    const string str_fft_type    = params_fft.get_string("FFT_type");
    const string str_source_type = params_source.get_string("source_type");

    const Bridge::VerboseLevel vl = vout.set_verbose_level(str_vlevel);


    // #### check conditions ####
    int NPEx = CommonParameters::NPEx();
    int NPEy = CommonParameters::NPEy();
    int NPEz = CommonParameters::NPEz();
    int NPEt = CommonParameters::NPEt();

    if (str_fft_type == "FFT_xyz_1dim") {
      if ((NPEt == 1) &&
          ((NPEx * NPEy == 1) || (NPEx * NPEz == 1) || (NPEy * NPEz == 1)))
      {
        // ok.
      } else {
        vout.general(vl, "%s: FFT_xyz_1dim supports parallelization in one direction among x, y, z. skip.\n", test_name.c_str());
        return EXIT_SKIP;
      }
    } else if (str_fft_type == "FFT_xyz_3dim") {
      if ((NPEt == 1) && (NPEx == 1) && (NPEy == 1) && (NPEz == 1))
      {
        // ok.
      } else {
        vout.general(vl, "%s: FFT_xyz_3dim supports parallelization only in z direction. skip.\n", test_name.c_str());
        return EXIT_SKIP;
      }
    }


    // #### object setup #####
    const unique_ptr<FFT> fft(FFT::New(str_fft_type));
    fft->set_parameters(params_fft);

    const unique_ptr<Source> source(Source::New(str_source_type));
    source->set_parameters(params_source);

    const unique_ptr<Index_lex> index(new Index_lex);
    const unique_ptr<Timer>     timer(new Timer(test_name));


    // ####  Execution main part  ####
    Field_F b;
    b.set(0.0);

    const int i_spin  = 0;
    const int i_color = 0;

    const int idx = i_color + Nc * i_spin;
    source->set(b, idx);
    vout.general(vl, "b.norm2(before FFT) = %e\n", b.norm2());


    timer->start();

    fft->fft(b);
    vout.general(vl, "b.norm2(after FFT) = %e\n", b.norm2());

    double result = 0.0;
    if (Communicator::nodeid() == 0) {
      const int i_site = index->site(0, 0, 0, 0);
      result = b.cmp_r(i_color, i_spin, i_site);
    }

    timer->report();


    if (do_check) {
      return Test::verify(result, expected_result);
    } else {
      vout.detailed(vl, "check skipped: expected_result not set.\n\n");
      return EXIT_SKIP;
    }
  }
} // namespace Test_FFT
#endif
