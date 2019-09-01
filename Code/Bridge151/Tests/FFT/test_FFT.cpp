/*!
        @file    $Id: test_FFT.cpp #$

        @brief

        @author  Yusuke Namekawa (namekawa)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate: 2013-01-22 13:51:53 #$

        @version $LastChangedRevision: 1593 $
*/

#include "Tests/test.h"

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
    int Lt = CommonParameters::Lt();

    //- local size
    int Nx   = CommonParameters::Nx();
    int Ny   = CommonParameters::Ny();
    int Nz   = CommonParameters::Nz();
    int Nvol = CommonParameters::Nvol();

    int Nc = CommonParameters::Nc();

    Parameters params_all = ParameterManager::read(filename_input);

    Parameters params_test   = params_all.lookup("Test_FFT");
    Parameters params_fft    = params_all.lookup("FFT");
    Parameters params_source = params_all.lookup("Source");

    const string str_vlevel = params_test.get_string("verbose_level");

    const bool   do_check        = params_test.is_set("expected_result");
    const double expected_result = do_check ? params_test.get_double("expected_result") : 0.0;

    const string str_fft_type    = params_fft.get_string("FFT_type");
    const string str_source_type = params_source.get_string("source_type");

    Bridge::VerboseLevel vl = vout.set_verbose_level(str_vlevel);


    // #### object setup #####
    unique_ptr<FFT> fft(FFT::New(str_fft_type));
    fft->set_parameters(params_fft);

    unique_ptr<Source> source(Source::New(str_source_type));
    source->set_parameters(params_source);

    unique_ptr<Index_lex> index(new Index_lex);
    unique_ptr<Timer>     timer(new Timer(test_name));


    // ####  Execution main part  ####
    Field_F b;
    b.set(0.0);

    int i_spin  = 0;
    int i_color = 0;

    int idx = i_color + Nc * i_spin;
    source->set(b, idx);
    vout.general(vl, "b.norm2(before FFT) = %e\n", b.norm2());


    timer->start();

    fft->fft(b);
    vout.general(vl, "b.norm2(after FFT) = %e\n", b.norm2());

    double result = 0.0;
    if (Communicator::nodeid() == 0) {
      int i_site = index->site(0, 0, 0, 0);
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
