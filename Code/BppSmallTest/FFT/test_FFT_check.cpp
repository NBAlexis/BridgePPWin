/*!
        @file    test_FFT_check.cpp

        @brief

        @author  Yusuke Namekawa (namekawa)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate: 2013-01-22 13:51:53 #$

        @version $LastChangedRevision: 1929 $
*/

#include "test.h"

#include "Tools/fft_xyz_1dim.h"
#include "Tools/fft_xyz_3dim.h"
#include "Tools/fft_3d.h"

#include "Field/field_F.h"

#include "IO/fieldIO_Binary.h"
#include "IO/fieldIO_Binary_Parallel.h"
#include "IO/io_format.h"
#include "Tools/randomNumbers_MT19937.h"
#include "Tools/file_utils.h"

#include "Parameters/commonParameters.h"
#include "Communicator/communicator.h"

//====================================================================
//! Test of Fast Fourier Transformation.

/*!
    Test of Fast Fourier Transformation.
                               [06 Jun 2015 Y.Namekawa]

    Tests FFT_xyz_1dim, FFT_xyz_3dim,
    FFT_3d_local, FFT_3d_parallel1, FFT_3d_parallel3
    as well as auto-mode.
    Generate random field, perform Forward and Backward FFTs,
    and compare results with reference data.
                               [07 Nov 2018 T.Aoyama]
 */

#ifdef USE_FFTWLIB
namespace Test_FFT {
  const std::string test_name = "FFT.fft_check";

  //- test-private parameters
  namespace {
    const std::string filename_input = "test_FFT_check.yaml";
  }

  //- prototype declaration
  int fft_check(void);

#ifdef USE_TESTMANAGER_AUTOREGISTER
  namespace {
    static const bool is_registered = TestManager::RegisterTest(
      test_name,
      fft_check
      );
  }
#endif

  //====================================================================
  int run_task(unique_ptr<FFT>& fft,
               const Field& src,
               const Field& ref,
               Field& res1,
               Field& res2,
               double check_precision,
               Bridge::VerboseLevel vl)
  {
    int err = 0;

    res1.set(0.0);
    res2.set(0.0);

    // ####  Execution main part  ####

    Timer time_fft_fw("FFT_Forward");

    time_fft_fw.start();
    fft->fft(res1, src, FFT::FORWARD);
    time_fft_fw.stop();

    Timer time_fft_bw("FFT_Backward");

    time_fft_bw.start();
    fft->fft(res2, res1, FFT::BACKWARD);
    time_fft_bw.stop();

    vout.general(vl, "source norm2  = %24.20e\n", src.norm2());
    vout.general(vl, "result norm2  = %24.20e\n", res1.norm2());
    vout.general(vl, "reverse norm2 = %24.20e\n", res2.norm2());

    double diff = 0.0;

    // forward check
    axpy(res1, -1.0, ref);

    diff = res1.norm();

    vout.general(vl, "diff (result - ref)  = %24.20e\n", diff);

    if (!(diff < check_precision)) ++err;

    // forward-backward reverse check
    axpy(res2, -1.0, src);

    diff = res2.norm();
    vout.general(vl, "diff (reverse - src) = %24.20e\n", diff);

    if (!(diff < check_precision)) ++err;

    time_fft_fw.report();
    time_fft_bw.report();

    return err;
  }


  //====================================================================
  int fft_check(void)
  {
    // ####  parameter setup  ####
    //- global lattice size
    const int Lx = CommonParameters::Lx();
    const int Ly = CommonParameters::Ly();
    const int Lz = CommonParameters::Lz();
    const int Lt = CommonParameters::Lt();

    //- local size
    const int Nx = CommonParameters::Nx();
    const int Ny = CommonParameters::Ny();
    const int Nz = CommonParameters::Nz();
    const int Nt = CommonParameters::Nt();

    const int Nvol = CommonParameters::Nvol();

    //- grid size
    const int NPEx = CommonParameters::NPEx();
    const int NPEy = CommonParameters::NPEy();
    const int NPEz = CommonParameters::NPEz();
    const int NPEt = CommonParameters::NPEt();


    Parameters params_all = ParameterManager::read(filename_input);

    // Parameters: Test section
    Parameters params_test = params_all.lookup("Test_FFT");

    int nex = params_test.get_int("external_dof");

    string str_vlevel      = params_test.get_string("verbose_level");
    bool   do_check        = params_test.is_set("expected_result");
    double expected_result = do_check ? params_test.get_double("expected_result") : 0.0;
    double check_precision = params_test.is_set("check_precision") ? params_test.get_double("check_precision") : 1.0e-8;

    // Parameters: FFT section
    Parameters params_fft = params_all.lookup("FFT");
    // string str_fft_type    = params_fft.get_string("FFT_type");

    Bridge::VerboseLevel vl = vout.set_verbose_level(str_vlevel);

    // #### object setup #####

    std::string filename_source = FileUtils::generate_filename("fft_src_%02dx%02dx%02dx%02d_%02d.dat", Lx, Ly, Lz, Lt, nex);

    std::string filename_reference = FileUtils::generate_filename("fft_ref_%02dx%02dx%02dx%02d_%02d.dat", Lx, Ly, Lz, Lt, nex);

    // vout.general(vl, "fft mode            = %s\n", str_fft_type.c_str());
    vout.general(vl, "external dof        = %d\n", nex);
    vout.general(vl, "source data file    = %s\n", filename_source.c_str());
    vout.general(vl, "reference data file = %s\n", filename_reference.c_str());
    vout.general(vl, "check precision     = %8.6e\n", check_precision);


#if 0
    unique_ptr<FieldIO> fio;
    if (Communicator::size() > 1) {
      fio.reset(new FieldIO_Binary_Parallel(IO_Format::Trivial));
    } else {
      fio.reset(new FieldIO_Binary(IO_Format::Trivial));
    }
#else
    unique_ptr<FieldIO> fio(new FieldIO_Binary(IO_Format::Trivial));
#endif


    // Field_F src, ref, res1, res2;
    int Nex = nex;

    Field_F src(Nvol, Nex);
    Field_F ref(Nvol, Nex);
    Field_F res1(Nvol, Nex);
    Field_F res2(Nvol, Nex);


    // find source or generate one randomly.

    std::ifstream fsrc(filename_source.c_str());
    if (fsrc) {
      vout.general(vl, "%s: source file exists.\n", test_name.c_str());

      fsrc.close();

      fio->read_file(&src, filename_source);
    } else {
      vout.general(vl, "%s: source file not exist. generate one.\n", test_name.c_str());

      RandomNumbers_MT19937 rng(0xdeadbeefUL);
      rng.gauss_lex_global(src);

      fio->write_file(&src, filename_source);
    }

    // find reference or calculate from source.

    std::ifstream fref(filename_reference.c_str());
    if (fref) {
      vout.general(vl, "%s: reference file exists.\n", test_name.c_str(), filename_reference.c_str());

      fref.close();

      fio->read_file(&ref, filename_reference);
    } else {
      vout.general(vl, "%s: reference data not yet exported.\n", test_name.c_str());

      unique_ptr<FFT> fft(FFT::New("auto"));
      fft->fft(ref, src, FFT::FORWARD);

      fio->write_file(&ref, filename_reference);
    }


    int err = 0;

    // test FFT_xyz_1dim
    if (
      (NPEt == 1) &&
      ((NPEx * NPEy == 1) || (NPEx * NPEz == 1) || (NPEy * NPEz == 1))
      ) {
      vout.general(vl, "test FFT_xyz_1dim.\n");

      unique_ptr<FFT> fft(FFT::New("FFT_xyz_1dim"));
      fft->set_parameters(params_fft);

      err += run_task(fft, src, ref, res1, res2, check_precision, vl);
    } else {
      vout.general(vl, "FFT_xyz_1dim is not available.\n");
    }

    // test FFT_xyz_3dim
    if (
      (NPEt == 1) && (NPEx == 1) && (NPEy == 1)
      ) {
      vout.general(vl, "test FFT_xyz_3dim.\n");

      unique_ptr<FFT> fft(FFT::New("FFT_xyz_3dim"));
      fft->set_parameters(params_fft);

      err += run_task(fft, src, ref, res1, res2, check_precision, vl);
    } else {
      vout.general(vl, "FFT_xyz_3dim is not available.\n");
    }

    // test FFT_3d_local
    if (NPEx * NPEy * NPEz == 1)
    {
      vout.general(vl, "test FFT_3d_local.\n");

      unique_ptr<FFT> fft(FFT::New("FFT_3d_local"));
      fft->set_parameters(params_fft);

      err += run_task(fft, src, ref, res1, res2, check_precision, vl);
    } else {
      vout.general(vl, "FFT_3d_local is not available.\n");
    }

    // test FFT_3d_parallel1d
#ifdef USE_MPI
    if (NPEx * NPEy == 1)
    {
      vout.general(vl, "test FFT_3d_parallel_1dim.\n");

      unique_ptr<FFT> fft(FFT::New("FFT_3d_parallel_1dim"));
      fft->set_parameters(params_fft);

      err += run_task(fft, src, ref, res1, res2, check_precision, vl);
    } else {
      vout.general(vl, "FFT_3d_parallel_1dim is not available.\n");
    }
#endif

    // test FFT_3d_parallel3d
#ifdef USE_MPI
    if (true)
    {
      vout.general(vl, "test FFT_3d_parallel_3dim.\n");

      unique_ptr<FFT> fft(FFT::New("FFT_3d_parallel_3dim"));
      fft->set_parameters(params_fft);

      err += run_task(fft, src, ref, res1, res2, check_precision, vl);
    } else {
      vout.general(vl, "FFT_3d_parallel_3dim is not available.\n");
    }
#endif

    // test FFT_3d auto
    if (true)
    {
      vout.general(vl, "test FFT auto.\n");

      unique_ptr<FFT> fft(FFT::New("auto"));
      fft->set_parameters(params_fft);

      err += run_task(fft, src, ref, res1, res2, check_precision, vl);
    } else {
      vout.general(vl, "FFT auto is not available.\n");
    }


    if (do_check) {
      return Test::verify(err, 0);
    } else {
      vout.detailed(vl, "check skipped: expected_result not set.\n\n");
      return EXIT_SKIP;
    }
  }
} // namespace Test_FFT
#endif
