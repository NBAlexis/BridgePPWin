/*!
        @file    test_Spectrum_2ptFunction_withFileIO_initial_guess.cpp

        @brief

        @author  Yusuke Namekawa  (namekawa)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/
#include "BppSmallTest.h"
#include "test.h"

#include "Field/shiftField_lex.h"

#include "Fopr/fopr_Smeared.h"

#include "IO/gaugeConfig.h"

#include "Measurements/Fermion/corr2pt_4spinor.h"
#include "Measurements/Fermion/fprop_Standard_lex.h"
#include "Measurements/Fermion/source.h"
#include "Measurements/Gauge/gaugeFixing.h"

//#include "Tools/filename.h"
#include "Tools/file_utils.h"
#include "Tools/randomNumberManager.h"
#include "Tools/gammaMatrixSet.h"

//====================================================================
//! Test of Spectroscopy with fileIO.

/*!
    This test class obtains the quark propagator for clover fermion
    and calculates typical hadron correlators.
    The quantum numbers of hadrons are specified with gamma matrices
    (GammaMatrix class instance) whose set is defined in a subclass
    of GammaMatrixSet class.
                                             [12 Apr 2012 H.Matsufuru]
    (Coding history will be recovered from trac.)
    Selectors are implemented.               [02 Feb 2013 Y.Namekawa]
    Smearing is implemented.                 [03 Mar 2013 Y.Namekawa]
    (Selectors are replaced with factories by Aoyama-san)
    Number_of_valence_quarks is implemented. [15 May 2013 Y.Namekawa]
    unique_ptr is introduced to avoid memory leaks
                                             [21 Mar 2015 Y.Namekawa]
    Add shift_origin of gauge config         [21 May 2017 Y.Namekawa]
    Modify to choose a solver element_type quark by quark
                                             [23 Apr 2018 Y.Namekawa]
    Employ initial guess for solver          [ 1 May 2018 Y.Namekawa]
 */

namespace Test_Spectrum {
  const std::string test_name = "Spectrum.Hadron2ptFunction_withFileIO_initial_guess";

  //- test-private parameters
  namespace {
    // const std::string filename_input  = "test_Spectrum_Clover_Hadron2ptFunction.yaml";
  }

  //- prototype declaration
  int hadron_2ptFunction_withFileIO_initial_guess(const std::string&);

  int calculate_quark_propagator_withFileIO_initial_guess(const std::string&);
  int calculate_hadron_correlator_withFileIO_initial_guess(const std::string&);

  //- hadron_2ptFunction for various fermions
  int hadron_2ptFunction_withFileIO_Clover_initial_guess()
  {
    return hadron_2ptFunction_withFileIO_initial_guess("test_Spectrum_Clover_Hadron2ptFunction_initial_guess.yaml");
  }


  // int hadron_2ptFunction_withFileIO_CloverGeneral()
  // {
  //   return hadron_2ptFunction_withFileIO("test_Spectrum_CloverGeneral_Hadron2ptFunction.yaml");
  // }


  // int hadron_2ptFunction_withFileIO_TMWilson()
  // {
  //   return hadron_2ptFunction_withFileIO("test_Spectrum_TMWilson_Hadron2ptFunction.yaml");
  // }


  // int hadron_2ptFunction_withFileIO_Wilson()
  // {
  //   return hadron_2ptFunction_withFileIO("test_Spectrum_Wilson_Hadron2ptFunction.yaml");
  // }


  // int hadron_2ptFunction_withFileIO_WilsonGeneral()
  // {
  //   return hadron_2ptFunction_withFileIO("test_Spectrum_WilsonGeneral_Hadron2ptFunction.yaml");
  // }


#ifdef USE_TESTMANAGER_AUTOREGISTER
  namespace {
#if defined(USE_GROUP_SU2)
    // Nc=2 is not available.
#else
    static const bool is_registered_Clover = TestManager::RegisterTest(
      "Spectrum.Clover.Hadron2ptFunction_withFileIO_initial_guess",
      hadron_2ptFunction_withFileIO_Clover_initial_guess
      );
#endif
  }
#endif

  //====================================================================
  int hadron_2ptFunction_withFileIO_initial_guess(const std::string& filename_input)
  {
    int result = 0;

    const unique_ptr<Timer> timer(new Timer(test_name));


    // ####  Execution main part  ####
    timer->start();

    result += calculate_quark_propagator_withFileIO_initial_guess(filename_input);
    result += calculate_hadron_correlator_withFileIO_initial_guess(filename_input);

    timer->report();

    return result;
  }


  //====================================================================
  int calculate_quark_propagator_withFileIO_initial_guess(const std::string& filename_input)
  {
    // ####  parameter setup  ####
    const int Nc   = CommonParameters::Nc();
    const int Nd   = CommonParameters::Nd();
    const int Ndim = CommonParameters::Ndim();
    const int Nvol = CommonParameters::Nvol();

    const Parameters params_all = ParameterManager::read(filename_input);

    const Parameters params_test = params_all.lookup("Test_Spectrum");

    const int N_quark = params_test.get_int("number_of_valence_quarks");

    const std::string str_gconf_status = params_test.get_string("gauge_config_status");
    const std::string str_gconf_read   = params_test.get_string("gauge_config_type_input");
    const std::string readfile         = params_test.get_string("config_filename_input");
    const vector<int> Nshift_origin    = params_test.get_int_vector("shift_origin");
    const std::string str_vlevel       = params_test.get_string("verbose_level");

    const std::string str_gfix_type  = params_all.lookup("GaugeFixing").get_string("gauge_fixing_type");
    const std::string str_proj_type  = params_all.lookup("Projection").get_string("projection_type");
    const std::string str_smear_type = params_all.lookup("Smear").get_string("smear_type");

    const Bridge::VerboseLevel vl = vout.set_verbose_level(str_vlevel);

    //- print input parameters
    vout.general(vl, "  gconf_status = %s\n", str_gconf_status.c_str());
    vout.general(vl, "  gconf_read   = %s\n", str_gconf_read.c_str());
    vout.general(vl, "  readfile     = %s\n", readfile.c_str());
    for (int mu = 0; mu < Ndim; ++mu) {
      vout.general(vl, "  shift_origin[%d] = %d\n", mu, Nshift_origin[mu]);
    }
    vout.general(vl, "  vlevel       = %s\n", str_vlevel.c_str());
    vout.general(vl, "  gfix_type    = %s\n", str_gfix_type.c_str());
    vout.general(vl, "  proj_type    = %s\n", str_proj_type.c_str());
    vout.general(vl, "  smear_type   = %s\n", str_smear_type.c_str());

    vector<Parameters> params_quark;
    char tmp[128];
    for (int iq = 0; iq < N_quark; ++iq) {
        memset(tmp, sizeof(char) * 128, 0);
        _itoa_s(iq + 1, tmp, 10);
        string qlabel = "Quark_" + string(tmp);
      params_quark.push_back(params_all.lookup(qlabel));
    }

    // NB. all str_gmset_type are supposed to be the same.
    std::string str_gmset_type = params_quark[0].lookup("Fopr").get_string("gamma_matrix_type");
    vout.general(vl, "  gmset_type  = %s\n", str_gmset_type.c_str());

    std::vector<std::string> savefile_base(N_quark);
    for (int iq = 0; iq < N_quark; ++iq) {
      savefile_base[iq] = params_quark[iq].get_string("temporary_filename_base");
    }

    for (int iq = 0; iq < N_quark; ++iq) {
      std::string str_solver_type = params_quark[iq].lookup("Solver").get_string("solver_type");
      std::string str_source_type = params_quark[iq].lookup("Source").get_string("source_type");
      vector<int> source_position = params_quark[iq].lookup("Source").get_int_vector("source_position");

      vout.general(vl, "  Quark_%d:\n", iq + 1);
      vout.general(vl, "    savefile_base[iq] = %s[%d]\n", savefile_base[iq].c_str(), iq);
      vout.general(vl, "    solver_type       = %s\n", str_solver_type.c_str());
      vout.general(vl, "    source_type       = %s\n", str_source_type.c_str());
      for (int mu = 0; mu < Ndim; ++mu) {
        vout.general(vl, "    source_position[%d] = %d\n", mu, source_position[mu]);
      }
    }
    vout.general(vl, "\n");


    //- input parameter check
    int err = 0;
    err += ParameterCheck::non_NULL(str_gconf_status);

    if (err) {
      vout.crucial(vl, "Error at %s: input parameters have not been set\n", test_name.c_str());
      exit(EXIT_FAILURE);
    }

    for (int iq = 0; iq < N_quark; ++iq) {
      std::string str_solver_type = params_quark[iq].lookup("Solver").get_string("solver_type");

      if (str_solver_type == "CG") {
        vout.crucial(vl, "Error at %s: CG can not be adopted. Use CGNE,CGNR, instead.\n", test_name.c_str());
        exit(EXIT_FAILURE);
      }
    }

    RandomNumberManager::initialize("Mseries", 1234567UL);


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

    {
      ShiftField_lex shift;
      for (int i_dir = 0; i_dir < Ndim; ++i_dir) {
        for (int i_shift = 0; i_shift < Nshift_origin[i_dir]; ++i_shift) {
          unique_ptr<Field_G> Ushift(new Field_G(Nvol, Ndim));
          shift.backward(*Ushift, *U, i_dir);
          copy(*U, *Ushift);
        }
      }
    }

    {
      unique_ptr<Projection> proj(Projection::New(str_proj_type));
      proj->set_parameters(params_all.lookup("Projection"));

      unique_ptr<Smear> smear(Smear::New(str_smear_type, proj));
      smear->set_parameters(params_all.lookup("Smear"));

      const unique_ptr<Director_Smear> dr_smear(new Director_Smear(smear));
      dr_smear->set_parameters(params_all.lookup("Director_Smear"));
      dr_smear->set_config(U);

      const int     Nsmear  = dr_smear->get_Nsmear();
      const Field_G *Usmear = (Field_G *)dr_smear->getptr_smearedConfig(Nsmear);

      copy(*U, *Usmear);
    }


    // ####  Gauge fixing  ####
    {
      unique_ptr<Field_G>           Ufix(new Field_G(Nvol, Ndim));
      const unique_ptr<GaugeFixing> gfix(GaugeFixing::New(str_gfix_type));
      gfix->set_parameters(params_all.lookup("GaugeFixing"));

      gfix->fix(*Ufix, *U);

      copy(*U, *Ufix);
    }


    // ####  object setup  #####
#ifdef LIB_CPP11
    std::vector<unique_ptr<Fopr> >   fopr(N_quark);
    std::vector<unique_ptr<Solver> > solver(N_quark);
    std::vector<unique_ptr<Fprop> >  fprop_lex(N_quark);
    std::vector<unique_ptr<Source> > source(N_quark);

    for (int iq = 0; iq < N_quark; ++iq) {
      std::string str_fopr_type = params_quark[iq].lookup("Fopr").get_string("fermion_type");
      fopr[iq].reset(Fopr::New(str_fopr_type, str_gmset_type));
      fopr[iq]->set_parameters(params_quark[iq].lookup("Fopr"));
      fopr[iq]->set_config(U);

      std::string str_solver_type = params_quark[iq].lookup("Solver").get_string("solver_type");
      solver[iq].reset(Solver::New(str_solver_type, fopr[iq]));
      solver[iq]->set_parameters(params_quark[iq].lookup("Solver"));

      fprop_lex[iq].reset(new Fprop_Standard_lex(solver[iq]));

      std::string str_source_type = params_quark[iq].lookup("Source").get_string("source_type");
      source[iq].reset(Source::New(str_source_type));
      source[iq]->set_parameters(params_quark[iq].lookup("Source"));
    }
#else
    std::vector<Fopr *>   fopr(N_quark);
    std::vector<Solver *> solver(N_quark);
    std::vector<Fprop *>  fprop_lex(N_quark);
    std::vector<Source *> source(N_quark);

    for (int iq = 0; iq < N_quark; ++iq) {
      std::string str_fopr_type = params_quark[iq].lookup("Fopr").get_string("fermion_type");
      fopr[iq] = Fopr::New(str_fopr_type, str_gmset_type);
      fopr[iq]->set_parameters(params_quark[iq].lookup("Fopr"));
      fopr[iq]->set_config(U);

      std::string str_solver_type = params_quark[iq].lookup("Solver").get_string("solver_type");
      solver[iq] = Solver::New(str_solver_type, fopr[iq]);
      solver[iq]->set_parameters(params_quark[iq].lookup("Solver"));

      fprop_lex[iq] = new Fprop_Standard_lex(solver[iq]);

      std::string str_source_type = params_quark[iq].lookup("Source").get_string("source_type");
      source[iq] = Source::New(str_source_type);
      source[iq]->set_parameters(params_quark[iq].lookup("Source"));
    }
#endif
    vout.general(vl, "\n");

    const unique_ptr<FieldIO> field_io(new FieldIO_Binary(IO_Format::Trivial));


    // ####  Execution main part  ####
    for (int iq = 0; iq < N_quark; ++iq) {
      vout.general(vl, "Solving quark propagator, flavor = %d:\n", iq + 1);
      vout.general(vl, "  color spin   Nconv      diff           diff2\n");

      for (int ispin = 0; ispin < Nd; ++ispin) {
        for (int icolor = 0; icolor < Nc; ++icolor) {
          int i_cd = icolor + Nc * ispin;

          Field_F b;  // b.set(0.0);
          source[iq]->set(b, i_cd);

          Field_F xq;
          int     Nconv;
          double  diff;

          std::string filename = FileUtils::generate_filename("%s_c%d_s%d.dat", savefile_base[iq].c_str(), (icolor + 1), (ispin + 1));
          field_io->read_file(&xq, filename);

          fprop_lex[iq]->invert_D(xq, b, Nconv, diff);

          Field_F y(b);
          fopr[iq]->set_mode("D");
          fopr[iq]->mult(y, xq);  // y  = fopr[iq]->mult(xq);
          axpy(y, -1.0, b);       // y -= b;
          double diff2 = y.norm2() / b.norm2();

          vout.general(vl, "   %2d   %2d   %6d   %12.4e   %12.4e\n",
                       icolor, ispin, Nconv, diff, diff2);

          field_io->write_file(&xq, filename);
        }
      }

      vout.general(vl, "\n");
    }

#ifdef LIB_CPP11
    // do nothing
#else
    // tidy-up
    for (int iq = 0; iq < N_quark; ++iq) {
      delete fopr[iq];
      delete solver[iq];
      delete fprop_lex[iq];
      delete source[iq];
    }
#endif


    return EXIT_SUCCESS;
  }


  //====================================================================
  int calculate_hadron_correlator_withFileIO_initial_guess(const std::string& filename_input)
  {
    // ####  parameter setup  ####
    const int Nc   = CommonParameters::Nc();
    const int Nd   = CommonParameters::Nd();
    const int Ndim = CommonParameters::Ndim();
    const int Nvol = CommonParameters::Nvol();

    const Parameters params_all = ParameterManager::read(filename_input);

    const Parameters params_test = params_all.lookup("Test_Spectrum");

    const int N_quark = params_test.get_int("number_of_valence_quarks");

    const std::string str_vlevel = params_test.get_string("verbose_level");

    const bool   do_check        = params_test.is_set("expected_result");
    const double expected_result = do_check ? params_test.get_double("expected_result") : 0.0;

    std::vector<std::string> savefile_base(N_quark);

    char tmp[128];
    for (int iq = 0; iq < N_quark; ++iq) {
        memset(tmp, sizeof(char) * 128, 0);
        _itoa_s(iq + 1, tmp, 10);
        string qlabel = "Quark_" + string(tmp);
      savefile_base[iq] = params_all
                          .lookup(qlabel)
                          .get_string("temporary_filename_base");
    }

    const Bridge::VerboseLevel vl = vout.set_verbose_level(str_vlevel);

    //- print input parameters
    vout.general(vl, "  vlevel      = %s\n", str_vlevel.c_str());

    // NB. all str_gmset_type are supposed to be the same.
    const std::string str_gmset_type = params_all.lookup("Quark_1").lookup("Fopr").get_string("gamma_matrix_type");
    vout.general(vl, "  gmset_type  = %s\n", str_gmset_type.c_str());

    for (int iq = 0; iq < N_quark; ++iq) {
      vout.general(vl, "  savefile_base[iq] = %s[%d]\n", savefile_base[iq].c_str(), iq);
    }

    // ####  Set up objects  #####
    unique_ptr<GammaMatrixSet> gmset(GammaMatrixSet::New(str_gmset_type));

    const unique_ptr<FieldIO> field_io(new FieldIO_Binary(IO_Format::Trivial));

    Corr2pt_4spinor corr(gmset);
    corr.set_parameters(params_all.lookup("Corr2pt_4spinor"));


    // ####  Execution main part  ####
    typedef std::vector<Field_F> PropagatorSet;
    std::vector<PropagatorSet> sq(N_quark);
    for (int iq = 0; iq < N_quark; ++iq) {
      sq[iq].resize(Nc * Nd);

      for (int i_cd = 0; i_cd < Nc * Nd; ++i_cd) {
        sq[iq][i_cd].set(0.0);
      }
    }


    for (int iq = 0; iq < N_quark; ++iq) {
      for (int ispin = 0; ispin < Nd; ++ispin) {
        for (int icolor = 0; icolor < Nc; ++icolor) {
          int i_cd = icolor + Nc * ispin;

          std::string filename = FileUtils::generate_filename("%s_c%d_s%d.dat", savefile_base[iq].c_str(), (icolor + 1), (ispin + 1));
          field_io->read_file(&sq[iq][i_cd], filename);
        }
      }
    }


    //- meson correlators
    std::vector<double> result(N_quark);

    vout.general(vl, "2-point correlator:\n");

    //- case(iq_1 == iq_2)
    for (int iq = 0; iq < N_quark; ++iq) {
      vout.general(vl, "Flavor combination = %d, %d\n", iq + 1, iq + 1);
      result[iq] = corr.meson_all(sq[iq], sq[iq]);
      vout.general(vl, "\n");
    }


    //- case(iq_1 < iq_2)
    for (int iq = 0; iq < N_quark; ++iq) {
      for (int jq = iq + 1; jq < N_quark; ++jq) {
        vout.general(vl, "Flavor combination = %d, %d\n", iq + 1, jq + 1);
        //double result_2 = corr.meson_all(sq[iq], sq[jq]);
        vout.general(vl, "\n");
      }
    }


    RandomNumberManager::finalize();

#ifdef LIB_CPP11
    // do nothing
#else
    // tidy-up
#endif

    if (do_check) {
      return Test::verify(result[0], expected_result);
    } else {
      vout.detailed(vl, "check skipped: expected_result not set.\n\n");
      return EXIT_SKIP;
    }
  }
} // namespace Test_Spectrum