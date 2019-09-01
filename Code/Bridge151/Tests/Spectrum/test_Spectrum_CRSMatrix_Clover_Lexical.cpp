/*!
        @file    $Id: test_Spectrum_CRSMatrix_Clover_Lexical.cpp #$

        @brief

        @author  Hideo Matsufuru  (matsufuru)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate: 2013-01-22 22:08:29 #$

        @version $LastChangedRevision: 1571 $
*/

#include "Tests/test.h"

#include "Fopr/fopr_CRS.h"
#include "Fopr/fopr_Clover.h"

#include "Measurements/Fermion/source.h"

#include "Solver/solver_CG.h"

#include "IO/gaugeConfig.h"

//====================================================================
//! Test of CRS matrix format.

/*!
    This class generates and tests the CRS matrix format data.
    Matrix data in CRS format are generated for the following
    fermion operators together with sample source and solution
    vectors.

    The egenrated data can be tested using the function
    CRSsolver() whether really solved result satisfy the linear
    equation.
                                [24 Dec 2011 H.Matsufuru]
    (Coding history will be recovered from trac.)
    YAML is implemented.        [14 Nov 2012 Y.Namekawa]
    unique_ptr is introduced to avoid memory leaks
                                [21 Mar 2015 Y.Namekawa]
 */

namespace Test_Spectrum_CRSMatrix {
  //- test-private parameters
  namespace {
    const std::string filename_input = "test_Spectrum_CRSMatrix_Clover_Lexical.yaml";
  }

  //- prototype declaration
  int clover_lex(void);

  int CRSsolver(
    const string& solution,
    const string& matrix,
    const string& source,
    double& result       /* return value */
    );

  void write_text(const Field& f, const string& filename);
  void read_text(Field& f, const string& filename);


#ifdef USE_TESTMANAGER_AUTOREGISTER
#ifdef USE_MPI
  // this test runs only in single-node environment.
#else
  namespace {
#if defined(USE_GROUP_SU2)
    // Nc=2 is not available.
#else
    static const bool is_registered = TestManager::RegisterTest(
      "Spectrum.CRSMatrix.Clover_Lexical",
      clover_lex
      );
#endif
  }
#endif
#endif

  //====================================================================
  int clover_lex(void)
  {
    if (CommonParameters::NPE() > 1) {
      // run in parallel: unsupported.
      vout.general("CRSMatrix test: node-parallel unsupported. skip.\n");
      return EXIT_SKIP;
    }

    // ####  parameter setup  ####
    int Nc   = CommonParameters::Nc();
    int Nd   = CommonParameters::Nd();
    int Ndim = CommonParameters::Ndim();
    int Nvol = CommonParameters::Nvol();

    Parameters params_all = ParameterManager::read(filename_input);

    Parameters params_test = params_all.lookup("Test_Spectrum_CRSMatrix");

    const string str_gconf_read = params_test.get_string("gauge_config_type_input");
    const string readfile       = params_test.get_string("config_filename_input");
    const string matrix_file    = params_test.get_string("matrix_output");
    const string source_file    = params_test.get_string("source_output");
    const string solution_file  = params_test.get_string("solution_output");
    const string str_vlevel     = params_test.get_string("verbose_level");

    const bool   do_check        = params_test.is_set("expected_result");
    const double expected_result = do_check ? params_test.get_double("expected_result") : 0.0;

    Bridge::VerboseLevel vl = vout.set_verbose_level(str_vlevel);

    //- print input parameters
    vout.general(vl, "  gconf_read      = %s\n", str_gconf_read.c_str());
    vout.general(vl, "  readfile        = %s\n", readfile.c_str());
    vout.general(vl, "  matrix_output   = %s\n", matrix_file.c_str());
    vout.general(vl, "  source_output   = %s\n", source_file.c_str());
    vout.general(vl, "  solution_output = %s\n", solution_file.c_str());
    vout.general(vl, "  vlevel          = %s\n", str_vlevel.c_str());
    vout.general(vl, "\n");

    //- input parameter check
    int err = 0;
    err += ParameterCheck::non_NULL(str_gconf_read);
    err += ParameterCheck::non_NULL(readfile);
    err += ParameterCheck::non_NULL(matrix_file);
    err += ParameterCheck::non_NULL(source_file);
    err += ParameterCheck::non_NULL(solution_file);

    if (err) {
      vout.crucial(vl, "Error: Test_Spectrum_CRSMatrix: input parameters have not been set\n");
      exit(EXIT_FAILURE);
    }


    // ####  Set up a gauge configuration  ####
    unique_ptr<Field_G>     U(new Field_G(Nvol, Ndim));
    unique_ptr<GaugeConfig> gconf_read(new GaugeConfig(str_gconf_read));
    gconf_read->read_file(U, readfile);


    // ####  object setup  #####
    double           kappa = 0.12;
    double           cSW   = 1.0;
    std::vector<int> boundary(Ndim);
    boundary[0] = -1;
    boundary[1] = -1;
    boundary[2] = -1;
    boundary[3] = -1;

    Parameters params_c;
    params_c.set_double("hopping_parameter", kappa);
    params_c.set_double("clover_coefficient", cSW);
    params_c.set_int_vector("boundary_condition", boundary);

    unique_ptr<Fopr> fopr_c(Fopr::New("Clover", "Dirac"));
    fopr_c->set_parameters(params_c);
    fopr_c->set_config(U);
    fopr_c->set_mode("D");

    unique_ptr<Fopr_CRS> fopr_crs(new Fopr_CRS(fopr_c.get()));
    fopr_crs->write_matrix(matrix_file);

    int                Niter     = 100;
    int                Nrestart  = 40;
    double             Stop_cond = 1.0e-28;
    unique_ptr<Solver> solver(new Solver_CG(fopr_c));
    solver->set_parameters(Niter, Nrestart, Stop_cond);

    std::vector<int> source_position(Ndim);
    source_position[0] = 0;
    source_position[1] = 0;
    source_position[2] = 0;
    source_position[3] = 0;

    unique_ptr<Source> source(Source::New("Local"));
    Parameters         params_source;
    params_source.set_int_vector("source_position", source_position);
    source->set_parameters(params_source);


    // ####  Execution main part  ####
    std::vector<Field_F> sq(Nc * Nd);
    Field_F              xq, b, b2;

    for (int ispin = 0; ispin < Nd; ++ispin) {
      for (int icolor = 0; icolor < Nc; ++icolor) {
        int idx = icolor + Nc * ispin;
        source->set(b, idx);

        write_text(b, source_file.c_str());

        fopr_c->set_mode("D");
        fopr_c->mult_dag(b2, b);

        int    Nconv;
        double diff_CG;
        fopr_c->set_mode("DdagD");
        solver->solve(xq, b2, Nconv, diff_CG);
        vout.general(vl, "  ispin = %2d  icolor = %2d  Nconv = %4d  diff = %12.6e\n",
                     ispin, icolor, Nconv, diff_CG);

        write_text(xq, solution_file.c_str());

        Field y(b);
        fopr_c->set_mode("D");
        fopr_c->mult(y, xq);
        axpy(y, -1.0, b); // y -= b;
        scal(y, -1.0);    //y *= -1.0;

        double yy = y.norm2();
        vout.general(vl, "    standard norm2 = %.8e\n", yy);

        sq[icolor + Nc * ispin] = xq;
      }
    }

#if 0
    // ### Solver for CRS version ###
    unique_ptr<Solver> solver_crs(new Solver_CG(fopr_crs2));
    solver_crs->set_parameters(Niter, Stop_cond);

    double yy = 0.0L;
    {
      int ispin = 0;
      {
        int icolor = 0;
        //for(int ispin = 0; ispin < Nd; ++ispin){
        // for(int icolor = 0; icolor < Nc; ++icolor){

        source->set(b, icolor, ispin);

        fopr_crs2->set_mode("Ddag");
        b2 = (Field_F)fopr_crs2->mult((Field)b);

        int    Nconv;
        double diff_CG;
        fopr_crs2->set_mode("DdagD");
        solver_crs->solve(xq, b2, Nconv, diff_CG);
        vout.general(vl, "  ispin = %2d  icolor = %2d  Nconv = %4d  diff = %12.6e\n",
                     ispin, icolor, Nconv, diff_CG);

        Field_F y(b);
        fopr_crs2->set_mode("D");
        y -= (Field_F)fopr_crs2->mult((Field)xq);
        yy = y.norm2();
        vout.general(vl, "    standard norm2 = %.8e\n", yy);

        sq[icolor + Nc * ispin] = xq;
      }
    }
#endif


    double result;
    CRSsolver(solution_file, matrix_file, source_file, result);


    if (do_check) {
      return Test::verify(result, expected_result);
    } else {
      vout.detailed(vl, "check skipped: expected_result not set.\n\n");
      return EXIT_SKIP;
    }
  }
} // namespace Test_Spectrum_CRSMatrix
