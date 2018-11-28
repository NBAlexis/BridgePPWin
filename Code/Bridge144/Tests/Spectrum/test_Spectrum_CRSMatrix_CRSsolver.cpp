/*!
        @file    $Id: test_Spectrum_CRSMatrix_CRSsolver.cpp #$

        @brief

        @author  Hideo Matsufuru  (matsufuru)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate: 2013-01-22 22:08:29 #$

        @version $LastChangedRevision: 1561 $
*/

#include "Tests/test.h"

#include "Fopr/fopr_CRS.h"
#include "Solver/solver_CG.h"

//====================================================================
//     // clover set
//     string solution = "solution4x8_clover.crs";
//     string matrix   = "matrix4x8_clover.crs";
//     string source   = "source4x8_clover.crs";

//     // 5d-overlap set
//     string solution = "solution4x4_ov5d.crs";
//     string matrix   = "matrix4x4_ov5d.crs";
//     string source   = "source4x4_ov5d.crs";

//     // domain-wall set
//     string solution = "solution4x4_dw.crs";
//     string matrix   = "matrix4x4_dw.crs";
//     string source   = "source4x4_dw.crs";


namespace Test_Spectrum_CRSMatrix
{
  void write_text(const Field& f, const string& filename);
  void read_text(Field& f, const string& filename);

  int CRSsolver(
    const string& solution,
    const string& matrix,
    const string& source,
    double& result       /* return value */
    )
  {
    Bridge::VerboseLevel vl = CommonParameters::Vlevel();

    // #####  Following part is common  #####

    // read source and solution vectors
    Field b, xq;

    read_text(b, source);
    read_text(xq, solution);

    // read CRS matrix
    unique_ptr<Fopr> fopr_crs(Fopr::New("CRS", matrix));


    // setup of CG solver
    int                Niter     = 100;
    int                Nrestart  = 40;
    double             Stop_cond = 1.0e-28;
    unique_ptr<Solver> solver(new Solver_CG(fopr_crs));
    solver->set_parameters(Niter, Nrestart, Stop_cond);

    Field b2(b);
    Field x(b);

    // setup of CGNE source
    fopr_crs->set_mode("D");
    fopr_crs->mult_dag(b2, b);

    // CGNE solver
    int    Nconv;
    double diff;
    fopr_crs->set_mode("DdagD");
    solver->solve(x, b2, Nconv, diff);
    vout.general(vl, "solver(CG):  Nconv = %4d  diff = %12.6e\n", Nconv, diff);

    // check
    axpy(x, -1.0, xq);
    double xx = x.norm2();
    vout.general(vl, "standard norm2 = %.8e\n", xx);

    result = xx;


    return EXIT_SUCCESS;
  }


//====================================================================
  void write_text(const Field& f, const string& filename)
  {
    // This function works only on single node.
    assert(CommonParameters::NPE() == 1);

    std::ofstream field_out(filename.c_str());
    if (!field_out.is_open()) {
      vout.crucial("Error: Failed to open the text file. %s\n", filename.c_str());
      exit(EXIT_FAILURE);
    }

    vout.general("Writing Field to %s\n", filename.c_str());

    field_out << f.nin() << std::endl;
    field_out << f.nvol() << std::endl;
    field_out << f.nex() << std::endl;

    field_out.setf(std::ios_base::scientific, std::ios_base::floatfield);
    field_out.precision(14);

    for (int j = 0, n = f.size(); j < n; ++j) {
      field_out << f.cmp(j) << std::endl;
    }

    field_out.close();

    vout.general("Writing Field finished.\n");
  }


//====================================================================
  void read_text(Field& f, const string& filename)
  {
    // This function works only on single node.
    assert(CommonParameters::NPE() == 1);

    std::ifstream field_in(filename.c_str());
    if (!field_in.is_open()) {
      vout.crucial("Error: Failed to open the text file. %s\n", filename.c_str());
      exit(EXIT_FAILURE);
    }

    vout.general("Reading Field from %s\n", filename.c_str());

    int Nin = 0, Nvol = 0, Nex = 0;
    field_in >> Nin;
    field_in >> Nvol;
    field_in >> Nex;

    f.reset(Nin, Nvol, Nex);

    double v;
    for (int j = 0, n = f.ntot(); j < n; ++j) {
      field_in >> v;
      f.set(j, v);
    }

    field_in.close();

    vout.general("Reading Field finished.\n");
  }
} // namespace Test_Spectrum_CRSMatrix
