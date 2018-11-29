/*
@file    $Id: testlist.h #$

@brief

@author  Tatsumi Aoyama  (aoym)
$LastChangedBy: aoym $

@date    $LastChangedDate: 2013-04-08 18:00:27 #$

@version $LastChangedRevision: 1593 $
*/

#ifndef TESTLIST_INCLUDED
#define TESTLIST_INCLUDED

//- prototype declarations

namespace Test_Eigensolver {
    int solve(void);
}
namespace Test_Eigensolver_Chebyshev {
    int solve_chebyshev(void);
}
namespace Test_Eigensolver_Clover_SF {
    int solve_SF(void);
}


#ifdef USE_FFTWLIB
namespace Test_FFT {
    int fft(void);
}
#endif


namespace Test_Gauge {
    int plaquette_lex(void);
    int plaquette_eo(void);
    int shift(void);
}


namespace Test_GradientFlow {
    int run_test_RK1(void);
    int run_test_RK2(void);
    int run_test_RK3(void);
    int run_test_RK4(void);
    int run_test_RK_adaptive(void);
}


namespace Test_HMC_Clover_Isochemical {
    int update_Nf2(void);
    int RHMC_Nf2p1(void);
}
namespace Test_HMC_Clover {
    //  int update_Nf2(void);
    int run_test(void);
    int run_test_HYP(void);
    int update_Nf2_eo(void);
    int RHMC_Nf2p1(void);
    int RHMC_Nf2p1_eo(void);
}
namespace Test_HMC_Clover_SF {
    int update_Nf2(void);
    int RHMC_Nf2p1(void);
}
namespace Test_HMC_Quenched {
    int leapfrog_Nf0(void);
    int update_Nf0(void);
}
namespace Test_HMC_Wilson {
    int leapfrog_Nf2(void);
    int update_Nf2(void);
}


namespace Test_HotStart {
    int determinant(void);
    int eigenvalue(void);
    int unitary(void);
}


namespace Test_IO_Data {
    int test_io_data_text(void);
}
namespace Test_IO_GaugeConfig {
    int test_io_gconf_binary(void);

#ifdef USE_MPI
    int test_io_gconf_binary_distributed(void);
    int test_io_gconf_binary_parallel(void);
#endif
    int test_io_gconf_fortran(void);
    int test_io_gconf_text(void);

#ifdef USE_LIMELIB
    int test_io_gconf_ILDG(void);

#ifdef USE_MPI
    int test_io_gconf_ILDG_parallel(void);
#endif
#endif
}


namespace Test_Mult {
    int mult_Clover(void);
    int mult_CloverGeneral(void);
    int mult_Clover_Isochemical(void);

    //- NB. Fopr_Clover_SF is implemented only in Chiral rep.
    // int mult_Clover_SF(void);

    //- NB. test_Mult_Wilson is implemented separately for beginners
    // int mult_Wilson(void);
    int mult_WilsonGeneral(void);
    int mult_Wilson_Isochemical(void);

    //- NB. Fopr_Wilson_SF is implemented only in Chiral rep.
    // int mult_Wilson_SF(void);
}
namespace Test_Mult_eo {
    int mult_Clover_eo(void);
    int mult_Wilson_eo(void);
}
namespace Test_Mult_Wilson {
    int mult(void);
}


namespace Test_PolyakovLoop {
    int polyakovloop(void);
}


namespace Test_QuarkNumSuscept {
    int quark_num_suscept(void);
}


namespace Test_RandomNumbers_Mseries {
    int uniform_calc_pi(void);
    int gaussian(void);
    int gaussian_field(void);
    int test_global(void);
}
namespace Test_RandomNumbers_MT19937 {
    int uniform_calc_pi(void);
    int gaussian_field(void);
    int test_global(void);
}
#ifdef USE_SFMTLIB
namespace Test_RandomNumbers_SFMT {
    int uniform_calc_pi(void);
    int gaussian_field(void);
    int test_global(void);
}
#endif


namespace Test_Rational {
    int approx(void);
    int inverse(void);
    int smeared_rational(void);
}


namespace Test_SF_fAfP {
    int boundary_meson_2ptFunction(void);
}


namespace Test_ShiftSolver {
    int solve(void);
}


namespace Test_Spectrum {
    int hadron_2ptFunction_Clover(void);
    int hadron_2ptFunction_CloverGeneral(void);

    //- NB. test_Spectrum_Wilson is implemented separately for beginners
    // int hadron_2ptFunction_Wilson(void);
    int hadron_2ptFunction_Wilson_WallSource(void);
    int hadron_2ptFunction_WilsonGeneral(void);

    int hadron_2ptFunction_withFileIO_Clover(void);
    int hadron_2ptFunction_withFileIO_CloverGeneral(void);
    int hadron_2ptFunction_withFileIO_Wilson(void);
    int hadron_2ptFunction_withFileIO_WilsonGeneral(void);

    int hadron_2ptFunction_eo_Clover(void);
    int hadron_2ptFunction_eo_Wilson(void);

    int hadron_2ptFunction_eo_withFileIO_Clover();
    int hadron_2ptFunction_eo_withFileIO_Wilson();
}

namespace Test_Spectrum_CRSMatrix {
    int CRSsolver(void);
    int clover_lex(void);
}
namespace Test_Spectrum_Wilson {
    int hadron_2ptFunction(void);
}


namespace Test_TopologicalCharge {
    int topological_charge(void);
}


namespace Test_WilsonLoop {
    int wilsonloop(void);
}
#endif
