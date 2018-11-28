/*
        @file    $Id: run_test.cpp #$

        @brief

        @author  Tatsumi Aoyama  (aoym)
                 $LastChangedBy: namekawa $

        @date    $LastChangedDate: 2013-04-08 18:00:27 #$

        @version $LastChangedRevision: 1605 $
*/

//#include "main.h"
#include <cstdlib>
#include "Tests/testlist.h"

//====================================================================
int run_test()
{
  //- put your exec test here

#if defined(USE_GROUP_SU2)
  // Nc=2 is not available.
#else
#ifdef USE_TEST_EIGENSOLVER
  Test_Eigensolver::solve();
  Test_Eigensolver_Chebyshev::solve_chebyshev();
  Test_Eigensolver_Clover_SF::solve_SF();
#endif

#ifdef USE_TEST_FFT
#ifdef USE_FFTWLIB
  Test_FFT::fft();
#endif
#endif

#ifdef USE_TEST_GAUGE
  Test_Gauge::plaquette_lex();
  Test_Gauge::plaquette_eo();
  Test_Gauge::shift();
#endif

#ifdef USE_TEST_GRADIENTFLOW
  Test_GradientFlow::run_test_RK1();
  Test_GradientFlow::run_test_RK2();
  Test_GradientFlow::run_test_RK3();
  Test_GradientFlow::run_test_RK4();
  Test_GradientFlow::run_test_RK_adaptive();
#endif

#ifdef USE_TEST_HMC
  // Test_HMC_Clover::update_Nf2();
  Test_HMC_Clover::run_test();
  Test_HMC_Clover::run_test_HYP();
  Test_HMC_Clover::update_Nf2_eo();
  Test_HMC_Clover::RHMC_Nf2p1();
  Test_HMC_Clover::RHMC_Nf2p1_eo();

  Test_HMC_Clover_Isochemical::update_Nf2();
  Test_HMC_Clover_Isochemical::RHMC_Nf2p1();

  Test_HMC_Clover_SF::update_Nf2();
  Test_HMC_Clover_SF::RHMC_Nf2p1();

  Test_HMC_Quenched::leapfrog_Nf0();
  Test_HMC_Quenched::update_Nf0();

  Test_HMC_Wilson::leapfrog_Nf2();
  Test_HMC_Wilson::update_Nf2();
#endif

#ifdef USE_TEST_HOTSTART
  Test_HotStart::determinant();
  Test_HotStart::eigenvalue();
  Test_HotStart::unitary();
#endif

#ifdef USE_TEST_IO
  Test_IO_Data::test_io_data_text();

  Test_IO_GaugeConfig::test_io_gconf_binary();
#ifdef USE_MPI
  Test_IO_GaugeConfig::test_io_gconf_binary_distributed();
  Test_IO_GaugeConfig::test_io_gconf_binary_parallel();
#endif
  Test_IO_GaugeConfig::test_io_gconf_fortran();
  Test_IO_GaugeConfig::test_io_gconf_text();

#ifdef USE_LIMELIB
  Test_IO_GaugeConfig::test_io_gconf_ILDG();
#ifdef USE_MPI
  Test_IO_GaugeConfig::test_io_gconf_ILDG_parallel();
#endif
#endif
#endif

#ifdef USE_TEST_MULT
  Test_Mult::mult_Clover();
  Test_Mult::mult_CloverGeneral();
  Test_Mult::mult_Clover_Isochemical();
  //- NB. Fopr_Clover_SF is implemented only in Chiral rep.
  // Test_Mult::mult_Clover_SF();
  //- NB. test_Mult_Wilson is implemented separately for beginners
  // Test_Mult::mult_Wilson();
  Test_Mult::mult_WilsonGeneral();
  Test_Mult::mult_Wilson_Isochemical();
  //- NB. Fopr_Wilson_SF is implemented only in Chiral rep.
  // Test_Mult::mult_Wilson_SF();

  Test_Mult_Wilson::mult();

  Test_Mult_eo::mult_Clover_eo();
  Test_Mult_eo::mult_Wilson_eo();
#endif

#ifdef USE_TEST_POLYAKOVLOOP
  Test_PolyakovLoop::polyakovloop();
#endif

#ifdef USE_TEST_QUARKNUMSUSCEPT
  Test_QuarkNumSuscept::quark_num_suscept();
#endif

#ifdef USE_TEST_RANDOMNUMBERS
  Test_RandomNumbers_Mseries::uniform_calc_pi();
  Test_RandomNumbers_Mseries::gaussian();
  Test_RandomNumbers_Mseries::gaussian_field();
  Test_RandomNumbers_Mseries::test_global();

  Test_RandomNumbers_MT19937::uniform_calc_pi();
  Test_RandomNumbers_MT19937::gaussian_field();
  Test_RandomNumbers_MT19937::test_global();

#ifdef USE_SFMTLIB
  Test_RandomNumbers_SFMT::uniform_calc_pi();
  Test_RandomNumbers_SFMT::gaussian_field();
  Test_RandomNumbers_SFMT::test_global();
#endif
#endif

#ifdef USE_TEST_RATIONAL
  Test_Rational::approx();
  Test_Rational::inverse();
  Test_Rational::smeared_rational();
#endif

#ifdef USE_TEST_SF_FAFP
#if defined(USE_OPENMP)
  // multithreading yet unsupported.
#else
  Test_SF_fAfP::boundary_meson_2ptFunction();
#endif
#endif

#ifdef USE_TEST_SHIFTSOLVER
  Test_ShiftSolver::solve();
#endif

#ifdef USE_TEST_SPECTRUM
  Test_Spectrum::hadron_2ptFunction_Clover();
  Test_Spectrum::hadron_2ptFunction_CloverGeneral();
  //- NB. test_Spectrum_Wilson is implemented separately for beginners
  // Test_Spectrum::hadron_2ptFunction_Wilson();
  Test_Spectrum::hadron_2ptFunction_Wilson_WallSource();
  Test_Spectrum::hadron_2ptFunction_WilsonGeneral();

  Test_Spectrum::hadron_2ptFunction_eo_withFileIO_Clover();
  //- NB. test_Spectrum_Wilson is implemented separately for beginners
  // Test_Spectrum::hadron_2ptFunction_eo_withFileIO_Wilson();

  Test_Spectrum::hadron_2ptFunction_withFileIO_Clover();
  Test_Spectrum::hadron_2ptFunction_withFileIO_CloverGeneral();
  //- NB. test_Spectrum_Wilson is implemented separately for beginners
  // Test_Spectrum::hadron_2ptFunction_withFileIO_Wilson();
  Test_Spectrum::hadron_2ptFunction_withFileIO_WilsonGeneral();

  Test_Spectrum::hadron_2ptFunction_eo_Clover();
  //- NB. test_Spectrum_Wilson is implemented separately for beginners
  // Test_Spectrum::hadron_2ptFunction_eo_Wilson();

  Test_Spectrum_Wilson::hadron_2ptFunction();

#ifdef USE_MPI
  // these tests run only in single-node environment.
#else
  Test_Spectrum_CRSMatrix::clover_lex();
#endif
#endif

#ifdef USE_TEST_TOPOLOGICALCHARGE
  Test_TopologicalCharge::topological_charge();
#endif

#ifdef USE_TEST_WILSONLOOP
  Test_WilsonLoop::wilsonloop();
#endif

//- #endif of #if defined(USE_GROUP_SU2)
#endif

  return EXIT_SUCCESS;
}
