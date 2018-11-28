//=============================================================================
// FILENAME : BridgeLib_Private.h
// 
// DESCRIPTION:
//
// REVISION:
//  [11/28/2018 nbale]
//=============================================================================


#pragma once

#ifndef _BRIDGELIB_PRIVATE_H_
#define _BRIDGELIB_PRIVATE_H_

#include "configure.h"
#include "defs.h"
#include "BridgeLibSetup.h"
#include "bridge_complex.h"

//First few defines
#include "Communicator/communicator.h"
#include "IO/bridgeIO.h"
#include "Parameters/commonParameters.h"

//=======================================================
//Communicator
#if COMMUNICATOR_USE_SINGLE
#include "Communicator//Single/channel.h"
#elif COMMUNICATOR_USE_MPI
#include "Communicator/MPI/channel.h"
#include "Communicator/MPI/layout.h"
#include "Communicator/MPI/communicator_mpi.h"
#endif

//=======================================================
//ResourceManager
#include "ResourceManager/threadManager_OpenMP.h"

//=======================================================
//Parameters
#include "Parameters/tinyxml2.h"
#include "Parameters/checker.h"
#include "Parameters/parameters.h"
#include "Parameters/parameterManager.h"
#include "Parameters/parameterManager_XML.h"
#include "Parameters/parameterManager_YAML.h"

//=======================================================
//Tools

//hh files
#include "Tools/location.hh"
#include "Tools/position.hh"
#include "Tools/stack.hh"

//other
#include "Tools/sorter.h"
#include "Tools/timer.h"
#include "Tools/factory.h"
#include "Tools/file_utils.h"
#include "Tools/filename.h" //this file is only use in test...

//randoms
#include "Tools/randomNumberManager.h"
#include "Tools/randomNumbers.h"
#include "Tools/randomNumbers_Mseries.h"
#include "Tools/randomNumbers_MT19937.h"

//math
#include "Tools/math_Sign_Zolotarev.h"
#include "Tools/math_Rational.h"

//gamma matrix
#include "Tools/gammaMatrix.h"
#include "Tools/gammaMatrixSet.h"
#include "Tools/gammaMatrixSet_Chiral.h"
#include "Tools/gammaMatrixSet_Dirac.h"

//SU(N)
#include "Tools/decompose_QR_Cmplx.h"
#include "Tools/decompose_LU_Cmplx.h"
#include "Tools/decompose_LUP_Cmplx.h"
#include "Tools/decompose_Hessenberg_Cmplx.h"
#include "Tools/eigen_QR_Cmplx.h"
#include "Tools/mat_SU_N.h"
#include "Tools/vec_SU_N.h"
#include "Tools/generatorSet_Mat_SU_N.h"

//the eval expr parser looks like not been used. In fact, we have GINAC working on Windows now, may replace to GINAC though

//=======================================================
//=======================================================
//Real Physics start here
//=======================================================
//=======================================================

//=======================================================
//Fields
#include "Field/field.h"
#include "Field/field_G.h"
#include "Field/field_F.h"
#include "Field/field_F_SF.h"
#include "Field/field_G_SF.h"
#include "Field/index_lex.h" //lex for lexicographical, (what is lexicographical)?
#include "Field/index_eo.h" //eo for even odd
#include "Field/shiftField_lex.h"
#include "Field/shiftField_eo.h"

//now we can include director, it need field_G and it is used in Fopr
#include "Tools/director.h"

//=======================================================
//Fopr
#include "Fopr/fopr.h"
#include "Fopr/fopr_eo.h"

#include "Fopr/fopr_Wilson.h"
#include "Fopr/fopr_WilsonGeneral.h"
#include "Fopr/fopr_Wilson_eo.h"
#include "Fopr/fopr_Wilson_Isochemical.h"
#include "Fopr/fopr_Wilson_SF.h"

#include "Fopr/fopr_CloverTerm.h"
#include "Fopr/fopr_CloverTerm_General.h"
#include "Fopr/fopr_CloverTerm_eo.h"
#include "Fopr/fopr_Clover.h"
#include "Fopr/fopr_CloverGeneral.h"
#include "Fopr/fopr_Clover_eo.h"
#include "Fopr/fopr_Clover_Isochemical.h"
#include "Fopr/fopr_Clover_SF.h"

#include "Fopr/fopr_Chebyshev.h"

#include "Fopr/fopr_CRS.h" //Parallel version is not implemented.

//Fopr rational need solver rational
//Fopr smear need smear

//=======================================================
//Solver
#include "Solver/solver.h"
#include "Solver/shiftsolver.h"
#include "Solver/shiftsolver_CG.h"
#include "Solver/solver_CG.h"
#include "Solver/solver_CGNE.h"
#include "Solver/solver_CGNR.h"
#include "Solver/solver_BiCGStab_Cmplx.h"
#include "Solver/solver_BiCGStab_L_Cmplx.h"
#include "Solver/solver_BiCGStab_DS_L_Cmplx.h"
#include "Solver/solver_BiCGStab_IDS_L_Cmplx.h"

//now add Fopr rational
#include "Fopr/fopr_Rational.h"
#include "Fopr/fopr_Rational_SF.h"

//=======================================================
//Measurements
#include "Measurements/Gauge/staple.h"
#include "Measurements/Gauge/staple_lex.h"
#include "Measurements/Gauge/staple_eo.h"
#include "Measurements/Gauge/staple_SF.h"
#include "Measurements/Gauge/wilsonLoop.h"
#include "Measurements/Gauge/polyakovLoop.h"
#include "Measurements/Gauge/fieldStrength.h"
#include "Measurements/Gauge/energyDensity.h"
#include "Measurements/Gauge/topologicalCharge.h"
#include "Measurements/Gauge/gaugeFixing.h"
#include "Measurements/Gauge/gaugeFixing_Coulomb.h"
#include "Measurements/Gauge/gaugeFixing_Landau.h"
#include "Measurements/Gauge/gaugeFixing_None.h"
#include "Measurements/Gauge/gradientFlow_RungeKutta.h"
#include "Measurements/Gauge/gradientFlow_RungeKutta_1st.h"
#include "Measurements/Gauge/gradientFlow_RungeKutta_2nd.h"
#include "Measurements/Gauge/gradientFlow_RungeKutta_3rd.h"
#include "Measurements/Gauge/gradientFlow_RungeKutta_4th.h"
#include "Measurements/Gauge/gradientFlow_AdaptiveRungeKutta.h"
#include "Measurements/Gauge/gradientFlow.h"

#include "Measurements/Fermion/contract_4spinor.h"
#include "Measurements/Fermion/corr2pt_4spinor.h"
#include "Measurements/Fermion/corr2pt_Wilson_SF.h"
#include "Measurements/Fermion/fprop.h"
#include "Measurements/Fermion/fprop_Standard_eo.h"
#include "Measurements/Fermion/fprop_Standard_lex.h"
#include "Measurements/Fermion/fprop_Wilson_Shift.h"
#include "Measurements/Fermion/source.h"
#include "Measurements/Fermion/source_Local.h"
#include "Measurements/Fermion/source_Wall.h"
//#include "Measurements/Fermion/source_Wall_SF.h" need smear
#include "Measurements/Fermion/source_MomentumWall.h"

#include "Measurements/Fermion/noiseVector.h"
#include "Measurements/Fermion/noiseVector_Z2.h"
#include "Measurements/Fermion/quarkNumberSusceptibility_Wilson.h"

//=======================================================
//Force
#include "Force/Gauge/force_G.h"
#include "Force/Gauge/force_G_Plaq.h"
#include "Force/Gauge/force_G_Plaq_SF.h"
#include "Force/Gauge/force_G_Rectangle.h"
#include "Force/Gauge/force_G_Rectangle_SF.h"

//=======================================================
//Action
#include "Action/action.h"
#include "Action/Gauge/action_G_Plaq.h"
#include "Action/Gauge/action_G_Rectangle.h"
#include "Action/Gauge/action_G_Plaq_SF.h"
#include "Action/Gauge/action_G_Rectangle_SF.h"
#include "Action/Fermion/action_F_Standard_lex.h"
#include "Action/Fermion/action_F_Standard_eo.h"
#include "Action/Fermion/action_F_Standard_SF.h"
#include "Action/Fermion/action_F_Rational.h"
#include "Action/Fermion/action_F_Rational_SF.h"
#include "Action/Fermion/action_F_Ratio_lex.h"
#include "Action/Fermion/action_F_Ratio_eo.h"


//=======================================================
//IO


#endif //#ifndef _BRIDGELIB_PRIVATE_H_
//=============================================================================
// END OF FILE
//=============================================================================