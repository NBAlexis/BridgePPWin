/*!
        @file    solver.cpp

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/
#include "BridgeLib_Private.h"
#include "solver.h"

#ifdef USE_FACTORY

#ifdef USE_FACTORY_AUTOREGISTER
#else
// setup factories for all subclasses

#include "solver_CG.h"
#include "solver_CGNE.h"
#include "solver_CGNR.h"
#include "solver_BiCGStab_Cmplx.h"
#include "solver_BiCGStab_L_Cmplx.h"
#include "solver_BiCGStab_DS_L_Cmplx.h"
#include "solver_BiCGStab_IDS_L_Cmplx.h"
#include "solver_GMRES_m_Cmplx.h"

bool Solver::init_factory()
{
  bool result = true;

  result &= Solver_CG::register_factory();
  result &= Solver_CGNE::register_factory();
  result &= Solver_CGNR::register_factory();
  result &= Solver_BiCGStab_Cmplx::register_factory();
  result &= Solver_BiCGStab_L_Cmplx::register_factory();
  result &= Solver_BiCGStab_DS_L_Cmplx::register_factory();
  result &= Solver_BiCGStab_IDS_L_Cmplx::register_factory();
  result &= Solver_GMRES_m_Cmplx::register_factory();

  return result;
}


#endif /* USE_FACTORY_AUTOREGISTER */

#endif /* USE_FACTORY */
