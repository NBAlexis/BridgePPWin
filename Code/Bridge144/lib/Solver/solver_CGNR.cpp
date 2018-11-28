#include "BridgeLib_Private.h"

/*!
        @file    $Id: solver_CGNR.cpp #$

        @brief

        @author  Tatsumi Aoyama  (aoym)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate: 2013-04-27 12:28:50 #$

        @version $LastChangedRevision: 1561 $
*/

#include "solver_CGNR.h"

using Bridge::vout;


#ifdef USE_FACTORY
namespace {
  Solver *create_object(Fopr *fopr)
  {
    return new Solver_CGNR(fopr);
  }


  bool init = Solver::Factory::Register("CGNR", create_object);
}
#endif


const std::string Solver_CGNR::class_name = "Solver_CGNR";

//====================================================================
void Solver_CGNR::set_parameters(const Parameters& params)
{
  return this->Solver_CG::set_parameters(params);
}


//====================================================================
void Solver_CGNR::solve(Field& xq, const Field& b, int& Nconv, double& diff)
{
  Fopr   *fopr     = this->Solver_CG::get_fopr();
  string mode_prev = fopr->get_mode();

  vout.detailed(m_vl, "%s: solver starts\n", class_name.c_str());

  if (mode_prev == "DdagD") {
    this->Solver_CG::solve(xq, b, Nconv, diff);  // fallback to CG solver
    return;
  }

  if (!((mode_prev == "D") || (mode_prev == "Ddag"))) {
    vout.crucial(m_vl, "Error at %s: unsupported mode for fopr %s.", class_name.c_str(), mode_prev.c_str());
    exit(EXIT_FAILURE);
  }

#pragma omp barrier
#pragma omp master
  {
    if ((m_b2.nin() != b.nin()) || (m_b2.nvol() != b.nvol()) ||
        (m_b2.nex() != b.nex())) {
      m_b2.reset(b.nin(), b.nvol(), b.nex());
    }
  }
#pragma omp barrier

  fopr->mult_dag(m_b2, b);  // b2 = fopr->mult_dag(b);

#pragma omp barrier
#pragma omp master
  {
    if (mode_prev == "D") {
      fopr->set_mode("DdagD");
    } else if (mode_prev == "Ddag") {
      fopr->set_mode("DDdag");
    }
  }
#pragma omp barrier

  this->Solver_CG::solve(xq, m_b2, Nconv, diff);

#pragma omp barrier
#pragma omp master
  {
    fopr->set_mode(mode_prev);
  }
#pragma omp barrier
}


//====================================================================
double Solver_CGNR::flop_count()
{
  //int    NPE = CommonParameters::NPE();
  double eps = CommonParameters::epsilon_criterion();

  double flop_fopr = this->Solver_CG::get_fopr()->flop_count();

  if (flop_fopr < eps) {
    vout.crucial(m_vl, "Warning at %s: no fopr->flop_count() is available, setting flop = 0.0.\n", class_name.c_str());
    return 0.0;
  }

  double flop_solver = this->Solver_CG::flop_count();

  double flop = flop_solver + flop_fopr;

  return flop;
}


//====================================================================
//============================================================END=====
