#include "BridgeLib_Private.h"

/*!
        @file    $Id: solver_CGNE.cpp #$

        @brief

        @author  Tatsumi Aoyama  (aoym)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate: 2013-04-27 12:28:50 #$

        @version $LastChangedRevision: 1561 $
*/

#include "solver_CGNE.h"

using Bridge::vout;


#ifdef USE_FACTORY
namespace {
  Solver *create_object(Fopr *fopr)
  {
    return new Solver_CGNE(fopr);
  }


  bool init = Solver::Factory::Register("CGNE", create_object);
}
#endif


const std::string Solver_CGNE::class_name = "Solver_CGNE";

//====================================================================
void Solver_CGNE::set_parameters(const Parameters& params)
{
  return this->Solver_CG::set_parameters(params);
}


//====================================================================
void Solver_CGNE::solve(Field& xq, const Field& b, int& Nconv, double& diff)
{
  Fopr   *fopr     = this->Solver_CG::get_fopr();
  string mode_prev = fopr->get_mode();

  vout.detailed(m_vl, "%s: solver starts\n", class_name.c_str());

  if (mode_prev == "DdagD") {
    this->Solver_CG::solve(xq, b, Nconv, diff);  // fallback to CG solver
    return;
  }

  if (!((mode_prev == "D") || (mode_prev == "Ddag"))) {
    vout.crucial(m_vl, "Error at %s: unsupported mode for fopr %s.",
                 class_name.c_str(), mode_prev.c_str());
    exit(EXIT_FAILURE);
  }

#pragma omp barrier
#pragma omp master
  {
    if ((m_y.nin() != b.nin()) || (m_y.nvol() != b.nvol()) ||
        (m_y.nex() != b.nex())) {
      m_y.reset(b.nin(), b.nvol(), b.nex());
    }

    if (mode_prev == "D") {
      fopr->set_mode("DDdag");
    } else if (mode_prev == "Ddag") {
      fopr->set_mode("DdagD");
    }
  }
#pragma omp barrier

  this->Solver_CG::solve(m_y, b, Nconv, diff);

#pragma omp barrier
#pragma omp master
  {
    fopr->set_mode(mode_prev);
  }
#pragma omp barrier

  fopr->mult_dag(xq, m_y);  // xq = fopr->mult_dag(y);

#pragma omp barrier
}


//====================================================================
double Solver_CGNE::flop_count()
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
