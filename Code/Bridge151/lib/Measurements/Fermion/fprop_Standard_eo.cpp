/*!
        @file    fprop_Standard_eo.cpp

        @brief

        @author  Satoru Ueda (sueda)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/
#include "BridgeLib_Private.h"
#include "fprop_Standard_eo.h"

const std::string Fprop_Standard_eo::class_name = "Fprop_Standard_eo";

//====================================================================
void Fprop_Standard_eo::set_config(Field *U)
{
  m_index->convertField(*m_Ueo, *U);

  m_solver->get_fopr()->set_config(U);
}


//====================================================================
void Fprop_Standard_eo::invert_D(Field& xq, const Field& b, int& Nconv, double& diff)
{
  const int Nin  = b.nin();
  const int Nvol = b.nvol();
  const int Nex  = b.nex();

  Field Be(Nin, Nvol / 2, Nex);
  Field bo(Nin, Nvol / 2, Nex);
  Field xe(Nin, Nvol / 2, Nex);

  int Nconv1 = 0;

  Fopr_eo *fopr = (Fopr_eo *)m_solver->get_fopr();

  fopr->set_mode("D");
  fopr->preProp(Be, bo, b);

#pragma omp parallel
  {
    m_solver->solve(xe, Be, Nconv1, diff);
  }

  fopr->postProp(xq, xe, bo);

  //- NB. #mult is doubled for even-odd
  Nconv = 2 * Nconv1;
}


//====================================================================
void Fprop_Standard_eo::invert_DdagD(Field& xq, const Field& b, int& Nconv, double& diff)
{
  const int Nin  = b.nin();
  const int Nvol = b.nvol();
  const int Nex  = b.nex();

  Field Be(Nin, Nvol / 2, Nex);
  Field bo(Nin, Nvol / 2, Nex);
  Field xe(Nin, Nvol / 2, Nex);

  int    Nconv1 = 0, Nconv2 = 0;
  double diff1 = 1.0, diff2 = 1.0;

  Fopr_eo *fopr = (Fopr_eo *)m_solver->get_fopr();

  fopr->set_mode("Ddag");
  fopr->preProp(Be, bo, b);

#pragma omp parallel
  {
    m_solver->solve(xe, Be, Nconv1, diff1);
  }

  fopr->postProp(xq, xe, bo);

  fopr->set_mode("D");
  fopr->preProp(Be, bo, xq);

#pragma omp parallel
  {
    m_solver->solve(xe, Be, Nconv2, diff2);
  }

  fopr->postProp(xq, xe, bo);

  //- NB. #mult is doubled for even-odd
  Nconv = 2 * (Nconv1 + Nconv2);

  //- rough estimate of diff
  diff = (diff1 + diff2) / 2.0;
}


//====================================================================
double Fprop_Standard_eo::flop_count()
{
  const int    NPE = CommonParameters::NPE();
  const double eps = CommonParameters::epsilon_criterion();

  //- NB1 Nin = 2 * Nc * Nd, Nex = 1  for field_F
  //- NB2 Nvol/2 for eo
  const int Nin  = 2 * CommonParameters::Nc() * CommonParameters::Nd();
  const int Nvol = CommonParameters::Nvol();
  const int Nex  = 1;

  const double flop_fopr = m_solver->get_fopr()->flop_count();

  if (flop_fopr < eps) {
    vout.crucial(m_vl, "Warning at %s: no fopr->flop_count() is available, setting flop = 0.0.\n", class_name.c_str());
    return 0.0;
  }

  const double flop_axpy = static_cast<double>(Nin * Nex * 2) * (Nvol / 2 * NPE);

  const double flop_preProp  = flop_fopr + flop_axpy;
  const double flop_solver   = m_solver->flop_count();
  const double flop_postProp = flop_fopr + flop_axpy;

  const double flop = flop_preProp + 2 * flop_solver + flop_postProp;

  return flop;
}


//====================================================================
//============================================================END=====
