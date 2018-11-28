#include "BridgeLib_Private.h"

/*!
        @file    $Id:: fprop_Standard_lex.cpp #$

        @brief

        @author  Satoru Ueda  (ueda)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 00:49:46 #$

        @version $LastChangedRevision: 1561 $
*/

#include "fprop_Standard_lex.h"

const std::string Fprop_Standard_lex::class_name = "Fprop_Standard_lex";

//====================================================================
void Fprop_Standard_lex::set_config(Field *U)
{
  m_solver->get_fopr()->set_config(U);
}


//====================================================================
void Fprop_Standard_lex::invert_D(Field& xq, const Field& b, int& Nconv, double& diff)
{
  m_solver->get_fopr()->set_mode("D");

#pragma omp parallel
  {
    m_solver->solve(xq, b, Nconv, diff);
  }
}


//====================================================================
void Fprop_Standard_lex::invert_DdagD(Field& xq, const Field& b, int& Nconv, double& diff)
{
  m_solver->get_fopr()->set_mode("DdagD");

#pragma omp parallel
  {
    m_solver->solve(xq, b, Nconv, diff);
  }
}


//====================================================================
double Fprop_Standard_lex::flop_count()
{
  return m_solver->flop_count();
}


//====================================================================
//============================================================END=====
