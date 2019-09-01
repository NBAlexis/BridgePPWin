/*!
        @file    fprop_Standard_lex.cpp

        @brief

        @author  Satoru Ueda  (ueda)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/
#include "BridgeLib_Private.h"
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
