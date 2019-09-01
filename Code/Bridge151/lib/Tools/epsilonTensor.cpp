/*!
        @file    epsilonTensor.cpp

        @brief

        @author  Yusuke Namekawa (namekawa)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/
#include "BridgeLib_Private.h"
#include "epsilonTensor.h"
#include <assert.h>

// const std::string Parameters::class_name = "EpsilonTensor";

//====================================================================
void EpsilonTensor::init()
{
  m_Nepsilon_3       = 3;
  m_Nepsilon_3_index = 6;

  m_epsilon_3_index.resize(m_Nepsilon_3_index);
  for (int n = 0; n < m_Nepsilon_3_index; ++n) {
    m_epsilon_3_index[n].resize(m_Nepsilon_3);
  }
  m_epsilon_3_index[0][0] = 0;
  m_epsilon_3_index[0][1] = 1;
  m_epsilon_3_index[0][2] = 2;

  m_epsilon_3_index[1][0] = 1;
  m_epsilon_3_index[1][1] = 2;
  m_epsilon_3_index[1][2] = 0;

  m_epsilon_3_index[2][0] = 2;
  m_epsilon_3_index[2][1] = 0;
  m_epsilon_3_index[2][2] = 1;

  m_epsilon_3_index[3][0] = 2;
  m_epsilon_3_index[3][1] = 1;
  m_epsilon_3_index[3][2] = 0;

  m_epsilon_3_index[4][0] = 1;
  m_epsilon_3_index[4][1] = 0;
  m_epsilon_3_index[4][2] = 2;

  m_epsilon_3_index[5][0] = 0;
  m_epsilon_3_index[5][1] = 2;
  m_epsilon_3_index[5][2] = 1;

  m_epsilon_3_value.resize(m_Nepsilon_3_index);
  m_epsilon_3_value[0] = 1;
  m_epsilon_3_value[1] = 1;
  m_epsilon_3_value[2] = 1;
  m_epsilon_3_value[3] = -1;
  m_epsilon_3_value[4] = -1;
  m_epsilon_3_value[5] = -1;
}


//====================================================================
int EpsilonTensor::epsilon_3_index(const int n, const int i) const
{
  //- range check
  assert(n >= 0);
  assert(n < m_Nepsilon_3_index);

  assert(i >= 0);
  assert(i < m_Nepsilon_3);

  return m_epsilon_3_index[n][i];
}


//====================================================================
int EpsilonTensor::epsilon_3_value(const int n) const
{
  //- range check
  assert(n >= 0);
  assert(n < m_Nepsilon_3_index);

  return m_epsilon_3_value[n];
}


//====================================================================
//============================================================END=====
