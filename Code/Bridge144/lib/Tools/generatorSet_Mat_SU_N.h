/*!
        @file    $Id:: generatorSet_Mat_SU_N.h #$

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 18:35:38 #$

        @version $LastChangedRevision: 1571 $
*/

#ifndef GENERATORSET_MAT_SU_N_INCLUDED
#define GENERATORSET_MAT_SU_N_INCLUDED

#include <cassert>

#include "Parameters/commonParameters.h"
#include "Communicator/communicator.h"

#include "mat_SU_N.h"
using namespace SU_N;

#include "IO/bridgeIO.h"
using Bridge::vout;

//! Set of SU(N) generators.

/*!
   This class defines generators of SU(Nc) group.
   General value of Nc is accepted.
   The normalization of generator is tr[T_a T_b] = \delta_{ab}/2.
                                       [12 Feb 2013 H.Matsufuru]
 */
class GeneratorSet_Mat_SU_N
{
  int m_Nc;
  int m_NcA;
  std::vector<Mat_SU_N *> m_Ta;
  Bridge::VerboseLevel    m_vl;

 public:
  GeneratorSet_Mat_SU_N() : m_vl(CommonParameters::Vlevel())
  {
    int Nc = CommonParameters::Nc();

    setup(Nc);
  }

  GeneratorSet_Mat_SU_N(int Nc) : m_vl(CommonParameters::Vlevel())
  {
    setup(Nc);
  }

  ~GeneratorSet_Mat_SU_N()
  {
    tidyup();
  }

  static const std::string class_name;

  void setup(const int Nc);

  void tidyup();

  Mat_SU_N get_generator(const int ica)
  {
    assert(ica < m_NcA);
    return *m_Ta[ica];
  }

  void print();
};
#endif
