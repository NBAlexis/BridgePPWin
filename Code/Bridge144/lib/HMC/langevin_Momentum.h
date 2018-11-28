/*!
        @file    $Id:: langevin_Momentum.h #$

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 18:35:38 #$

        @version $LastChangedRevision: 1571 $
*/


#ifndef LANGEVIN_MOMENTUM_INCLUDED
#define LANGEVIN_MOMENTUM_INCLUDED

#include "Field/field_G.h"
#include "Parameters/parameters.h"
#include "Tools/generatorSet_Mat_SU_N.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

class RandomNumbers;

//! Langevin part of HMC for conjugate momentum to link variable.

/*!
   This class implements the Langevin part of the standartd HMC
   for conjugate momentum to the link variable.
   The conjugate momentum field is Hermitian and traceless,
   and set to random values with Gaussian distribution.
   This class was separated from HMC_Leapfrog (or equivalently
   HMC_General) so as to implement for general value of Nc.
                                     [13 Feb 2013 H.Matsufuru]
*/
class Langevin_Momentum
{
 private:
  RandomNumbers        *m_rand;
  Bridge::VerboseLevel m_vl;

 public:
  //! Constructor requires a pointer to random number generator.
  Langevin_Momentum(RandomNumbers *rand) :
    m_vl(CommonParameters::Vlevel())
  {
    m_rand = rand;
  }

  ~Langevin_Momentum() {}

 private:
  // non-copyable
  Langevin_Momentum(const Langevin_Momentum&);
  Langevin_Momentum& operator=(const Langevin_Momentum&);

 public:
  //! Setting conjugate momenta and returns kinetic part of Hamiltonian.
  double set_iP(Field_G& iP);

 private:
  //! Implementation for SU(3)
  double set_iP_SU3(Field_G& iP);

  //! Implementation for general value of Nc for SU(Nc) link variables.
  double set_iP_general_SU_N(Field_G& iP);

  //! Alternative of set_iP_SU3() for checking set_iP_general_SU_N().
  double set_iP_SU3_alt(Field_G& iP);
};
#endif
