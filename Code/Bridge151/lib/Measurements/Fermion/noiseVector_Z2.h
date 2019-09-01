/*!
        @file    noiseVector_Z2.h

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/

#ifndef NOISEVECTOR_Z2_INCLUDED
#define NOISEVECTOR_Z2_INCLUDED

#include "noiseVector.h"

#include "Tools/randomNumbers.h"

//! Z2 Noise vector for a trace calculation.

/*!
    Z2 Noise vector.
                                     [30 Aug 2012 H.Matsufuru]
 */

class BAPI NoiseVector_Z2 : public NoiseVector
{
 public:
  static const std::string class_name;

 private:
  RandomNumbers *m_rand;

 public:
  NoiseVector_Z2(RandomNumbers *rand)
    : NoiseVector(), m_rand(rand) {}

  NoiseVector_Z2(unique_ptr<RandomNumbers>& rand)
    : NoiseVector(), m_rand(rand.get()) {}

  void set(Field& v);
};
#endif
