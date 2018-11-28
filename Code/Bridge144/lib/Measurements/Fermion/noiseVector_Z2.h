/*!
        @file    $Id:: noiseVector_Z2.h #$

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 00:49:46 #$

        @version $LastChangedRevision: 1561 $
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

class NoiseVector_Z2 : public NoiseVector
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
