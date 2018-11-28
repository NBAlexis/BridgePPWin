/*!
        @file    $Id:: noiseVector.h #$

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 00:49:46 #$

        @version $LastChangedRevision: 1561 $
*/


#ifndef NOISEVECTOR_INCLUDED
#define NOISEVECTOR_INCLUDED

#include "Field/field.h"

#include "IO/bridgeIO.h"

//! Base class for noise vector generator.

/*!
    This is the base class of noise vector generator for
    trace calculations.
    This class only defines the interface.
                                     [2 Sep 2012 H.Matsufuru]
 */

class NoiseVector
{
 protected:
  Bridge::VerboseLevel m_vl;

 public:

  NoiseVector()
    : m_vl(CommonParameters::Vlevel()) {}

  virtual ~NoiseVector() {}

 private:
  // non-copyable
  NoiseVector(const NoiseVector&);
  NoiseVector& operator=(const NoiseVector&);

 public:

  void set_parameter_verboselevel(const Bridge::VerboseLevel vl) { m_vl = vl; }

  //! setting a noise vector.
  virtual void set(Field& v) = 0;
};
#endif
