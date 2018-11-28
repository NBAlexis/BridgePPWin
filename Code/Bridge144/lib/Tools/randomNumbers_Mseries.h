/*!
        @file    $Id:: randomNumbers_Mseries.h #$

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 18:35:38 #$

        @version $LastChangedRevision: 1571 $
*/

#ifndef RANDOMNUMBERS_MSERIES_INCLUDED
#define RANDOMNUMBERS_MSERIES_INCLUDED

#include <assert.h>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>

#include "randomNumbers.h"
#include "Communicator/communicator.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//! Random number generator base on M-series.

/*!
    This class generates the M-series random numbers.
    The original version in Fortran was written by
             J.Makino and O.Miyamura (Ver.3.0 21 July 1991).
    Public version is available under GNU GPL:
     Shinji Hioki, QCDMPI http://insam.sci.hiroshima-u.ac.jp/QCDMPI/
    which implements
     Jun Makino, "Lagged-Fibonacci random number generators on parallel
     computers", Parallel Computing, 20 (1994) 1357-1367.

    An instance is created with a given integer number which is used
    to set the initial random numbers.
                                          [23 Jul 2012 H.Matsufuru]
 */

class RandomNumbers_Mseries : public RandomNumbers
{
  //  static const double Fnorm = 4.656612870908988e-10;
  static const double Fnorm;  //!< initialized in .cpp file.

 public:
  static const std::string class_name;

 private:
  static const int Np = 521, Nq = 32;
  int              w[Np];
  int              jr, kr;

  double sq2r;
  double pi, pi2;

 public:
  RandomNumbers_Mseries(const int ndelay)
  {
    initset(ndelay);
  }

  RandomNumbers_Mseries(const std::string& filename)
  {
    read_file(filename);
  }

  double get()
  {
    w[jr] = w[jr] ^ w[kr];
    double rw = w[jr] * Fnorm;
    jr = jr + 1;
    if (jr >= Np) jr = jr - Np;
    kr = kr + 1;
    if (kr >= Np) kr = kr - Np;
    return rw;
  }

  void write_file(const std::string&);
  void read_file(const std::string&);

  void reset(unsigned long seed);

 private:

  void initset(const int ndelay);

  void delay3(const int ndelay);
};
#endif
