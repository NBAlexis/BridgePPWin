/*!
        @file    randomNumbers_Mseries.h

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
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
    This class BAPI generates the M-series random numbers.
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

class BAPI RandomNumbers_Mseries : public RandomNumbers
{
  //  static const double Fnorm = 4.656612870908988e-10;
  static const double Fnorm;  //!< initialized in .cpp file.

 public:
  static const std::string class_name;

 private:
  static const int Np = 521;
  static const int Nq = 32;
  int w[Np];
  int jr, kr;

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

#ifdef USE_FACTORY
 private:
  static RandomNumbers *create_object_with_int(const int& iseed)
  {
    return new RandomNumbers_Mseries(iseed);
  }

  static RandomNumbers *create_object_with_file(const std::string& filename)
  {
    return new RandomNumbers_Mseries(filename);
  }

 public:
  static bool register_factory()
  {
    bool init1 = RandomNumbers::Factory_int::Register("Mseries", create_object_with_int);
    bool init2 = RandomNumbers::Factory_file::Register("Mseries", create_object_with_file);

    return init1 && init2;
  }
#endif
};
#endif
