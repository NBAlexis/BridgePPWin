/*!
        @file    randomNumbers_MT19937.h

        @brief

        @author  Satoru Ueda  (sueda)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/

#ifndef RANDOMNUMBERS_MT19937_INCLUDED
#define RANDOMNUMBERS_MT19937_INCLUDED

//==========================================================
//  Random number generator.
//     MT19937 random number generator.
//     Original version in C was written by
//     Takuji Nishimura and Makoto Matsumoto.
//     See URL http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/MT2002/mt19937ar.html
//
//                                 25 Dec 2004  H.Matsufuru
//==========================================================

#include <cmath>
#include <string>
#include <vector>
#include <assert.h>
#include <iostream>
#include <fstream>
#include <vector>
using std::vector;
#include <string>
using std::string;

#include "randomNumbers.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

class BAPI RandomNumbers_MT19937 : public RandomNumbers
{
  static const std::string class_name;

 public:
  RandomNumbers_MT19937(int s);
  //RandomNumbers_MT19937(unsigned long s = 5489UL);
  RandomNumbers_MT19937(unsigned long s);
  RandomNumbers_MT19937(std::vector<unsigned long>& key);
  RandomNumbers_MT19937(const std::string& filename);

  ~RandomNumbers_MT19937() {}

  double get() { return randDouble2(); }

  void write_file(const std::string&);
  void read_file(const std::string&);

  void reset(unsigned long seed);

 private:
  void init(unsigned long s);
  void init(unsigned long s, std::vector<unsigned long>& key);

  void nextState() const;
  unsigned long twist(unsigned long u, unsigned long v) const;

  // return [0,0xffffffff]
  unsigned long randInt32() const;

  // return [0,0x7fffffff]
  long randInt31() const;

  // return [0,1] random number
  double randDouble1() const;

  // return [0,1) random number
  double randDouble2() const;

  // return (0,1) random number
  double randDouble3() const;

  // return [0,1) random number with 53-bit resolution
  double randRes53() const;

  enum { N=624, M=397 };

  mutable int m_left;
  mutable unsigned long m_state[N];
  mutable unsigned long *m_next;

#ifdef USE_FACTORY
 private:
  static RandomNumbers *create_object_with_int(const int& iseed)
  {
    return new RandomNumbers_MT19937(iseed);
  }

  static RandomNumbers *create_object_with_file(const std::string& filename)
  {
    return new RandomNumbers_MT19937(filename);
  }

 public:
  static bool register_factory()
  {
    bool init1 = RandomNumbers::Factory_int::Register("MT19937", create_object_with_int);
    bool init2 = RandomNumbers::Factory_file::Register("MT19937", create_object_with_file);

    return init1 && init2;
  }
#endif
};
#endif
