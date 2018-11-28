/*!
        @file    $Id:: randomNumbers_SFMT.h #$

        @brief

        @author  T.Aoyama  (aoym)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 00:49:46 #$

        @version $LastChangedRevision: 1561 $
*/

#ifndef RANDOMNUMBERS_SFMT_INCLUDED
#define RANDOMNUMBERS_SFMT_INCLUDED

#ifdef USE_SFMTLIB

// #define __STDC_CONSTANT_MACROS

#include <string>
#include <cassert>

#include "randomNumbers.h"

#include <SFMT.h>

#ifdef HAVE_NTL
#define ENABLE_SFMT_JUMP
#include <SFMT-jump-alt.h>
#endif

#include "IO/bridgeIO.h"
using Bridge::vout;

class RandomNumbers_SFMT : public RandomNumbers
{
  static const std::string class_name;

 public:
  RandomNumbers_SFMT(const int s);
  RandomNumbers_SFMT(const std::string& filename) { read_file(filename); }

  ~RandomNumbers_SFMT() {}

  double get();
  void get_block(double *v, const size_t n);

  void write_file(const std::string& filename);
  void read_file(const std::string& filename);

  void reset(unsigned long seed);

#ifdef ENABLE_SFMT_JUMP
  virtual void uniform_lex_global(Field&);
  virtual void gauss_lex_global(Field&);
#endif

 private:

  template<typename T>
  void generate_global_jump(Field& f);

  sfmt_t m_state;
};
#endif /* USE_SFMTLIB */
#endif /* RANDOMNUMBERS_SFMT_INCLUDED */
