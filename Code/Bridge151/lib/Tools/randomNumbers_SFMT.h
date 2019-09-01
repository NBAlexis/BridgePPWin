/*!
        @file    randomNumbers_SFMT.h

        @brief

        @author  T.Aoyama  (aoym)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/

#ifndef RANDOMNUMBERS_SFMT_INCLUDED
#define RANDOMNUMBERS_SFMT_INCLUDED

#ifdef USE_SFMTLIB

// #define __STDC_CONSTANT_MACROS

#include <string>
#include <cassert>

#include "randomNumbers.h"

#include <SFMT.h>

#ifdef USE_NTLLIB
#define ENABLE_SFMT_JUMP
#include <SFMT-jump-alt.h>
#endif

#include "IO/bridgeIO.h"
using Bridge::vout;

class BAPI RandomNumbers_SFMT : public RandomNumbers
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

#ifdef USE_FACTORY
 private:
  static RandomNumbers *create_object_with_int(const int& iseed)
  {
    return new RandomNumbers_SFMT(iseed);
  }

  static RandomNumbers *create_object_with_file(const std::string& filename)
  {
    return new RandomNumbers_SFMT(filename);
  }

 public:
  static bool register_factory()
  {
    bool init1 = RandomNumbers::Factory_int::Register("SFMT", create_object_with_int);
    bool init2 = RandomNumbers::Factory_file::Register("SFMT", create_object_with_file);

    return init1 && init2;
  }
#endif
};
#endif /* USE_SFMTLIB */
#endif /* RANDOMNUMBERS_SFMT_INCLUDED */
