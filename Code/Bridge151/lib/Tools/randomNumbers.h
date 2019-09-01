/*!
        @file    randomNumbers.h

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/

#ifndef RANDOMNUMBERS_INCLUDED
#define RANDOMNUMBERS_INCLUDED

#include "Field/field.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

#ifdef USE_FACTORY
#include "Tools/factory.h"
#endif


//! Base class BAPI of random number generators.

/*!
   This class BAPI defines the interface of random number
   generator, and implements common methods.
   Practical methods to generate random numbers are
   defined in subclasses.
   This class BAPI also implements Gaussian random number and
   method to set a global field of Gaussian random numbers
   and cut it out to the local field for the own node
   (gauss_lex_global()) which is useful in HMC etc.
                                 [25 Dec 2011 H.Matsufuru]
   U1,Z2 are added               [11 Jan 2017 Y.Namekawa]
   Factory is introduced         [ 2 Feb 2017 Y.Namekawa]
 */

class BAPI RandomNumbers
{
 public:
  static const std::string class_name;

 protected:
  Bridge::VerboseLevel m_vl;

 public:

  RandomNumbers()
    : m_vl(CommonParameters::Vlevel()) {}

  virtual ~RandomNumbers() {}

 private:
  // non-copyable
  RandomNumbers(const RandomNumbers&);
  RandomNumbers& operator=(const RandomNumbers&);

 public:
  void set_parameter_verboselevel(const Bridge::VerboseLevel vl) { m_vl = vl; }

  virtual double get() = 0;

  void gauss(double& rand1, double& rand2);


  virtual void lex_global(const std::string&, Field&);

  //! gaussian random number defined on global lattice.
  virtual void gauss_lex_global(Field&);

  //! uniform random number defined on global lattice.
  virtual void uniform_lex_global(Field&);

  //! U(1) random number defined on global lattice.
  virtual void U1_lex_global(Field&);

  //! Z(2) random number defined on global lattice.
  virtual void Z2_lex_global(Field&);


  //! gaussian noise for even-odd perconditioned field (S.UEDA)
  virtual void gauss_eo_global(Field&);


  //! save and load random number status.
  virtual void read_file(const std::string&)  = 0;
  virtual void write_file(const std::string&) = 0;

  //! reset state with new seed.
  virtual void reset(unsigned long seed) = 0;


 protected:
  class BAPI rand_gauss_even
  {
   public:
    rand_gauss_even(Field& f, RandomNumbers *rand_gauss)
      : m_rand_gauss(rand_gauss), m_ptr(f.ptr(0)), m_block(f.nin()) {}
    void operator()(const bool do_fill);
    size_t block_size() const;

   private:
    RandomNumbers *m_rand_gauss;
    double *m_ptr;
    size_t m_block;
  };

  class BAPI rand_gauss_odd
  {
   public:
    rand_gauss_odd(Field& f, RandomNumbers *rand_gauss)
      : m_rand_gauss(rand_gauss), m_ptr(f.ptr(0)), m_block(f.nin()) {}
    void operator()(const bool do_fill);
    size_t block_size() const;

   private:
    RandomNumbers *m_rand_gauss;
    double *m_ptr;
    size_t m_block;
  };

  class BAPI rand_uniform
  {
   public:
    rand_uniform(Field& f, RandomNumbers *rand_gauss)
      : m_rand_gauss(rand_gauss), m_ptr(f.ptr(0)), m_block(f.nin()) {}
    void operator()(const bool do_fill);
    size_t block_size() const;

   private:
    RandomNumbers *m_rand_gauss;
    double *m_ptr;
    size_t m_block;
  };

 private:
  template<typename InnerGenerator>
  void generate_global(Field& f);

#ifdef USE_FACTORY
 public:
  typedef RandomNumbers *(*ProductCreator_int)(const int& iseed);
  typedef RandomNumbers *(*ProductCreator_file)(const std::string& filename);

  typedef FactoryTemplate<RandomNumbers, ProductCreator_int>    Factory_int;
  typedef FactoryTemplate<RandomNumbers, ProductCreator_file>   Factory_file;

  static RandomNumbers *New(const IdentifierType& subtype, const int& iseed)
  {
    ProductCreator_int p = Factory_int::Find(subtype);

    return p ? (*p)(iseed) : 0;
  }

  static RandomNumbers *New(const IdentifierType& subtype, const std::string& filename)
  {
    ProductCreator_file p = Factory_file::Find(subtype);

    return p ? (*p)(filename) : 0;
  }

#ifdef USE_FACTORY_AUTOREGISTER
#else
  static bool init_factory();
#endif
#endif
};
#endif
