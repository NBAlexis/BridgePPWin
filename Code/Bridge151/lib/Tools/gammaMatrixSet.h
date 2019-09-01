/*!
        @file    gammaMatrixSet.h

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/

#ifndef GAMMAMATRIXSET_INCLUDED
#define GAMMAMATRIXSET_INCLUDED

#include <cassert>

#include "Parameters/commonParameters.h"
#include "gammaMatrix.h"

#ifdef USE_FACTORY
#include "factory.h"
#endif

//! Set of Gamma Matrices: basis class BAPI.

/*!
   This class BAPI defines a set of gamma matrices.
   Present implementation is applicable only to Ndim=4.
   Just possible specied of gamma matrices are enumerated in
   this class BAPI, and practical form is given in subclass by
   implementing virtual function init_GM().
                                        [4 Feb 2012 H.Matsufuru]
 */

class BAPI GammaMatrixSet
{
 protected:
  int m_Nspecies;
  std::vector<GammaMatrix> m_gm;

  Bridge::VerboseLevel m_vl;

 public:
  enum GMspecies
  {
    UNITY, GAMMA1, GAMMA2, GAMMA3, GAMMA4, GAMMA5,
    GAMMA51, GAMMA52, GAMMA53, GAMMA54,
    GAMMA15, GAMMA25, GAMMA35, GAMMA45,
    SIGMA12, SIGMA23, SIGMA31,
    SIGMA41, SIGMA42, SIGMA43,
    CHARGECONJG
  };

  GammaMatrixSet()
    : m_vl(CommonParameters::Vlevel())
  {
    int Nd = CommonParameters::Nd();

    assert(Nd == 4);
    m_Nspecies = Nd * (Nd + 1) + 1;    // must be 21.
    m_gm.resize(m_Nspecies);
  }

  virtual ~GammaMatrixSet() {}

 private:
  // non-copyable
  GammaMatrixSet(const GammaMatrixSet&);
  GammaMatrixSet& operator=(const GammaMatrixSet&);

 public:
  virtual void init_GM() = 0;

  GammaMatrix get_GM(GMspecies spec)
  {
    assert(spec < m_Nspecies);
    return m_gm[spec];
  }

  virtual void print() = 0;

#ifdef USE_FACTORY
 public:
  typedef GammaMatrixSet *(*ProductCreator)();
  typedef FactoryTemplate<GammaMatrixSet, ProductCreator> Factory;

  static GammaMatrixSet *New(const IdentifierType& subtype)
  {
    ProductCreator p = Factory::Find(subtype);

    return p ? (*p)() : 0;
  }

#ifdef USE_FACTORY_AUTOREGISTER
#else
  static bool init_factory();
#endif
#endif
};
#endif
