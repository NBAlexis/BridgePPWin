/*!
        @file    gammaMatrixSet_Dirac.h

        @brief

        @author  Hideo Matsufuru  (matsufuru)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/

#ifndef GAMMAMATRIXSET_DIRAC_INCLUDED
#define GAMMAMATRIXSET_DIRAC_INCLUDED

#include "gammaMatrixSet.h"


//! Set of Gamma Matrix: Dirac representation.

/*!
                                        [4 Feb 2012 H.Matsufuru]
 */
class BAPI GammaMatrixSet_Dirac : public GammaMatrixSet {
 public:
  static const std::string class_name;

 public:
  GammaMatrixSet_Dirac()
  {
    init_GM();
  }

  void print();

 private:
  void init_GM();

#ifdef USE_FACTORY
 private:
  static GammaMatrixSet *create_object()
  {
    return new GammaMatrixSet_Dirac();
  }

 public:
  static bool register_factory()
  {
    return GammaMatrixSet::Factory::Register("Dirac", create_object);
  }
#endif
};
#endif
