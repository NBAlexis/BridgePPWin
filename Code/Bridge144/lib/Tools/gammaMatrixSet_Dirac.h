/*!
        @file    $Id:: gammaMatrixSet_Dirac.h #$

        @brief

        @author  Hideo Matsufuru  (matsufuru)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 00:49:46 #$

        @version $LastChangedRevision: 1561 $
*/

#ifndef GAMMAMATRIXSET_DIRAC_INCLUDED
#define GAMMAMATRIXSET_DIRAC_INCLUDED

#include "gammaMatrixSet.h"


//! Set of Gamma Matrix: Dirac representation.

/*!
                                        [4 Feb 2012 H.Matsufuru]
 */
class GammaMatrixSet_Dirac : public GammaMatrixSet {
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
};
#endif
