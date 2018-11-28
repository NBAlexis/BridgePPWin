/*!
        @file    $Id:: gammaMatrixSet_Chiral.h #$

        @brief

        @author  Hideo Matsufuru  (matsufuru)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 00:49:46 #$

        @version $LastChangedRevision: 1561 $
*/

#ifndef GAMMAMATRIXSET_CHIRAL_INCLUDED
#define GAMMAMATRIXSET_CHIRAL_INCLUDED

#include "gammaMatrixSet.h"


//! Set of Gamma Matrix: chiral representation.

/*!
                                        [4 Feb 2012 H.Matsufuru]
 */

class GammaMatrixSet_Chiral : public GammaMatrixSet
{
 public:
  static const std::string class_name;

 public:
  GammaMatrixSet_Chiral()
  {
    init_GM();
  }

  void print();

 private:
  void init_GM();
};
#endif
