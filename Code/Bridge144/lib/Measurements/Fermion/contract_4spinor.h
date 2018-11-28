/*!
        @file    $Id:: contract_4spinor.h #$

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 18:35:38 #$

        @version $LastChangedRevision: 1571 $
*/

#ifndef CONTRACT_4SPINOR_INCLUDED
#define CONTRACT_4SPINOR_INCLUDED

#include "Field/field_F.h"
#include "Tools/gammaMatrix.h"

#include "bridge_complex.h"

//! Contraction of hadron for 4-spinor fermion.

/*!
    This class calculates contraction of hadron for 4-spinor fermion.
                                 [15 Mar 2012 H.Matsufuru]
    (Coding history will be recovered from trac.)
    Add momentum of sink.        [21 Mar 2015 Y.Namekawa]
 */


//! contraction for meson at a given time t.
void contract_at_t(dcomplex& corr,
                   const GammaMatrix& gm_sink,
                   const Field_F& f1, const Field_F& f2,
                   const int time);

//! contraction for meson at a given time t.
void contract_at_t(dcomplex& corr,
                   const std::vector<int>& momentum_sink,
                   const GammaMatrix& gm_sink,
                   const std::vector<int>& source_position,
                   const Field_F& f1, const Field_F& f2,
                   const int time);

//! contraction for baryon (Nc=3 case only) at a given time t.
void contract_at_t(dcomplex& corr,
                   const GammaMatrix& gm_sink, const int i_alpha,
                   const Field_F& f1, const Field_F& f2, const Field_F& f3,
                   const int time);
#endif /* CONTRACT_H_INCLUDED */
