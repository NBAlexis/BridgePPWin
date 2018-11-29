/*!
        @file    $Id:: field_F_SF.h #$

        @brief

        @author  Yusuke Taniguchi
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 00:49:46 #$

        @version $LastChangedRevision: 1561 $
*/


#ifndef FIELD_F_SF_INCLUDED
#define FIELD_F_SF_INCLUDED

#include "Parameters/commonParameters.h"
#include "Communicator/communicator.h"
//#include "field_F.h"
//#include "field_G.h"
#include "Tools/vec_SU_N.h"

//! A class BAPI generated to add a function for the SF.

/*!
  <ul>
  <li>A function to set the boundary field.
  <li>This class BAPI does not contain the filed object but manipulate it.
  <li>[3 Apr 2012 Y.Taniguchi]
  </ul>
 */
class BAPI Field_F_SF {
 private:
  int m_Nc2;  // num of the double color elements
  int m_Nvol; // lattice volume
  int m_Nd;   // num of the spinor elements
  int m_Nin;  // internal d.o.f.
  int Svol;

 public:
  Field_F_SF() :
    m_Nvol(CommonParameters::Nvol()),
    m_Nd(CommonParameters::Nd())
  {
    int Nc = CommonParameters::Nc();

    m_Nc2 = 2 * Nc;
    int Nt = CommonParameters::Nt();
    Svol = m_Nvol / Nt;
  }

/*!
  Set the boundary field to zero: \f$f(t=0,\vec{x})=0\f$
 */
  void set_boundary_zero(Field& f)
  {
    if (Communicator::ipe(3) == 0) {
      for (int site = 0; site < Svol; ++site) {
        for (int s = 0; s < m_Nd; ++s) {
          for (int cc = 0; cc < m_Nc2; ++cc) {
            f.set(cc + m_Nc2 * s, site, 0, 0.0);
          }
        }
      }
    }
  }
};
#endif
