/*!
        @file    shiftField_eo.h

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/

#ifndef SHIFTFIELD_EO_INCLUDED
#define SHIFTFIELD_EO_INCLUDED

#include <iostream>

#include "field.h"
#include "index_eo.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//! Methods to shift the even-odd field.

/*!
    This class BAPI defines the methods that shift a given Field
    instance with even-odd site index in specified direction.
    The field must have half the size of the numbers of sites
    in a node, and field on even(odd) sites are shifted to
    a field on odd (even) sites.
    Thus the method is called with an argument specifying
    one of them: ieo (0: even<-odd, 1: odd<-even).
    The forward shift means, e.g. in mu-direction,
    v(site) = w(site-\hat{mu}), where v is the shifted field
    (output, the first argument) and w the original field
    (input, the second argument).
    At present, there is no implementation for a field which
    has both the even and odd entries (and should be written).
    Names of private functions might be still confusing.
                                   [25 Dec 2011 H.Matsufuru]
 */

class BAPI ShiftField_eo {
 public:
  static const std::string class_name;

 private:
  int m_Nx, m_Ny, m_Nz, m_Nt, m_Nx2;
  Index_eo m_index_eo;

  Bridge::VerboseLevel m_vl;

 public:
  ShiftField_eo() :
    m_Nx(CommonParameters::Nx()),
    m_Ny(CommonParameters::Ny()),
    m_Nz(CommonParameters::Nz()),
    m_Nt(CommonParameters::Nt()),
    m_Nx2(CommonParameters::Nx() / 2),
    m_vl(CommonParameters::Vlevel())
  {
  }

 private:
  // non-copyable
  ShiftField_eo(const ShiftField_eo&);
  ShiftField_eo& operator=(const ShiftField_eo&);

 public:
  void forward_h(Field&, const Field&, const int mu, const int ieo);
  void backward_h(Field&, const Field&, const int mu, const int ieo);

  //                                    direction, ieo (0: e<-o, 1: o<-e)

  void forward_h(Field&, const Field&, const int boundary_condition, const int mu,
                 const int ieo);
  void backward_h(Field&, const Field&, const int boundary_condition, const int mu,
                  const int ieo);

  //                boundary condition, direction, ieo (0: e<-o, 1: o<-e)

  void forward(Field&, const Field&, const int mu);
  void backward(Field&, const Field&, const int mu);

  void forward(Field&, const Field&, const int boundary_condition, const int mu);
  void backward(Field&, const Field&, const int boundary_condition, const int mu);



 private:
  void up_xh(Field&, const Field&, const int, const int);
  void up_yh(Field&, const Field&, const int, const int);
  void up_zh(Field&, const Field&, const int, const int);
  void up_th(Field&, const Field&, const int, const int);

  void dn_xh(Field&, const Field&, const int, const int);
  void dn_yh(Field&, const Field&, const int, const int);
  void dn_zh(Field&, const Field&, const int, const int);
  void dn_th(Field&, const Field&, const int, const int);
};
#endif
