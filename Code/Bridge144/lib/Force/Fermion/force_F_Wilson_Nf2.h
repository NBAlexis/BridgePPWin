/*!
        @file    $Id:: force_F_Wilson_Nf2.h #$

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 18:35:38 #$

        @version $LastChangedRevision: 1571 $
*/

#ifndef FORCE_F_WILSON_NF2_INCLUDED
#define FORCE_F_WILSON_NF2_INCLUDED

#include "force_F.h"
#include "Fopr/fopr_Wilson.h"

#include "tensorProd.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//! Force for the standard Wilson fermion operator

/*!
    This class calculates the force of the standard Wilson fermion.
    The gamma matrix representation is given as control string
    "Dirac"(default) or "Chiral" at the construction, which is
    used to construct the Fopr_Wilson instance.
                                     [23 Dec 2011 H.Matusfuru]
    (Coding history will be recovered from trac.)
    YAML is implemented.             [14 Nov 2012 Y.Namekawa]
 */


class Force_F_Wilson_Nf2 : public Force
{
 public:
  static const std::string class_name;

 private:
  double           m_kappa;
  std::vector<int> m_boundary;
  Fopr_Wilson      *m_fopr_w;
  Field_F          m_psf;
  std::string      m_repr;

 public:
  Force_F_Wilson_Nf2()
    : Force()
  {
    m_repr   = "Dirac";
    m_fopr_w = new Fopr_Wilson(m_repr);
    m_boundary.resize(CommonParameters::Ndim());
  }

  Force_F_Wilson_Nf2(std::string repr)
    : Force()
  {
    m_repr   = repr;
    m_fopr_w = new Fopr_Wilson(m_repr);
    m_boundary.resize(CommonParameters::Ndim());
  }

  ~Force_F_Wilson_Nf2()
  {
    delete m_fopr_w;
  }

  void set_parameters(const Parameters& params);

  void set_parameters(const double kappa, const std::vector<int> bc);

  void set_config(Field *U)
  {
    m_U = (Field_G *)U;
    m_fopr_w->set_config(U);
  }

  void force_udiv(Field& force, const Field& eta);
  void force_udiv1(Field& force, const Field& zeta, const Field& eta);

 private:
  void force_udiv1_impl(Field_G& force, const Field_F& zeta, const Field_F& eta);
};
#endif
