/*!
        @file    $Id:: force_F_Wilson_eo.h #$

        @brief

        @author  UEDA, Satoru (sueda)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 18:35:38 #$

        @version $LastChangedRevision: 1571 $
*/

#ifndef FORCE_F_WILSON_EO_INCLUDED
#define FORCE_F_WILSON_EO_INCLUDED

#include "force_F.h"
#include "Fopr/fopr_Wilson_eo.h"

#include "Tools/gammaMatrix.h"

#include "tensorProd.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//! Force for the Wilson fermion operator with even-odd precondition

/*!
    This class calculates the force of the standard Wilson fermion.
    The gamma matrix representation is given as control string
    "Dirac"(default) or "Chiral" at the construction, which is
    used to construct the Fopr_Wilson instance.
                                     [19 June 2012 S.UEDA]
    (Coding history will be recovered from trac.)
    YAML is implemented.             [14 Nov 2012 Y.Namekawa]
 */


class BAPI Force_F_Wilson_eo : public Force
{
 public:
  static const std::string class_name;

 private:
  Field_G *m_Ueo;

  double           m_kappa;
  std::vector<int> m_boundary;
  Fopr_Wilson_eo   *m_fopr_w;
  Field_F          m_psf;
  std::string      m_repr;

  Index_eo m_index;

 public:

  Force_F_Wilson_eo()
    : Force()
  {
    m_repr   = "Dirac";
    m_fopr_w = new Fopr_Wilson_eo(m_repr);
    m_boundary.resize(CommonParameters::Ndim());
  }

  Force_F_Wilson_eo(std::string repr)
    : Force()
  {
    m_repr   = repr;
    m_fopr_w = new Fopr_Wilson_eo(m_repr);
    m_boundary.resize(CommonParameters::Ndim());
  }

  ~Force_F_Wilson_eo()
  {
    delete m_fopr_w;
  }

  void set_parameters(const Parameters& params);
  void set_parameters(double kappa, const std::vector<int> bc);

  void set_config(Field *U)
  {
    m_U = (Field_G *)U;
    m_fopr_w->set_config(U);
  }

  //! return the force field, differentiated by the gauge field U form the fermion field on the even site.
  //! force_udiv = force_udiv1(eta, H eta) + force_udiv1(H eta, eta)
  void force_udiv(Field& force, const Field& eta);

  //! eta and zeta are the fermion fields on the even-odd site
  //! eta_o = M_oe eta_e
  //! zeta_e = g_5 (1- M_eo M_oe) eta_e
  //! zeta_o = M_oe zeta_o
  void force_udiv1(Field& force, const Field& zeta, const Field& eta);

 private:
  void force_udiv1_impl(Field_G& force, const Field_F& zeta, const Field_F& eta);
};
#endif
