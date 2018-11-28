/*!
        @file    $Id:: force_F_Wilson_Nf2_Isochemical.h #$

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 18:35:38 #$

        @version $LastChangedRevision: 1571 $
*/

#ifndef FORCE_F_WILSON_NF2_ISOCHEMICAL_INCLUDED
#define FORCE_F_WILSON_NF2_ISOCHEMICAL_INCLUDED

#include "force_F.h"
#include "Fopr/fopr_Wilson_Isochemical.h"

#include "tensorProd.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//! Force for the Wilson fermion operator with isospin chemical potential.

/*!
    This class calculates the force of the standard Wilson fermion with
    isospin chemical potential with two flavors.
                                     [24 Aug 2011 H.Matusfuru]
    (Coding history will be recovered from trac.)
    YAML is implemented.             [14 Nov 2012 Y.Namekawa]
 */


class Force_F_Wilson_Nf2_Isochemical : public Force
{
 public:
  static const std::string class_name;

 private:
  double                  m_kappa;     //!< hopping parameter
  double                  m_mu;        //!< isospin chemical potential
  double                  m_exp_mu;    //!< exp(mu)
  std::vector<int>        m_boundary;
  Fopr_Wilson_Isochemical *m_fopr_w;
  Field_F                 m_psf;

  std::string m_repr;
  std::string m_mode;

 public:

  Force_F_Wilson_Nf2_Isochemical()
    : Force()
  {
    m_repr   = "Dirac";
    m_fopr_w = new Fopr_Wilson_Isochemical(m_repr);
    m_boundary.resize(CommonParameters::Ndim());
  }

  Force_F_Wilson_Nf2_Isochemical(std::string repr)
    : Force()
  {
    m_repr   = repr;
    m_fopr_w = new Fopr_Wilson_Isochemical(m_repr);
    m_boundary.resize(CommonParameters::Ndim());
  }

  ~Force_F_Wilson_Nf2_Isochemical()
  {
    delete m_fopr_w;
  }

  void set_parameters(const Parameters& params);

  void set_parameters(const double kappa, const double mu, const std::vector<int> bc);

  void set_config(Field *U)
  {
    m_U = (Field_G *)U;
    m_fopr_w->set_config(U);
  }

  void set_mode(const std::string& mode)
  {
    m_mode = mode;
  }

  void force_udiv(Field& force, const Field& eta);
  void force_udiv1(Field& force, const Field& zeta, const Field& eta);

 private:
  void force_udiv1_impl(Field_G& force, const Field_F& zeta, const Field_F& eta);
};
#endif
