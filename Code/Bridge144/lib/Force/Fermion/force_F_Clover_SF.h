/*!
        @file    $Id:: force_F_Clover_SF.h #$

        @brief

        @author  Yusuke Taniguchi
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 18:35:38 #$

        @version $LastChangedRevision: 1571 $
*/

#ifndef FORCE_F_CLOVER_SF_INCLUDED
#define FORCE_F_CLOVER_SF_INCLUDED

#include "Fopr/fopr_Clover_SF.h"
#include "force_F_Wilson_SF.h"
#include "Field/field_G_SF.h"
#include "Measurements/Gauge/staple_SF.h"

#include "tensorProd.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//! Force calculation for clover quark action with SF BC.

/*!
    At present, only the Dirac representation for gamma-matrix
    is available.
                                [28 Dec 2011 H.Matsufuru]
    (Coding history will be recovered from trac.)
    YAML is implemented.        [14 Nov 2012 Y.Namekawa]
 */


class BAPI Force_F_Clover_SF : public Force
{
 public:
  static const std::string class_name;

 private:
  double           m_kappa;      //!< hopping parameter
  double           m_cSW;        //!< clover coefficient
  std::vector<int> m_boundary;   //!< boundary conditions

  int     m_Ndim;
  Field_G *m_Cud;                  //!< for force calculation

  Fopr_Clover_SF    *m_fopr_c;
  Force_F_Wilson_SF *m_force_w;
  Force_F_Clover_SF *m_imp;

  //! SF boundary condition at t=0
  double m_phi[3];
  //! SF boundary condition at t=Nt
  double m_phipr[3];

  //! In order to set the boundary field.
  Field_F_SF set_zero;
  Field_G_SF set_wk;

 public:

  Force_F_Clover_SF()
  {
    m_fopr_c  = new Fopr_Clover_SF;
    m_force_w = new Force_F_Wilson_SF;

    int Nvol = CommonParameters::Nvol();
    m_Ndim = CommonParameters::Ndim();

    m_boundary.resize(m_Ndim);
    m_Cud = new Field_G(Nvol, m_Ndim * m_Ndim);
  }

  ~Force_F_Clover_SF()
  {
    delete m_Cud;
    delete m_force_w;
    delete m_fopr_c;
  }

  void set_parameters(const Parameters& params);

  //! Setting parameters of clover fermion.
  void set_parameters(double kappa, double cSW, const std::vector<int> bc,
                      double *phi, double *phipr);

  //! Setting gauge configuration
  void set_config(Field *U)
  {
    m_U = (Field_G *)U;
    m_fopr_c->set_config(m_U);
    m_force_w->set_config(U);
    set_component();
  }

  //! For recursive calculation of smeared force.
  void force_udiv(Field& force, const Field& eta);

  //! For recursive calculation of smeared force.
  void force_udiv1(Field& force, const Field& zeta, const Field& eta);


 private:

  //! Core implemetation of clover force calculation.
  void force_udiv1_impl(Field_G& force, const Field_F& zeta, const Field_F& eta);

  //! Set building components for force calculation.
  void set_component();

  int index_dir(int mu, int nu)
  {
    return mu + m_Ndim * nu;
  }
};
#endif
