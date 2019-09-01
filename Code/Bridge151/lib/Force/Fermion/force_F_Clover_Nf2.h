/*!
        @file    force_F_Clover_Nf2.h

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/

#ifndef FORCE_F_CLOVER_NF2_INCLUDED
#define FORCE_F_CLOVER_NF2_INCLUDED

#include "Fopr/fopr_Clover.h"

#include "force_F_Wilson_Nf2.h"
#include "force_F_CloverTerm.h"

#include "Measurements/Gauge/staple_lex.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//! Force calculation for clover quark action.

/*!
    At present, only the Dirac representation for gamma-matrix
    is available.
                                [28 Dec 2011 H.Matsufuru]
    (Coding history will be recovered from trac.)
    YAML is implemented.        [14 Nov 2012 Y.Namekawa]
 */


class BAPI Force_F_Clover_Nf2 : public Force
{
 public:
  static const std::string class_name;

 private:
  double m_kappa;                //!< hopping parameter
  double m_cSW;                  //!< clover coefficient
  std::vector<int> m_boundary;   //!< boundary conditions

  std::string m_repr;            //!< gamma matrix representation

//  Field_G            *m_U;         //!< pointer to gauge field
  int m_Ndim;                      //!< spacetime dimension
  Field_G *m_Cud;                  //!< for force calculation
  Fopr_Clover *m_fopr_c;           //!< fermion operator
  Force_F_Wilson_Nf2 *m_force_w;   //!< Wilson fermion force
  Force_F_CloverTerm *m_force_csw; //!< Clover term force

 public:

  Force_F_Clover_Nf2()
    : Force()
  {
    init("Dirac");    //!< default gamma matrix representation
  }

  //! Construction with gamma matrix representation
  Force_F_Clover_Nf2(const std::string repr)
    : Force()
  {
    init(repr);
  }

  //! Deconstructor
  ~Force_F_Clover_Nf2()
  {
    tidyup();
  }

  //! Setting parameters of clover fermion force.
  void set_parameters(const Parameters& params);

  //! Setting parameters with clover fermion parameters.
  void set_parameters(const double kappa, const double cSW,
                      const std::vector<int> bc);

  //! Setting gauge configuration
  void set_config(Field *U)
  {
    m_U = (Field_G *)U;
    m_fopr_c->set_config(U);
    m_force_w->set_config(U);
    m_force_csw->set_config(U);
    set_component();
  }

  //! For recursive calculation of smeared force.
  void force_udiv(Field& force, const Field& eta);

  //! For recursive calculation of smeared force.
  void force_udiv1(Field& force, const Field& zeta, const Field& eta);

 private:
  void init(const std::string);
  void tidyup();

  //! Core implemetation of clover force calculation.
  void force_udiv1_impl(Field_G& force, const Field_F& zeta, const Field_F& eta);

  //! Set building components for force calculation.
  void set_component();

  int index_dir(const int mu, const int nu)
  {
    return mu + m_Ndim * nu;
  }
};
#endif
