/*!
        @file    $Id:: force_F_Clover_Nf2_Isochemical.h #$

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 18:35:38 #$

        @version $LastChangedRevision: 1571 $
*/

#ifndef FORCE_F_CLOVER_NF2_ISOCHEMICAL_INCLUDED
#define FORCE_F_CLOVER_NF2_ISOCHEMICAL_INCLUDED

#include "Fopr/fopr_Clover_Isochemical.h"

#include "force_F_Wilson_Nf2_Isochemical.h"
#include "force_F_CloverTerm.h"

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


class BAPI Force_F_Clover_Nf2_Isochemical : public Force
{
 public:
  static const std::string class_name;

 private:
  double           m_kappa;      //!< hopping parameter
  double           m_cSW;        //!< clover coefficient
  double           m_mu;         //!< Isospin chemical potential
  std::vector<int> m_boundary;   //!< boundary conditions

  std::string m_repr;            //!< gamma matrix representation
  std::string m_mode;            //!< mult mode

//  Field_G *m_U;                //!< pointer to gauge field
//  Field_G* m_Cud;              //!< for force calculation

  Fopr_Clover_Isochemical        *m_fopr_c;    //!< fermion operator
  Force_F_Wilson_Nf2_Isochemical *m_force_w;   //!< Wilson fermion force
  Force_F_CloverTerm             *m_force_csw; //!< Clover term force

  int m_Ndim;                                  //!< spacetime dimension

 public:
  //! Constructor
  Force_F_Clover_Nf2_Isochemical()
    : Force()
  {
    init("Dirac");  //!< default gamma matrix representation
  }

  Force_F_Clover_Nf2_Isochemical(std::string repr)
    : Force()
  {
    init(repr);
  }

  //! Deconstructor
  ~Force_F_Clover_Nf2_Isochemical()
  {
    tidyup();
  }

  //! Setting parameters of clover fermion force.
  void set_parameters(const Parameters& params);

  //! Setting parameters with clover fermion parameters.
  void set_parameters(const double kappa, const double cSW, const double mu,
                      const std::vector<int> bc);

  void set_mode(const std::string& mode)
  {
    m_mode = mode;
    m_fopr_c->set_mode(mode);
    m_force_w->set_mode(mode);
  }

  //! Setting gauge configuration
  void set_config(Field *U)
  {
    m_U = (Field_G *)U;
    m_fopr_c->set_config(U);
    m_force_w->set_config(U);
    m_force_csw->set_config(U);
    //    set_component();
  }

  //! For recursive calculation of smeared force.
  void force_udiv(Field& force, const Field& eta);

  //! For recursive calculation of smeared force.
  void force_udiv1(Field& force, const Field& zeta, const Field& eta);


 private:
  void init(std::string repr);
  void tidyup();

  //! Core implemetation of clover force calculation.
  void force_udiv1_impl(Field_G& force, const Field_F& zeta, const Field_F& eta);

  //! Set building components for force calculation.
  //  void set_component();

  int index_dir(int mu, int nu)
  {
    return mu + m_Ndim * nu;
  }
};
#endif
