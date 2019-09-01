/*!
        @file    force_F_CloverTerm.h

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/

#ifndef FORCE_F_CLOVERTERM_INCLUDED
#define FORCE_F_CLOVERTERM_INCLUDED

#include "force_F.h"
#include "tensorProd.h"

#include "Field/shiftField_lex.h"
#include "Fopr/fopr_CloverTerm.h"
#include "Measurements/Gauge/staple_lex.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//! Force calculation for clover term of clover fermion.

/*!
    This class BAPI calculates contribution clover term to the clover
    fermion action.
                                [28 Dec 2011 H.Matsufuru]
    (Coding history will be recovered from trac.)
    YAML is implemented.        [14 Nov 2012 Y.Namekawa]
 */


class BAPI Force_F_CloverTerm : public Force
{
 public:
  static const std::string class_name;

 private:
  double m_kappa;              //!< hopping parameter
  double m_cSW;                //!< clover coefficient
  std::vector<int> m_boundary; //!< boundary conditions

  std::string m_repr;          //!< gamma matrix representation

//  Field_G         *m_U;          //!< pointer to gauge field
  int m_Ndim;                    //!< spacetime dimension
  Field_G *m_Cud;                //!< for force calculation
  Fopr_CloverTerm *m_fopr_csw;   //!< fermion operator

 public:
  //! Constructor
  Force_F_CloverTerm()
    : Force()
  {
    init("Dirac");    //!< default gamma matrix representation
  }

  Force_F_CloverTerm(const std::string repr)
    : Force()
  {
    init(repr);    //!< default gamma matrix representation
  }

  //! Deconstructor
  ~Force_F_CloverTerm()
  {
    tidyup();
  }

  //! Setting parameters of clover fermion force.
  void set_parameters(const Parameters& params);

  //! Setting parameters with clover fermion parameters.
  void set_parameters(const double kappa, const double cSW, const std::vector<int> bc);

  //! Setting gauge configuration
  void set_config(Field *U)
  {
    m_U = (Field_G *)U;
    m_fopr_csw->set_config(m_U);
    //    m_forcew->set_config(U);
    set_component();
  }

  //! Force determination for clover fermion.
  void force_core(Field& force, const Field& eta);  // override default

//   //! Force determination for clover fermion.
//   void force_core1(Field& force, const Field& zeta, const Field& eta);

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
