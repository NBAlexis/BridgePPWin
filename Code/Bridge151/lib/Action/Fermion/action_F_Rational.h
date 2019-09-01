/*!
        @file    action_F_Rational.h

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/

#ifndef ACTION_F_RATIONAL_INCLUDED
#define ACTION_F_RATIONAL_INCLUDED

#include "Action/action.h"

#include "Fopr/fopr_Rational.h"
#include "Force/Fermion/force_F_Rational.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//! action class BAPI for RHMC, with externally constructed Fopr_Rational.

/*!
    For the class BAPI, Fopr and Force objects are instantiated outside
    the class BAPI and specified at the construction.
    This class BAPI just provides the framework of rational actions.
                                        [28 Dec 2011 H.Matsufuru]
    unique_ptr is introduced to avoid memory leaks
                                        [21 Mar 2015 Y.Namekawa]
 */


class BAPI Action_F_Rational : public Action
{
 public:
  static const std::string class_name;

 private:
  std::string m_label;    // label of action

  Fopr *m_fopr_langev;
  Fopr *m_fopr_H;
  Force *m_fopr_force_MD;

  Field *m_U;

  Field m_psf;

 public:
  //! constructor requires pointers to Fopr and Force instances.
  Action_F_Rational(
    Fopr *fopr_langev, Fopr *fopr_H, Force *fopr_force_MD)
    : Action(),
    m_fopr_langev(fopr_langev), m_fopr_H(fopr_H), m_fopr_force_MD(fopr_force_MD)
  {}

  Action_F_Rational(
    unique_ptr<Fopr>& fopr_langev, unique_ptr<Fopr>& fopr_H, unique_ptr<Force>& fopr_force_MD)
    : Action(),
    m_fopr_langev(fopr_langev.get()), m_fopr_H(fopr_H.get()), m_fopr_force_MD(fopr_force_MD.get())
  {}

  //! destructor. constructed instances are deconstructed in tydyup().
  ~Action_F_Rational() {}

  void set_parameters(const Parameters& params);

  //! set the label of action.
  void set_label(const std::string label)
  {
    m_label = label;
    vout.detailed(m_vl, "  label: %s\n", m_label.c_str());
  }

  //! returns the label of action.
  std::string get_label()
  {
    return m_label;
  }

  //! setting gauge configuration.
  void set_config(Field *U)
  {
    m_U = U;
    m_fopr_langev->set_config(U);
    m_fopr_H->set_config(U);
    m_fopr_force_MD->set_config(U);
  }

  //! Langevin step called at the beginning of HMC.
  double langevin(RandomNumbers *);

  //! calculation of Hamiltonian.
  double calcH();

  //! returns the force for updating conjugate momentum.
  //const Field force();
  void force(Field&);

 private:
};
#endif
