/*!
        @file    $Id:: action_F_Rational_SF.h #$

        @brief

        @author  Yusuke Taniguchi
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 18:35:38 #$

        @version $LastChangedRevision: 1571 $
*/

#ifndef ACTION_FRATIONAL_SF_INCLUDED
#define ACTION_FRATIONAL_SF_INCLUDED

#include "Action/action.h"

#include "Fopr/fopr_Rational.h"
#include "Force/Fermion/force_F_Rational.h"

#include "Field/field_G_SF.h"
#include "Field/field_F_SF.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//! action class for RHMC, with externally constructed Fopr_Rational.

/*!
    For the class, Fopr and Force objects are instantiated outside
    the class and specified at the construction.
    This class just provides the framework of rational actions.
                                         [28 Dec 2011 H.Matsufuru]
    unique_ptr is introduced to avoid memory leaks
                                         [21 Mar 2015 Y.Namekawa]
 */


class Action_F_Rational_SF : public Action {
 public:
  static const std::string class_name;

 private:
  std::string m_label;   // label of action

  Fopr  *m_fopr_langev;
  Fopr  *m_fopr_H;
  Force *m_fopr_force_MD;

  Field *m_U;

  Field m_psf;

 public:
  //! constructor requires pointers to Fopr and Force instances.
  Action_F_Rational_SF(Fopr *fopr_langev, Fopr *fopr_H,
                       Force *fopr_force_MD)
  {
    m_fopr_langev   = fopr_langev;
    m_fopr_H        = fopr_H;
    m_fopr_force_MD = fopr_force_MD;
    setup();
  }

  Action_F_Rational_SF(unique_ptr<Fopr>& fopr_langev, unique_ptr<Fopr>& fopr_H,
                       unique_ptr<Force>& fopr_force_MD)
  {
    m_fopr_langev   = fopr_langev.get();
    m_fopr_H        = fopr_H.get();
    m_fopr_force_MD = fopr_force_MD.get();
    setup();
  }

  //! destructor. constructed instances are deconstructed in tydyup().
  ~Action_F_Rational_SF()
  {
  }

  //! setting parameters and creating class instances.
  void set_parameters(const Parameters& params);

  //! set the label of action.
  void set_label(std::string label)
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
  //! creating instances. called from set_parameters().
  void setup();

  //! destruct class instances constructed in setup()
  void tidyup();
};
#endif
