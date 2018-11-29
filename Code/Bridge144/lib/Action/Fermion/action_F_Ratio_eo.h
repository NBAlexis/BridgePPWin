/*!
        @file    $Id:: action_F_Ratio_eo.h #$

        @brief

        @author  Yusuke Namekawa (namekawa)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 18:35:38 #$

        @version $LastChangedRevision: 1571 $
*/

#ifndef ACTION_F_RATIO_EO_INCLUDED
#define ACTION_F_RATIO_EO_INCLUDED

#include "Action/action.h"

#include "Force/Fermion/force_F.h"
#include "Measurements/Fermion/fprop.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//! HMC action for Hasenbusch preconditioned fermions.

/*!
    This class BAPI is used to define an fermion action used in HMC
    which is given as a ratio of two fermion operators.
    Two sets of fermion and Force operators and given at the
    construction.
                                        [05 Dec 2011 H.Matsufuru]
    (Coding history will be recovered from trac.)
    Even-odd is implemented.            [03 Mar 2013 Y.Namekawa]
    unique_ptr is introduced to avoid memory leaks
                                        [21 Mar 2015 Y.Namekawa]
 */

class BAPI Action_F_Ratio_eo : public Action
{
 public:
  static const std::string class_name;

 private:
  Field *m_U;

  Fopr        *m_fopr_prec;       // preconditioner
  Force       *m_fopr_prec_force; // force of preconditioner
  Fopr        *m_fopr;            // dynamical fermion
  Force       *m_fopr_force;      // force of dynamical fermion
  Field       m_psf;              // pseudofermion field
  std::string m_label;            // label of action

  Fprop *m_fprop_H_prec;
  Fprop *m_fprop_MD;
  Fprop *m_fprop_H;

 public:
  Action_F_Ratio_eo(
    Fopr *fopr_prec, Force *fopr_prec_force,
    Fopr *fopr, Force *fopr_force,
    Fprop *fprop_H_prec,
    Fprop *fprop_MD, Fprop *fprop_H)
    : Action(),
      m_fopr_prec(fopr_prec), m_fopr_prec_force(fopr_prec_force),
      m_fopr(fopr), m_fopr_force(fopr_force),
      m_fprop_H_prec(fprop_H_prec),
      m_fprop_MD(fprop_MD), m_fprop_H(fprop_H)
  {
    set_parameters();
  }

  Action_F_Ratio_eo(
    unique_ptr<Fopr>& fopr_prec, unique_ptr<Force>& fopr_prec_force,
    unique_ptr<Fopr>& fopr, unique_ptr<Force>& fopr_force,
    unique_ptr<Fprop>& fprop_H_prec,
    unique_ptr<Fprop>& fprop_MD, unique_ptr<Fprop>& fprop_H)
    : Action(),
      m_fopr_prec(fopr_prec.get()), m_fopr_prec_force(fopr_prec_force.get()),
      m_fopr(fopr.get()), m_fopr_force(fopr_force.get()),
      m_fprop_H_prec(fprop_H_prec.get()),
      m_fprop_MD(fprop_MD.get()), m_fprop_H(fprop_H.get())
  {
    set_parameters();
  }

  ~Action_F_Ratio_eo() {}

  void set_parameters(const Parameters&);
  void set_parameters();

  void set_label(std::string label)
  {
    m_label = label;
    vout.detailed(m_vl, "  label: %s\n", m_label.c_str());
  }

  std::string get_label()
  {
    return m_label;
  }

  void set_config(Field *U);

  double langevin(RandomNumbers *);
  double calcH();

  void force(Field&);
};
#endif
