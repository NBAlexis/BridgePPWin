/*!
        @file    $Id:: action_F_Standard_eo.h #$

        @brief

        @author  UEDA, Satoru  (sueda)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 18:35:38 #$

        @version $LastChangedRevision: 1571 $
*/

#ifndef ACTION_F_STANDARD_EO_INCLUDED
#define ACTION_F_STANDARD_EO_INCLUDED

#include "Action/action.h"

#include "Force/Fermion/force_F.h"
#include "Measurements/Fermion/fprop.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//! Standard even-odd preconditioned fermion action for HMC.

/*!
    This class is used to define an action used in HMC.
    Fermion and Force operators and given at the construction.
                                        [19 Jun 2012 S.UEDA]
    (Coding history will be recovered from trac.)
    Modify this code to work.           [03 Mar 2013 Y.Namekawa]
    unique_ptr is introduced to avoid memory leaks
                                        [21 Mar 2015 Y.Namekawa]
 */

class Action_F_Standard_eo : public Action
{
 public:
  static const std::string class_name;

 private:
  Fopr        *m_fopr;
  Force       *m_fopr_force;
  Field       m_psf;
  std::string m_label;

  Fprop *m_fprop_MD;
  Fprop *m_fprop_H;

  Field *m_U;


 public:
  Action_F_Standard_eo(
    Fopr *fopr, Force *fopr_force,
    Fprop *fprop_MD, Fprop *fprop_H)
    : Action(),
      m_fopr(fopr), m_fopr_force(fopr_force),
      m_fprop_MD(fprop_MD), m_fprop_H(fprop_H)
  {
    set_parameters();
  }

  Action_F_Standard_eo(
    unique_ptr<Fopr>& fopr, unique_ptr<Force>& fopr_force,
    unique_ptr<Fprop>& fprop_MD, unique_ptr<Fprop>& fprop_H)
    : Action(),
      m_fopr(fopr.get()), m_fopr_force(fopr_force.get()),
      m_fprop_MD(fprop_MD.get()), m_fprop_H(fprop_H.get())
  {
    set_parameters();
  }

  ~Action_F_Standard_eo() {}

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
