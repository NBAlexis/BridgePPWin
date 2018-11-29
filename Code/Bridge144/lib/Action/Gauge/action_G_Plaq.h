/*!
        @file    $Id:: action_G_Plaq.h #$

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 18:35:38 #$

        @version $LastChangedRevision: 1571 $
*/

#ifndef ACTION_G_PLAQ_INCLUDED
#define ACTION_G_PLAQ_INCLUDED

#include "Action/action.h"
#include "Force/Gauge/force_G_Plaq.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//! HMC action class BAPI for plaquette gauge action.

/*!
    Standard plaquette gauge action.
                             [28 Dec 2011 H.Matsufuru]
    (Coding history will be recovered from trac.)
    YAML is implemented.     [14 Nov 2012 Y.Namekawa]
 */


class BAPI Action_G_Plaq : public Action
{
 public:
  static const std::string class_name;

 private:
  double      m_beta;
  std::string m_label;

  Field_G    *m_U;
  Staple_lex m_staple;
  Force_G    *m_force_G;

 public:
  Action_G_Plaq()
    : Action()
  {
    m_force_G = Force_G::New("Force_G_Plaq");
  }

  ~Action_G_Plaq()
  {
    delete m_force_G;
  }

  void set_parameters(const Parameters& params);
  void set_parameters(double beta);

  void set_label(std::string label)
  {
    m_label = label;
    vout.detailed(m_vl, "  label: %s\n", m_label.c_str());
  }

  std::string get_label()
  {
    return m_label;
  }

  void set_config(Field *U)
  {
    m_U = (Field_G *)U;
  }

  double langevin(RandomNumbers *);

  double calcH();

  void force(Field&);
};
#endif
