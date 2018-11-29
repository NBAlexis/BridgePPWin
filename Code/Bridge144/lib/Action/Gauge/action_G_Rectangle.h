/*!
        @file    $Id:: action_G_Rectangle.h #$

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 18:35:38 #$

        @version $LastChangedRevision: 1571 $
*/

#ifndef ACTION_G_RECTANGLE_INCLUDED
#define ACTION_G_RECTANGLE_INCLUDED

#include "Action/action.h"
#include "Force/Gauge/force_G_Rectangle.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//! HMC action class BAPI for rectangular gauge action.

/*!
    Gauge action with plaquette and rectangular Wilson loops.
    Iwasaki, Luscher-Weisz, DBW2 are examples of this type
    of action.
                                   [28 Dec 2011 H.Matsufuru]
    (Coding history will be recovered from trac.)
    YAML is implemented.           [14 Nov 2012 Y.Namekawa]
 */


class BAPI Action_G_Rectangle : public Action
{
 public:
  static const std::string class_name;

 private:
  double      m_beta;
  double      m_c_plaq;
  double      m_c_rect;
  std::string m_label;

  Field_G        *m_U;
  Staple_lex     m_staple;
  ShiftField_lex m_shift;
  Force_G        *m_force_G;

 public:
  Action_G_Rectangle()
    : Action()
  {
    m_force_G = Force_G::New("Force_G_Rectangle");
  }

  ~Action_G_Rectangle()
  {
    delete m_force_G;
  }

  void set_parameters(const Parameters& params);
  void set_parameters(double beta, double c_plaq, double c_rect);

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
