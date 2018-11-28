/*!
        @file    $Id:: force_G_Rectangle.h #$

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 18:35:38 #$

        @version $LastChangedRevision: 1571 $
*/

#ifndef FORCE_G_RECTANGLE_INCLUDED
#define FORCE_G_RECTANGLE_INCLUDED

#include "force_G.h"

#include "Measurements/Gauge/staple_lex.h"
#include "Field/shiftField_lex.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//! HMC force class for rectangular gauge action.

/*!
    Gauge action with plaquette and rectangular Wilson loops.
    Iwasaki, Luscher-Weisz, DBW2 are examples of this type
    of action.
                                   [28 Dec 2011 H.Matsufuru]
    (Coding history will be recovered from trac.)
    YAML is implemented.           [14 Nov 2012 Y.Namekawa]
 */

class Force_G_Rectangle : public Force_G
{
 public:
  static const std::string class_name;

 private:
  //- NB. m_U has been defined in force_G.h
  // Field_G *m_U;

  double m_beta;
  double m_c_plaq;
  double m_c_rect;

  Staple_lex     m_staple;
  ShiftField_lex m_shift;

 public:
  Force_G_Rectangle()
    : Force_G() {}

  ~Force_G_Rectangle() {}

  void set_parameters(const Parameters& params);
  void set_parameters(double beta, double c_plaq, double c_rect);

  //- NB. set_config has been defined in force_G.h
  // void set_config(Field *U);

  void force_core(Field&);
};
#endif
