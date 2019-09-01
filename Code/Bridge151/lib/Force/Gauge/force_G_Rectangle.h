/*!
        @file    force_G_Rectangle.h

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/

#ifndef FORCE_G_RECTANGLE_INCLUDED
#define FORCE_G_RECTANGLE_INCLUDED

#include "force_G.h"

#include "Measurements/Gauge/staple_lex.h"
#include "Field/shiftField_lex.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//! HMC force class BAPI for rectangular gauge action.

/*!
    Gauge action with plaquette and rectangular Wilson loops.
    Iwasaki, Luscher-Weisz, DBW2 are examples of this element_type
    of action.
                                   [28 Dec 2011 H.Matsufuru]
    (Coding history will be recovered from trac.)
    YAML is implemented.           [14 Nov 2012 Y.Namekawa]
 */

class BAPI Force_G_Rectangle : public Force_G
{
 public:
  static const std::string class_name;

 private:
  //- NB. m_U has been defined in force_G.h
  // Field_G *m_U;

  double m_beta;
  double m_c_plaq;
  double m_c_rect;

  Staple_lex m_staple;
  ShiftField_lex m_shift;

 public:
  Force_G_Rectangle()
    : Force_G() {}

  ~Force_G_Rectangle() {}

  void set_parameters(const Parameters& params);
  void set_parameters(const double beta,
                      const double c_plaq, const double c_rect);

  //- NB. set_config has been defined in force_G.h
  // void set_config(Field *U);

  void force_core(Field&);

#ifdef USE_FACTORY
 private:
  static Force_G *create_object()
  {
    return new Force_G_Rectangle();
  }

 public:
  static bool register_factory()
  {
    return Force_G::Factory::Register("Force_G_Rectangle", create_object);
  }
#endif
};
#endif
