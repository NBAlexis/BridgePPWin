/*!
        @file    $Id:: force_G_Plaq.h #$

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 18:35:38 #$

        @version $LastChangedRevision: 1571 $
*/

#ifndef FORCE_G_PLAQ_INCLUDED
#define FORCE_G_PLAQ_INCLUDED

#include "force_G.h"

#include "Measurements/Gauge/staple_lex.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//! HMC force class for plaquette gauge action.

/*!
    Standard plaquette gauge action.
                             [28 Dec 2011 H.Matsufuru]
    (Coding history will be recovered from trac.)
    YAML is implemented.     [14 Nov 2012 Y.Namekawa]
 */

class Force_G_Plaq : public Force_G
{
 public:
  static const std::string class_name;

 private:
  //- NB. m_U has been defined in force_G.h
  // Field_G *m_U;

  double     m_beta;
  Staple_lex m_staple;

 public:
  Force_G_Plaq() : Force_G() {}

  ~Force_G_Plaq() {}

  void set_parameters(const Parameters& params);
  void set_parameters(const double beta);

  //- NB. set_config has been defined in force_G.h
  // void set_config(Field *U);

  void force_core(Field& force);
};
#endif
