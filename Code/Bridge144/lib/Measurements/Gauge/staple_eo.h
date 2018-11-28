/*!
        @file    $Id:: staple_eo.h #$

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 18:35:38 #$

        @version $LastChangedRevision: 1571 $
*/

#ifndef STAPLE_EO_INCLUDED
#define STAPLE_EO_INCLUDED

#include "staple.h"

#include "Field/shiftField_eo.h"

#include "IO/bridgeIO.h"

//! Staple construction.

/*!
    While the originial version was written by J.Noaki,
    the present version is completely modified by H.Matsufuru
    except for the interface.
                                    [28 Dec 2011 H.Matsufuru]
    Factory is introduced.          [24 Jan 2017 Y.Namekawa]
 */

class Staple_eo : public Staple
{
 public:
  static const std::string class_name;

 private:
  std::string m_filename_output;

  Field_G       m_v, m_w;
  ShiftField_eo m_shift;

 public:
  Staple_eo() : Staple()
  {
    m_filename_output = "stdout";
  }

  ~Staple_eo() {}

 public:
  void set_parameters(const Parameters& params);

  void upper(Field_G&, const Field_G&, const int, const int);
  void lower(Field_G&, const Field_G&, const int, const int);
  double plaq_s(const Field_G&);
  double plaq_t(const Field_G&);
  double plaquette(const Field_G&);

  void staple(Field_G&, const Field_G&, const int);
};
#endif
