/*!
        @file    $Id:: staple_lex.h #$

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 18:35:38 #$

        @version $LastChangedRevision: 1571 $
*/

#ifndef STAPLE_LEX_INCLUDED
#define STAPLE_LEX_INCLUDED

#include "staple.h"

#include "Field/shiftField_lex.h"

#include "IO/bridgeIO.h"

//! Staple construction.

/*!
    This class constructs the staple.
    While the originial version was written by J.Noaki,
    the present version is completely modified by H.Matsufuru
    except for the interface.
                                    [28 Dec 2011 H.Matsufuru]
    Thread-parallelized.
    void version of upper and lower functions added; these are
    faster than the versions returning Field_G object.
                                    [28 Sep 2013 H.Matsufuru]
    Add parameters for output.      [27 Jun 2016 Y.Namekawa]
    Factory is introduced.          [24 Jan 2017 Y.Namekawa]
 */

class Staple_lex : public Staple
{
 public:
  static const std::string class_name;

 private:
  std::string m_filename_output;

  ShiftField_lex *m_shift;
  Field_G        m_staple, m_v, m_w; //!< temporary fields.

 public:
  Staple_lex() : Staple()
  {
    m_filename_output = "stdout";
    m_shift           = new ShiftField_lex;
  }

  ~Staple_lex()
  {
    delete m_shift;
  }

  //! setting parameters.
  void set_parameters(const Parameters& params);

  //! constructs upper staple in mu-nu plane.
  void upper(Field_G&, const Field_G&, const int mu, const int nu);

  //! constructs lower staple in mu-nu plane.
  void lower(Field_G&, const Field_G&, const int mu, const int nu);

  //! constructs staple in mu-direction (summing up nu-direction).
  void staple(Field_G&, const Field_G&, const int mu);

  //! calculates plaquette value.
  double plaquette(const Field_G&);

  //! calculates spatial plaquette value.
  double plaq_s(const Field_G&);

  //! calculates temporal plaquette value.
  double plaq_t(const Field_G&);
};
#endif
