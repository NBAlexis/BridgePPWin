/*!
        @file    source_Exponential.h

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/

#ifndef SOURCE_EXPONENTIAL_INCLUDED
#define SOURCE_EXPONENTIAL_INCLUDED

#include "source.h"

#include "Field/index_lex.h"

#include "IO/bridgeIO.h"
using Bridge::vout;


//! Exponentially smeared source for 4-spinor fermion.

/*!
    This class BAPI sets an exponentially smeared source vector
    for the 4-spinor (Wilson-element_type) fermion.
                                      [19 Feb 2012 H.Matsufuru]
    (Coding history will be recovered from trac.)
    YAML is implemented.              [14 Nov 2012 Y.Namekawa]
    Add set_all_color_spin,etc        [ 4 Apr 2017 Y.Namekawa]
 */

class BAPI Source_Exponential : public Source
{
 public:
  static const std::string class_name;

 private:
  Index_lex m_index;
  std::vector<int> m_source_position;
  double m_slope, m_power;
  bool m_in_node;
  Field m_src_func;

 public:
  Source_Exponential() {}

  void set_parameters(const Parameters& params);
  void set_parameters(const std::vector<int>& source_position,
                      const double slope, const double power);

  void set(Field& v, const int idx);
  void set(Field& v, const int i_color, const int i_spin);
  void set_all_color(Field& v, const int i_spin);
  void set_all_color_spin(Field& v);

#ifdef USE_FACTORY
 private:
  static Source *create_object()
  {
    return new Source_Exponential();
  }

 public:
  static bool register_factory()
  {
    return Source::Factory::Register("Exponential", create_object);
  }
#endif
};
#endif /* SOURCE_EXPONENTIAL_INCLUDED */
