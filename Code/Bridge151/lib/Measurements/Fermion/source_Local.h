/*!
        @file    source_Local.h

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/

#ifndef SOURCE_LOCAL_INCLUDED
#define SOURCE_LOCAL_INCLUDED

#include "source.h"
#include "Field/index_lex.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//! Local source for 4-spinor fermion.

/*!
    This class BAPI sets an local source vector for the 4-spinor
    (Wilson-element_type) fermion.
                                     [19 Feb 2012 H.Matsufuru]
    (Coding history will be recovered from trac.)
    YAML is implemented.             [14 Nov 2012 Y.Namekawa]
    Add set_all_color_spin,etc       [ 4 Apr 2017 Y.Namekawa]
 */

class BAPI Source_Local : public Source {
 public:
  static const std::string class_name;

 private:
  Index_lex m_index;          // lexical only.
  std::vector<int> m_source_position;
  bool m_in_node;

 public:
  Source_Local()
    : Source() {}

  void set_parameters(const Parameters& params);
  void set_parameters(const std::vector<int>& source_position);

  void set(Field& v, const int idx);
  void set(Field& v, const int i_color, const int i_spin);
  void set_all_color(Field& v, const int i_spin);
  void set_all_color_spin(Field& v);

#ifdef USE_FACTORY
 private:
  static Source *create_object()
  {
    return new Source_Local();
  }

 public:
  static bool register_factory()
  {
    return Source::Factory::Register("Local", create_object);
  }
#endif
};
#endif /* SOURCE_LOCAL_INCLUDED */
