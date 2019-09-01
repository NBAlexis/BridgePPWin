/*!
        @file    source_Wall.h

        @brief

        @author  Yusuke Namekawa (namekawa)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/

#ifndef SOURCE_WALL_INCLUDED
#define SOURCE_WALL_INCLUDED

#include "source.h"
#include "Field/index_lex.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//! Wall source for 4-spinor fermion.

/*!
    This class BAPI sets a wall source vector
    for the 4-spinor (Wilson-element_type) fermion.
                                [02 Feb 2013 Y.Namekawa]
    Add set_all_color_spin,etc  [ 4 Apr 2017 Y.Namekawa]
 */

class BAPI Source_Wall : public Source
{
 public:
  static const std::string class_name;

 private:
  Index_lex m_index;
  std::vector<int> m_source_position;
  bool m_in_node;

 public:
  Source_Wall() {}

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
    return new Source_Wall();
  }

 public:
  static bool register_factory()
  {
    return Source::Factory::Register("Wall", create_object);
  }
#endif
};
#endif /* SOURCE_WALL_INCLUDED */
