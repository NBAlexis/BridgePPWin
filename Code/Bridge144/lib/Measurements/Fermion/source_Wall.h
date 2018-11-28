/*!
        @file    $Id:: source_Wall.h #$

        @brief

        @author  Yusuke Namekawa (namekawa)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 18:35:38 #$

        @version $LastChangedRevision: 1571 $
*/

#ifndef SOURCE_WALL_INCLUDED
#define SOURCE_WALL_INCLUDED

#include "source.h"
#include "Field/index_lex.h"

#include "IO/bridgeIO.h"
using Bridge::vout;


//! Wall source for 4-spinor fermion.

/*!
    This class sets a wall source vector
    for the 4-spinor (Wilson-type) fermion.
                          [02 Feb 2013 Y.Namekawa]
 */

class Source_Wall : public Source
{
 public:
  static const std::string class_name;

 public:
  Source_Wall() {}

  void set_parameters(const Parameters& params);
  void set_parameters(const std::vector<int>& source_position);

  void set(Field& v, int j);

 private:
  Index_lex        m_index;
  std::vector<int> m_source_position;
  bool             m_in_node;
};
#endif /* SOURCE_WALL_INCLUDED */
