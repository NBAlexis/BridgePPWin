/*!
        @file    $Id:: source_Local.h #$

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 18:35:38 #$

        @version $LastChangedRevision: 1571 $
*/

#ifndef SOURCE_LOCAL_INCLUDED
#define SOURCE_LOCAL_INCLUDED

#include "source.h"
#include "Field/index_lex.h"

#include "IO/bridgeIO.h"
using Bridge::vout;


//! Local source for 4-spinor fermion.

/*!
    This class sets an local source vector for the 4-spinor
    (Wilson-type) fermion.
                                     [19 Feb 2012 H.Matsufuru]
    (Coding history will be recovered from trac.)
    YAML is implemented.             [14 Nov 2012 Y.Namekawa]
 */

class Source_Local : public Source {
 public:
  static const std::string class_name;

 public:

  Source_Local()
    : Source() {}

  void set_parameters(const Parameters& params);
  void set_parameters(const std::vector<int>& source_position);

  void set(Field& v, int j);

 private:
  Index_lex        m_index;   // lexical only.
  std::vector<int> m_source_position;
  bool             m_in_node;
};
#endif /* SOURCE_LOCAL_INCLUDED */
