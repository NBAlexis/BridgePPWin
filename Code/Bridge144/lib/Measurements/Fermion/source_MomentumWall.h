/*!
        @file    $Id:: source_MomentumWall.h #$

        @brief   Momentum wall source

        @author  Noriyoshi Ishii (ishii)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 18:35:38 +0900 #$

        @version $LastChangedRevision: 1571 $
*/

#ifndef SOURCE_MOMENTUM_WALL_INCLUDED
#define SOURCE_MOMENTUM_WALL_INCLUDED

#include "source.h"
#include "Field/index_lex.h"

#include "IO/bridgeIO.h"
using Bridge::vout;


//! Momentum wall source for 4-spinor fermion.

/*!
    This class sets a momentum wall source vector
    for the 4-spinor (Wilson-type) fermion.
                          [29 Jan 2014 N.Ishii]
    In order to commit to Bridge++ main stream,
    comments etc... are modified.
                          [19 Feb 2014 S.Ueda]
 */

class Source_MomentumWall : public Source
{
 public:
  static const std::string class_name;

 public:
  Source_MomentumWall() {}

  void set_parameters(const Parameters& params);
  void set_parameters(const std::vector<int>& source_position,
                      const std::vector<int>& source_momentum);

  void set(Field& v, int j);

 private:
  Index_lex        m_index;
  std::vector<int> m_source_position;
  std::vector<int> m_source_momentum;
  bool             m_in_node;
};
#endif /* SOURCE_MOMENTUM_WALL_INCLUDED */
