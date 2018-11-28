/*!
        @file    $Id:: source_Exponential.h #$

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 18:35:38 #$

        @version $LastChangedRevision: 1571 $
*/

#ifndef SOURCE_EXPONENTIAL_INCLUDED
#define SOURCE_EXPONENTIAL_INCLUDED

#include "source.h"
#include "Field/index_lex.h"

#include "IO/bridgeIO.h"
using Bridge::vout;


//! Exponentially smeared source for 4-spinor fermion.

/*!
    This class sets an exponentially smeared source vector
    for the 4-spinor (Wilson-type) fermion.
                                      [19 Feb 2012 H.Matsufuru]
    (Coding history will be recovered from trac.)
    YAML is implemented.              [14 Nov 2012 Y.Namekawa]
 */

class Source_Exponential : public Source
{
 public:
  static const std::string class_name;

 public:

  Source_Exponential() {}

  void set_parameters(const Parameters& params);
  void set_parameters(const std::vector<int>& source_position,
                      const double slope, const double power);

  void set(Field& v, int j);

 private:
  Index_lex        m_index;
  std::vector<int> m_source_position;
  double           m_slope, m_power;
  bool             m_in_node;
  Field            m_src_func;
};
#endif /* SOURCE_EXPONENTIAL_INCLUDED */
