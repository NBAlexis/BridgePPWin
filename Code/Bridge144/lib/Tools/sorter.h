/*!
        @file    $Id: sorter.h #$

        @brief

        @author  Tatsumi Aoyama  (aoym)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 00:49:46 #$

        @version $LastChangedRevision: 1561 $
*/

#ifndef SORTER_INCLUDED
#define SORTER_INCLUDED

#include <cmath>
#include <algorithm>
#include <cassert>
#include <vector>
#include <string>

//! @class Sorter
//    sort utility class.
//
//  usage:
//    Sorter sorter(sort_type);
//      sort_type =
//        abs_ascending  : absolute ascending order   |v_0| <= |v_1| <= ...
//        abs_descending : absolute descending order  |v_0| >= |v_1| >= ...
//        ascending      : arithmetic ascending order  v_0 <= v_1 <= ...
//        descending     : arithmetic descending order v_0 >= v_1 >= ...
//    sorter.sort(vec);  to sort vec : vector<double>. destructive.
//    sorter.sort_index(vec);  returns index of vec in sorted order.
//    sorter.comp(a, b);  compare two values a, b with the sort condition.
//

class Sorter
{
 public:

  Sorter(const std::string& type);

  virtual ~Sorter();

  //! sort an array of values; v is sorted on exit.
  void sort(std::vector<double>& v);

  //! sort first n elements of an array of values; v is sorted on exit.
  void sort(std::vector<double>& v, const size_t nelem);

  //! sort an array and return list of index; v is sorted on exit.
  std::vector<int> sort_index(std::vector<double>& v);

  //! sort an array and return list of index. first n elemnts are considered.
  std::vector<int> sort_index(std::vector<double>& v, const size_t nelem);

  //! call sort condition.
  bool comp(const double lhs, const double rhs);

 private:

  // sort order
  struct by_abs_ascend;
  struct by_abs_descend;
  struct by_ascend;
  struct by_descend;

  struct by_order;
  struct proxy;

  typedef std::pair<int, double>   pair_t;

  by_order *m_order;

  // non-copyable
  Sorter(const Sorter&);
  Sorter& operator=(const Sorter&);
};
#endif /* SORTER_INCLUDED */
