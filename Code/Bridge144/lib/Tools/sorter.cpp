#include "BridgeLib_Private.h"

/*!
        @file    $Id:: sorter.cpp #$

        @brief

        @author  Tatsumi Aoyama  (aoym)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 00:49:46 #$

        @version $LastChangedRevision: 1561 $
*/

#include "sorter.h"

//====================================================================
//! base class for sort ordering
struct Sorter::by_order
{
  virtual ~by_order() {}

  virtual
  bool operator()(const double lhs, const double rhs) const = 0;

  bool operator()(const pair_t& lhs, const pair_t& rhs) const
  {
    return operator()(lhs.second, rhs.second);
  }
};

//====================================================================
//! absolute ascending order
//  i.e. |v_0| <= |v_1| <= ... <= |v_{n-1}|
struct Sorter::by_abs_ascend : public Sorter::by_order
{
  bool operator()(const double lhs, const double rhs) const
  { return fabs(lhs) < fabs(rhs); }
};

//====================================================================
//! absolute descending order
//  i.e. |v_0| >= |v_1| >= ... >= |v_{n-1}|
struct Sorter::by_abs_descend : public Sorter::by_order
{
  bool operator()(const double lhs, const double rhs) const
  {
    return fabs(lhs) > fabs(rhs);
  }
};

//====================================================================
//! ascending order
//  i.e. v_0 <= v_1 <= ... <= v_{n-1}
struct Sorter::by_ascend : public Sorter::by_order
{
  virtual
  bool operator()(const double lhs, const double rhs) const
  {
    return lhs < rhs;
  }
};

//====================================================================
//! descending order
//  i.e. v_0 >= v_1 >= ... >= v_{n-1}
struct Sorter::by_descend : public Sorter::by_order
{
  virtual
  bool operator()(const double lhs, const double rhs) const
  {
    return lhs > rhs;
  }
};

//====================================================================
//! proxy object to pass to stl sort algorithm
struct Sorter::proxy
{
 private:
  const Sorter::by_order& m_order;

 public:
  proxy(const Sorter::by_order& order) : m_order(order) {}

  bool operator()(const double lhs, const double rhs) const
  {
    return m_order(lhs, rhs);
  }

  bool operator()(const pair_t& lhs, const pair_t& rhs) const
  {
    return m_order(lhs, rhs);
  }
};

//====================================================================
//! constructor with sort ordering as a string arg
Sorter::Sorter(const std::string& type) : m_order(0)
{
  if (type == "abs_ascending") m_order = new by_abs_ascend;
  if (type == "abs_descending") m_order = new by_abs_descend;
  if (type == "ascending") m_order = new by_ascend;
  if (type == "descending") m_order = new by_descend;
}


//====================================================================
//! destructor
Sorter::~Sorter()
{
  if (m_order) delete m_order;
}


//====================================================================
//! sort an array of values; v is sorted on exit.
void Sorter::sort(std::vector<double>& v)
{
  assert(m_order != NULL);

  return std::sort(v.begin(), v.end(), proxy(*m_order));
}


//====================================================================
//! sort first n elments of an array of values; v is sorted on exit.
void Sorter::sort(std::vector<double>& v, const size_t nelem)
{
  assert(m_order != NULL);

  return std::sort(v.begin(), v.begin() + nelem, proxy(*m_order));
}


//====================================================================
//! sort an array and return list of index.
std::vector<int> Sorter::sort_index(std::vector<double>& v)
{
  assert(m_order != NULL);

  std::vector<std::pair<int, double> > w(v.size());

  for (size_t i = 0, n = v.size(); i < n; ++i) {
    w[i].first  = i;
    w[i].second = v[i];
  }

  std::sort(w.begin(), w.end(), proxy(*m_order));

  std::vector<int> idx(v.size());

  for (size_t i = 0, n = v.size(); i < n; ++i) {
    idx[i] = w[i].first;
    v[i]   = w[i].second;
  }

  return idx;
}


//====================================================================
//! sort an array and return list of index. first n elements are condidered.
std::vector<int> Sorter::sort_index(std::vector<double>& v, const size_t nelem)
{
  assert(m_order != NULL);

  std::vector<std::pair<int, double> > w(v.size());

  for (size_t i = 0, n = v.size(); i < n; ++i) {
    w[i].first  = i;
    w[i].second = v[i];
  }

  std::sort(w.begin(), w.begin() + nelem, proxy(*m_order));

  std::vector<int> idx(v.size());

  for (size_t i = 0, n = v.size(); i < n; ++i) {
    idx[i] = w[i].first;
    v[i]   = w[i].second;
  }

  return idx;
}


//====================================================================
//! call sort condition.
bool Sorter::comp(const double lhs, const double rhs)
{
  assert(m_order != NULL);
  return m_order->operator()(lhs, rhs);
}


//==========================================================
//==================================================END=====
