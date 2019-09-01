/*!
        @file    epsilonTensor.h

        @brief

        @author  Yusuke Namekawa (namekawa)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/

#ifndef EPSILONTENSOR_INCLUDED
#define EPSILONTENSOR_INCLUDED

#include <vector>

//! Epsilon tensor utility

/*!
    EpsilonTensor gives index and value of a totally antisymmetric
    tensor, exported from corr2pt_4spinor. It becomes const, rather
    than functions, thanks to Aoyama-san.
                                   [20 Apr 2017 Y.Namekawa]
 */

/*
namespace EpsilonTensor
{
  const int m_epsilon_3_index[][3] = {
    { 0, 1, 2 },
    { 1, 2, 0 },
    { 2, 0, 1 },
    { 2, 1, 0 },
    { 1, 0, 2 },
    { 0, 2, 1 },
  };

  const int m_epsilon_3_value[] = {
     1,
     1,
     1,
    -1,
    -1,
    -1,
  };
}
*/

class BAPI EpsilonTensor
{
  // public:
  //   static const std::string class_name;
  //
  // protected:
  //   Bridge::VerboseLevel m_vl;

 private:
  int m_Nepsilon_3;
  int m_Nepsilon_3_index;

  std::vector<std::vector<int> > m_epsilon_3_index;
  std::vector<int> m_epsilon_3_value;

 public:
  EpsilonTensor()
  {
    init();
  }

  virtual ~EpsilonTensor() {}

 private:
  // non-copyable
  EpsilonTensor(const EpsilonTensor&);
  EpsilonTensor& operator=(const EpsilonTensor&);

 public:
  int epsilon_3_index(const int n, const int i) const;
  int epsilon_3_value(const int n) const;

 private:
  void init();
};
#endif
