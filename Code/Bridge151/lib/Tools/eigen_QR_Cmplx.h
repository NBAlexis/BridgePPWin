/*!
        @file    eigen_QR_Cmplx.h

        @brief

        @author  Satoru Ueda  (sueda)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/

#ifndef EIGEN_QR_CMPLX_INCLUDED
#define EIGEN_QR_CMPLX_INCLUDED

#include <cfloat>
#include <valarray>
#include <string>

#include "decompose_Hessenberg_Cmplx.h"

class BAPI Eigen_QR_Cmplx
{
 public:
  static const std::string class_name;

 public:
  Eigen_QR_Cmplx(const size_t N) : N(N), N2(2 * N), size(N2 * N),
    m_mat(size), m_q(size)
  {
  }

  std::valarray<double> solve(const double *matrix);
  void get_R(double *r);
  void get_Q(double *q);

 private:
  size_t N;
  size_t N2;
  size_t size;

  std::valarray<double> m_mat;
  std::valarray<double> m_q;
  inline size_t re(const size_t i, const size_t j)
  {
    return N2 * i + 2 * j;
  }

  inline size_t im(const size_t i, const size_t j)
  {
    return N2 * i + 2 * j + 1;
  }

  inline size_t re(const size_t i)
  {
    return 2 * i;
  }

  inline size_t im(const size_t i)
  {
    return 2 * i + 1;
  }

  void qr_step(const int rank);
};
#endif
