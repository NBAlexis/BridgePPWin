/*!
        @file    decompose_QR_Cmplx.h

        @brief

        @author  Satoru Ueda  (sueda)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/

#ifndef DECOMPOSE_QR_CMPLX_INCLUDED
#define DECOMPOSE_QR_CMPLX_INCLUDED

#include <cmath>
#include <cfloat>
#include <valarray>

/*
  class BAPI for QR decomposition with Householder transformation
  The notation of Householder transformation is the same as that of GSL.
  H = I - tau vv^\dag
  H^\dag x = - |x| sign(Re(x_1)) e_1
 */

class BAPI Decompose_QR_Cmplx
{
 public:
  Decompose_QR_Cmplx(const size_t& Nin) : N(Nin), N2(2 * N), size(N * N2),
    m_qr(size), m_tau(N2)
  {
  }

  void set_matrix(const double *mat);
  void get_R(double *mat);
  void solve(double *vec);
  void get_inverse(double *mat);
  void mult_inverse(double *mat);
  void get_Q(double *mat);
  void get_Qu(double *mat);
  void mult_Q(double *mat);

 private:
  size_t N;
  size_t N2;
  size_t size;

  // store R and householder vector
  std::valarray<double> m_qr;
  std::valarray<double> m_tau;

  // matrix index [i,j]
  inline size_t re(const size_t i, const size_t j)
  {
    return N2 * i + 2 * j;
  }

  inline size_t im(const size_t i, const size_t j)
  {
    return N2 * i + 2 * j + 1;
  }

  // vector index [i]
  inline size_t re(const size_t i)
  {
    return 2 * i;
  }

  inline size_t im(const size_t i)
  {
    return 2 * i + 1;
  }
};
#endif
