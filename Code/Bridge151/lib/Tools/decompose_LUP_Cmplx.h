/*!
        @file    decompose_LUP_Cmplx.h

        @brief

        @author  Satoru Ueda  (sueda)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/

#ifndef DECOMPOSE_LUP_CMPLX_INCLUDED
#define DECOMPOSE_LUP_CMPLX_INCLUDED

#include <valarray>
#include <cfloat>
#include <cmath>

#include "bridge_complex.h"

#include "IO/bridgeIO.h"

class BAPI Decompose_LUP_Cmplx
{
 public:
  Decompose_LUP_Cmplx(size_t N) : N(static_cast<int>(N)), N2(2 * static_cast<int>(N)), size(static_cast<int>(N) * N2),
    m_lu(size), m_pivot(N) {}

  void set_matrix(const double *mat);

  // solve Ax = b: vec -> A^{-1} vec
  void solve(double *vec);

  // M -> A^{-1}
  void get_inverse(double *mat_inv);

  // M -> A^{-1} * M
  void mult_inverse(double *mat);

  // return det(A)
  dcomplex determinant();

 private:
  int N;
  int N2;
  int size;
  std::valarray<double> m_lu;
  // pivot index
  std::valarray<int> m_pivot;
  // # of pivot
  // even ->  1;
  // odd  -> -1;
  int m_sign;

  inline size_t re(int i, int j)
  {
    return N2 * i + 2 * j;
  }

  inline size_t im(int i, int j)
  {
    return N2 * i + 2 * j + 1;
  }

  inline size_t re(int i)
  {
    return 2 * i;
  }

  inline size_t im(int i)
  {
    return 2 * i + 1;
  }
};
#endif
