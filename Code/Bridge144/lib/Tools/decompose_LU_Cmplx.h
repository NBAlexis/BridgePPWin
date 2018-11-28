/*!
        @file    $Id:: decompose_LU_Cmplx.h #$

        @brief

        @author  Satoru Ueda  (sueda)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 00:49:46 #$

        @version $LastChangedRevision: 1561 $
*/

#ifndef DECOMPOSE_LU_CMPLX_INCLUDED
#define DECOMPOSE_LU_CMPLX_INCLUDED

#include <cfloat>
#include <cmath>
#include <valarray>

#include "bridge_complex.h"

#include "IO/bridgeIO.h"

class Decompose_LU_Cmplx
{
 public:
  Decompose_LU_Cmplx(size_t N) : N((int)N), N2(2 * (int)N), size((int)N * N2), m_lu(size) {}

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
