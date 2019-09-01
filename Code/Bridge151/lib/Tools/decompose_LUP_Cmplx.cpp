/*!
        @file    decompose_LUP_Cmplx.cpp

        @brief

        @author  Satoru Ueda  (sueda)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/
#include "BridgeLib_Private.h"
#include "decompose_LUP_Cmplx.h"

//====================================================================
void Decompose_LUP_Cmplx::set_matrix(const double *mat)
{
  for (int i = 0; i < size; ++i) {
    m_lu[i] = mat[i];
  }

  // scale factor
  std::valarray<double> scale(N);
  for (int i = 0; i < N; ++i) {
    scale[i] = 1;
  }
  double epsilon = DBL_EPSILON;
  m_sign = 1;
  // define scale;
  // abort singlar matrix
  for (int i = 0; i < N; ++i) {
    double max_norm = 0;
    for (int j = 0; j < N; ++j) {
      double norm2 = m_lu[re(i, j)] * m_lu[re(i, j)] + m_lu[im(i, j)] * m_lu[im(i, j)];
      if (norm2 > max_norm) {
        max_norm = norm2;
      }
      if (max_norm < epsilon) {
        // assert
      }
      scale[i] = max_norm;
    }
  }

  for (int r = 0; r < N - 1; ++r) {
    // pivot
    double max = 0;
    for (int i = r; i < N; ++i) {
      double norm_ir = m_lu[re(i, r)] * m_lu[re(i, r)] + m_lu[im(i, r)] * m_lu[im(i, r)];
      if (norm_ir > scale[i] * max) {
        max        = norm_ir;
        m_pivot[r] = i;
      }
    }

    // pivot (a_{rj} <-> a_{pj})
    if (m_pivot[r] != r) {
      int p = m_pivot[r];
      for (int j = 0; j < N; ++j) {
        std::swap(m_lu[re(r, j)], m_lu[re(p, j)]);
        std::swap(m_lu[im(r, j)], m_lu[im(p, j)]);
      }
      std::swap(scale[r], scale[p]);
      m_sign *= -1;
    }
    double inv_rr_re = m_lu[re(r, r)] / max;
    double inv_rr_im = -m_lu[im(r, r)] / max;
    // a_{ir} /= a_{rr}
    // a_{ij} -= a_{ir} a_{rj} (j > r)
    for (int i = r + 1; i < N; ++i) {
      double ir_re = m_lu[re(i, r)];
      double ir_im = m_lu[im(i, r)];
      m_lu[re(i, r)] = ir_re * inv_rr_re - ir_im * inv_rr_im;
      m_lu[im(i, r)] = ir_re * inv_rr_im + ir_im * inv_rr_re;
      ir_re          = m_lu[re(i, r)];
      ir_im          = m_lu[im(i, r)];

      for (int j = r + 1; j < N; ++j) {
        double& rj_re = m_lu[re(r, j)];
        double& rj_im = m_lu[im(r, j)];
        m_lu[re(i, j)] -= ir_re * rj_re - ir_im * rj_im;
        m_lu[im(i, j)] -= ir_re * rj_im + ir_im * rj_re;
      }
    }
  }

  m_pivot[N - 1] = N - 1;
}


//====================================================================
void Decompose_LUP_Cmplx::solve(double *vec)
{
  // Pivot
  for (int i = 0; i < N; ++i) {
    std::swap(vec[re(i)], vec[re(m_pivot[i])]);
    std::swap(vec[im(i)], vec[im(m_pivot[i])]);
  }

  // Forward substitution
  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < i; ++j) {
      vec[re(i)] -= m_lu[re(i, j)] * vec[re(j)] - m_lu[im(i, j)] * vec[im(j)];
      vec[im(i)] -= m_lu[re(i, j)] * vec[im(j)] + m_lu[im(i, j)] * vec[re(j)];
    }
  }

  // Back substitution
  for (int i = N - 1; i >= 0; --i) {
    for (int j = i + 1; j < N; ++j) {
      vec[re(i)] -= m_lu[re(i, j)] * vec[re(j)] - m_lu[im(i, j)] * vec[im(j)];
      vec[im(i)] -= m_lu[re(i, j)] * vec[im(j)] + m_lu[im(i, j)] * vec[re(j)];
    }
    double inv_norm = 1 / (m_lu[re(i, i)] * m_lu[re(i, i)] + m_lu[im(i, i)] * m_lu[im(i, i)]);
    double real     = vec[re(i)];
    double imag     = vec[im(i)];

    vec[re(i)] = (real * m_lu[re(i, i)] + imag * m_lu[im(i, i)]) * inv_norm;
    vec[im(i)] = (imag * m_lu[re(i, i)] - real * m_lu[im(i, i)]) * inv_norm;
  }
}


//====================================================================
void Decompose_LUP_Cmplx::get_inverse(double *mat)
{
  for (int i = 0; i < size; ++i) {
    mat[i] = 0;
  }
  for (int i = 0; i < N; ++i) {
    mat[re(i, i)] = 1;
  }

  mult_inverse(mat);
}


//====================================================================
void Decompose_LUP_Cmplx::mult_inverse(double *mat)
{
  // Pivot
  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < N; ++j) {
      std::swap(mat[re(i, j)], mat[re(m_pivot[i], j)]);
      std::swap(mat[im(i, j)], mat[im(m_pivot[i], j)]);
    }
  }

  // Forward substitution
  for (int i = 0; i < N; ++i) {
    for (int r = 0; r < i; ++r) {
      double& luir_re = m_lu[re(i, r)];
      double& luir_im = m_lu[im(i, r)];
      for (int j = 0; j < N; ++j) {
        mat[re(i, j)] -= luir_re * mat[re(r, j)] - luir_im * mat[im(r, j)];
        mat[im(i, j)] -= luir_re * mat[im(r, j)] + luir_im * mat[re(r, j)];
      }
    }
  }

  // Back substitution
  for (int i = N - 1; i >= 0; --i) {
    for (int r = i + 1; r < N; ++r) {
      double& luir_re = m_lu[re(i, r)];
      double& luir_im = m_lu[im(i, r)];
      for (int j = 0; j < N; ++j) {
        mat[re(i, j)] -= luir_re * mat[re(r, j)] - luir_im * mat[im(r, j)];
        mat[im(i, j)] -= luir_re * mat[im(r, j)] + luir_im * mat[re(r, j)];
      }
    }
    for (int j = 0; j < N; ++j) {
      double inv  = 1 / (m_lu[re(i, i)] * m_lu[re(i, i)] + m_lu[im(i, i)] * m_lu[im(i, i)]);
      double real = mat[re(i, j)];
      double imag = mat[im(i, j)];

      mat[re(i, j)] = (real * m_lu[re(i, i)] + imag * m_lu[im(i, i)]) * inv;
      mat[im(i, j)] = (imag * m_lu[re(i, i)] - real * m_lu[im(i, i)]) * inv;
    }
  }
}


//====================================================================
dcomplex Decompose_LUP_Cmplx::determinant()
{
  double real = m_sign;
  double imag = 0;

  for (int i = 0; i < N; ++i) {
    double old_re = real;
    double old_im = imag;
    real = old_re * m_lu[re(i, i)] - old_im * m_lu[im(i, i)];
    imag = old_re * m_lu[im(i, i)] + old_im * m_lu[re(i, i)];
  }

  return cmplx(real, imag);
}


//====================================================================

/*
void Decompose_LU_Cmplx::copy_LU(double* l, double* u)
{
  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < i; ++j) {
      l[re(i,j)] = m_lu[re(i,j)];
      l[im(i,j)] = m_lu[im(i,j)];
      u[re(i,j)] = 0;
      u[im(i,j)] = 0;
    }
    l[re(i,i)] = 1;
    l[im(i,i)] = 0;
    u[re(i,i)] = m_lu[re(i,i)];
    u[im(i,i)] = m_lu[im(i,i)];
    for (int j = i+1; j < N; ++j) {
      l[re(i,j)] = 0;
      l[im(i,j)] = 0;
      u[re(i,j)] = m_lu[re(i,j)];
      u[im(i,j)] = m_lu[im(i,j)];
    }
  }

}
*/

//====================================================================
