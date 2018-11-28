#include "BridgeLib_Private.h"

/*!
        @file    $Id:: decompose_QR_Cmplx.cpp #$

        @brief

        @author  Satoru Ueda  (sueda)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 00:49:46 #$

        @version $LastChangedRevision: 1561 $
*/

#include "decompose_QR_Cmplx.h"

//====================================================================
void Decompose_QR_Cmplx::set_matrix(const double *mat)
{
  // double epsilon = DBL_EPSILON;

  // matrix copy
  for (int i = 0; i < size; ++i) {
    m_qr[i] = mat[i];
  }

  // Householder trans for r col.
  for (int r = 0; r < N; ++r) {
    // diagnal part
    double& alpha_re = m_qr[re(r, r)];
    double& alpha_im = m_qr[im(r, r)];

    // make Householder vector v and factor tau.
    // H = I - tau vv^\dag
    // H^\dag x = (beta, 0, ...)^T
    // nu = |x| sign(Re(x_1))

    // tau = sigma / nu
    // beta = - nu
    // Re(sigma) = Re(x_1) + nu
    // Im(sigma) = Im(x_1)
    // v = (1, x_2/sigma, ..., x_n/sigma);


    double nu = 0;
    for (int i = r; i < N; ++i) {
      nu += m_qr[re(i, r)] * m_qr[re(i, r)] + m_qr[im(i, r)] * m_qr[im(i, r)];
    }

    /*
    if (nu < epsilon) {
      // assertion?
    }
    */

    nu           = std::sqrt(nu) * ((alpha_re < 0) ? -1 : 1);
    m_tau[re(r)] = alpha_re + nu;
    m_tau[im(r)] = alpha_im;

    double norm2_sigma  = m_tau[re(r)] * m_tau[re(r)] + m_tau[im(r)] * m_tau[im(r)];
    double inv_sigma_re = m_tau[re(r)] / norm2_sigma;
    double inv_sigma_im = -m_tau[im(r)] / norm2_sigma;

    m_tau[re(r)] /= nu;
    m_tau[im(r)] /= nu;

    // m_qr[r,r] = -nu
    // m_qr[r,i] /= sigma
    alpha_re = -nu;
    alpha_im = 0;
    for (int i = r + 1; i < N; ++i) {
      double ir_re = m_qr[re(i, r)];
      double ir_im = m_qr[im(i, r)];

      m_qr[re(i, r)] = ir_re * inv_sigma_re - ir_im * inv_sigma_im;
      m_qr[im(i, r)] = ir_im * inv_sigma_re + ir_re * inv_sigma_im;
    }

    // H^\dag a_j
    for (int j = r + 1; j < N; ++j) {
      // va = <v|a_j>
      double va_re = m_qr[re(r, j)];
      double va_im = m_qr[im(r, j)];
      for (int i = r + 1; i < N; ++i) {
        va_re += m_qr[re(i, r)] * m_qr[re(i, j)] + m_qr[im(i, r)] * m_qr[im(i, j)];
        va_im += m_qr[re(i, r)] * m_qr[im(i, j)] - m_qr[im(i, r)] * m_qr[re(i, j)];
      }
      // beta = tau^* * va
      // a_j -= bata * v
      double beta_re = m_tau[re(r)] * va_re + m_tau[im(r)] * va_im;
      double beta_im = m_tau[re(r)] * va_im - m_tau[im(r)] * va_re;
      m_qr[re(r, j)] -= beta_re;
      m_qr[im(r, j)] -= beta_im;
      for (int i = r + 1; i < N; ++i) {
        double vi_re = m_qr[re(i, r)];
        double vi_im = m_qr[im(i, r)];

        m_qr[re(i, j)] -= beta_re * vi_re - beta_im * vi_im;
        m_qr[im(i, j)] -= beta_re * vi_im + beta_im * vi_re;
      }
    }
  }
}


//====================================================================
void Decompose_QR_Cmplx::get_R(double *mat)
{
  for (int i = 0; i < size; ++i) {
    mat[i] = m_qr[i];
  }
  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < i; ++j) {
      mat[re(i, j)] = mat[im(i, j)] = 0;
    }
  }
}


//====================================================================
void Decompose_QR_Cmplx::solve(double *vec)
{
  // vec -> Q^\dag vec
  for (int r = 0; r < N; ++r) {
    // va = <v|vec>
    double va_re = vec[re(r)];
    double va_im = vec[im(r)];
    for (int i = r + 1; i < N; ++i) {
      va_re += m_qr[re(i, r)] * vec[re(i)] + m_qr[im(i, r)] * vec[im(i)];
      va_im += m_qr[re(i, r)] * vec[im(i)] - m_qr[im(i, r)] * vec[re(i)];
    }

    // beta = tau^* * va
    // a_j -= bata * v
    double beta_re = m_tau[re(r)] * va_re + m_tau[im(r)] * va_im;
    double beta_im = m_tau[re(r)] * va_im - m_tau[im(r)] * va_re;
    vec[re(r)] -= beta_re;
    vec[im(r)] -= beta_im;
    for (int i = r + 1; i < N; ++i) {
      double vi_re = m_qr[re(i, r)];
      double vi_im = m_qr[im(i, r)];

      vec[re(i)] -= beta_re * vi_re - beta_im * vi_im;
      vec[im(i)] -= beta_re * vi_im + beta_im * vi_re;
    }
  }

  // vec -> R^{-1} vec
  // Back substitution
  for (int i = (int)N - 1; i >= 0; --i) {
    for (int j = i + 1; j < N; ++j) {
      vec[re(i)] -= m_qr[re(i, j)] * vec[re(j)] - m_qr[im(i, j)] * vec[im(j)];
      vec[im(i)] -= m_qr[re(i, j)] * vec[im(j)] + m_qr[im(i, j)] * vec[re(j)];
    }
    double inv = 1 / m_qr[re(i, i)];

    vec[re(i)] *= inv;
    vec[im(i)] *= inv;
  }
}


//====================================================================
void Decompose_QR_Cmplx::get_inverse(double *mat)
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
void Decompose_QR_Cmplx::mult_inverse(double *mat)
{
  // H^\dag a_j
  for (int r = 0; r < N; ++r) {
    for (int j = 0; j < N; ++j) {
      // va = <v|a_j>
      double va_re = mat[re(r, j)];
      double va_im = mat[im(r, j)];
      for (int i = r + 1; i < N; ++i) {
        va_re += m_qr[re(i, r)] * mat[re(i, j)] + m_qr[im(i, r)] * mat[im(i, j)];
        va_im += m_qr[re(i, r)] * mat[im(i, j)] - m_qr[im(i, r)] * mat[re(i, j)];
      }

      // beta = tau^* * va
      // a_j -= bata * v
      double beta_re = m_tau[re(r)] * va_re + m_tau[im(r)] * va_im;
      double beta_im = m_tau[re(r)] * va_im - m_tau[im(r)] * va_re;
      mat[re(r, j)] -= beta_re;
      mat[im(r, j)] -= beta_im;
      for (int i = r + 1; i < N; ++i) {
        double vi_re = m_qr[re(i, r)];
        double vi_im = m_qr[im(i, r)];

        mat[re(i, j)] -= beta_re * vi_re - beta_im * vi_im;
        mat[im(i, j)] -= beta_re * vi_im + beta_im * vi_re;
      }
    }
  }

  // Back substitution
  for (int i = (int)N - 1; i >= 0; --i) {
    for (int r = i + 1; r < N; ++r) {
      double& qrir_re = m_qr[re(i, r)];
      double& qrir_im = m_qr[im(i, r)];
      for (int j = 0; j < N; ++j) {
        mat[re(i, j)] -= qrir_re * mat[re(r, j)] - qrir_im * mat[im(r, j)];
        mat[im(i, j)] -= qrir_re * mat[im(r, j)] + qrir_im * mat[re(r, j)];
      }
    }
    for (int j = 0; j < N; ++j) {
      double inv = 1 / m_qr[re(i, i)];

      mat[re(i, j)] *= inv;
      mat[im(i, j)] *= inv;
    }
  }
}


//====================================================================
void Decompose_QR_Cmplx::get_Q(double *mat)
{
  for (int i = 0; i < size; ++i) {
    mat[i] = 0;
  }
  for (int i = 0; i < N; ++i) {
    mat[re(i, i)] = 1;
  }

  mult_Q(mat);
}


//====================================================================
void Decompose_QR_Cmplx::get_Qu(double *mat)
{
  get_Q(mat);
  std::valarray<double> sign(N);
  for (int i = 0; i < N; ++i) {
    sign[i] = (m_qr[re(i, i)] < 0) ? -1 : 1;
  }
  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < N; ++j) {
      mat[re(i, j)] *= sign[j];
      mat[im(i, j)] *= sign[j];
    }
  }
}


//====================================================================
void Decompose_QR_Cmplx::mult_Q(double *mat)
{
  for (int r = (int)N - 1; r >= 0; --r) {
    for (int j = 0; j < N; ++j) {
      // vm = <v|m_j>
      double vm_re = mat[re(r, j)];
      double vm_im = mat[im(r, j)];
      for (int i = r + 1; i < N; ++i) {
        vm_re += m_qr[re(i, r)] * mat[re(i, j)] + m_qr[im(i, r)] * mat[im(i, j)];
        vm_im += m_qr[re(i, r)] * mat[im(i, j)] - m_qr[im(i, r)] * mat[re(i, j)];
      }

      // beta = tau * va
      // m_j -= beta * v
      double beta_re = m_tau[re(r)] * vm_re - m_tau[im(r)] * vm_im;
      double beta_im = m_tau[re(r)] * vm_im + m_tau[im(r)] * vm_re;
      mat[re(r, j)] -= beta_re;
      mat[im(r, j)] -= beta_im;
      for (int i = r + 1; i < N; ++i) {
        double vi_re = m_qr[re(i, r)];
        double vi_im = m_qr[im(i, r)];
        mat[re(i, j)] -= beta_re * vi_re - beta_im * vi_im;
        mat[im(i, j)] -= beta_re * vi_im + beta_im * vi_re;
      }
    }
  }
}


//====================================================================
