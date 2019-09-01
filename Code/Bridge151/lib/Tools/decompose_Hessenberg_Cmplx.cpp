/*!
        @file    decompose_Hessenberg_Cmplx.cpp

        @brief

        @author  Satoru Ueda  (sueda)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/
#include "BridgeLib_Private.h"
#include "decompose_Hessenberg_Cmplx.h"

//====================================================================
void Decompose_Hessenberg_Cmplx::set_matrix(const double *mat)
{
  // double epsilon = DBL_EPSILON;

  // matrix copy
  for (int i = 0; i < size; ++i) {
    m_H[i] = mat[i];
  }

  // Householder trans for r col.
  for (int r = 0; r < N - 1; ++r) {
    // sub diagnal part
    double& alpha_re = m_H[re(r + 1, r)];
    double& alpha_im = m_H[im(r + 1, r)];

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
    for (int i = r + 1; i < N; ++i) {
      nu += m_H[re(i, r)] * m_H[re(i, r)] + m_H[im(i, r)] * m_H[im(i, r)];
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

    // m_H[r+1,r] = -nu
    // m_H[i,r] /= sigma
    alpha_re = -nu;
    alpha_im = 0;
    for (int i = r + 2; i < N; ++i) {
      double ir_re = m_H[re(i, r)];
      double ir_im = m_H[im(i, r)];

      m_H[re(i, r)] = ir_re * inv_sigma_re - ir_im * inv_sigma_im;
      m_H[im(i, r)] = ir_im * inv_sigma_re + ir_re * inv_sigma_im;
    }

    // H^\dag a_j
    for (int j = r + 1; j < N; ++j) {
      // va = <v|a_j>
      double va_re = m_H[re(r + 1, j)];
      double va_im = m_H[im(r + 1, j)];
      for (int i = r + 2; i < N; ++i) {
        va_re += m_H[re(i, r)] * m_H[re(i, j)] + m_H[im(i, r)] * m_H[im(i, j)];
        va_im += m_H[re(i, r)] * m_H[im(i, j)] - m_H[im(i, r)] * m_H[re(i, j)];
      }
      // beta = tau^* * va
      // a_j -= bata * v
      double beta_re = m_tau[re(r)] * va_re + m_tau[im(r)] * va_im;
      double beta_im = m_tau[re(r)] * va_im - m_tau[im(r)] * va_re;
      m_H[re(r + 1, j)] -= beta_re;
      m_H[im(r + 1, j)] -= beta_im;
      for (int i = r + 2; i < N; ++i) {
        double vi_re = m_H[re(i, r)];
        double vi_im = m_H[im(i, r)];

        m_H[re(i, j)] -= beta_re * vi_re - beta_im * vi_im;
        m_H[im(i, j)] -= beta_re * vi_im + beta_im * vi_re;
      }
    }
    // AH = (H^\dag A^\dag)^\dag
    // A = (a_0, ...,a_{n-1})^T
    for (int i = 0; i < N; ++i) {
      // va = <v|a_i^*> v=H.col(r)
      double va_re = m_H[re(i, r + 1)];
      double va_im = -m_H[im(i, r + 1)];
      for (int j = r + 2; j < N; ++j) {
        va_re += m_H[re(j, r)] * m_H[re(i, j)] - m_H[im(j, r)] * m_H[im(i, j)];
        va_im -= m_H[re(j, r)] * m_H[im(i, j)] + m_H[im(j, r)] * m_H[re(i, j)];
      }
      // beta = tau^* * va
      // a_i^* -= bata * v
      double beta_re = m_tau[re(r)] * va_re + m_tau[im(r)] * va_im;
      double beta_im = m_tau[re(r)] * va_im - m_tau[im(r)] * va_re;
      m_H[re(i, r + 1)] -= beta_re;
      m_H[im(i, r + 1)] += beta_im;
      for (int j = r + 2; j < N; ++j) {
        double vj_re = m_H[re(j, r)];
        double vj_im = m_H[im(j, r)];

        m_H[re(i, j)] -= beta_re * vj_re - beta_im * vj_im;
        m_H[im(i, j)] += beta_re * vj_im + beta_im * vj_re;
      }
    }
  }
}


//====================================================================
void Decompose_Hessenberg_Cmplx::get_Hessenberg(double *mat)
{
  for (int i = 0; i < size; ++i) {
    mat[i] = m_H[i];
  }
  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < i - 1; ++j) {
      mat[re(i, j)] = mat[im(i, j)] = 0;
    }
  }
}


//====================================================================
void Decompose_Hessenberg_Cmplx::get_Q(double *mat)
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
void Decompose_Hessenberg_Cmplx::mult_Q(double *mat)
{
  for (int r = static_cast<int>(N) - 2; r >= 0; --r) {
    for (int j = 0; j < static_cast<int>(N); ++j) {
      // vm = <v|m_j>
      double vm_re = mat[re(r + 1, j)];
      double vm_im = mat[im(r + 1, j)];
      for (int i = r + 2; i < N; ++i) {
        vm_re += m_H[re(i, r)] * mat[re(i, j)] + m_H[im(i, r)] * mat[im(i, j)];
        vm_im += m_H[re(i, r)] * mat[im(i, j)] - m_H[im(i, r)] * mat[re(i, j)];
      }
      // beta = tau * va
      // m_j -= beta * v
      double beta_re = m_tau[re(r)] * vm_re - m_tau[im(r)] * vm_im;
      double beta_im = m_tau[re(r)] * vm_im + m_tau[im(r)] * vm_re;
      mat[re(r + 1, j)] -= beta_re;
      mat[im(r + 1, j)] -= beta_im;
      for (int i = r + 2; i < N; ++i) {
        double vi_re = m_H[re(i, r)];
        double vi_im = m_H[im(i, r)];
        mat[re(i, j)] -= beta_re * vi_re - beta_im * vi_im;
        mat[im(i, j)] -= beta_re * vi_im + beta_im * vi_re;
      }
    }
  }
}


//====================================================================
