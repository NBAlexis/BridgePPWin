#include "BridgeLib_Private.h"

/*!
        @file    $Id:: eigen_QR_Cmplx.cpp #$

        @brief

        @author  Satoru Ueda  (sueda)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 00:49:46 #$

        @version $LastChangedRevision: 1561 $
*/

#include "eigen_QR_Cmplx.h"

const std::string Eigen_QR_Cmplx::class_name = "Eigen_QR_Cmplx";

//====================================================================
std::valarray<double> Eigen_QR_Cmplx::solve(const double *matrix)
{
  int    Niter     = 1000;
  double Stop_cond = DBL_EPSILON;

  std::valarray<double> eigen_val(N2);

  for (int i = 0; i < size; ++i) {
    m_mat[i] = matrix[i];
    m_q[i]   = 0;
  }
  for (int i = 0; i < N; ++i) {
    m_q[i * N2 + i * 2] = 1;
  }

  // change to hessenberg matrix
  Decompose_Hessenberg_Cmplx hessenberg(N);
  hessenberg.set_matrix(&m_mat[0]);
  hessenberg.get_Hessenberg(&m_mat[0]);
  hessenberg.get_Q(&m_q[0]);

  // # of unfound eigen value
  int Neigen = N;
  for (int i_iter = 0; i_iter < Niter; ++i_iter) {
    int r1 = Neigen - 1;

    // ****
    // ****
    //  *ab
    //   cd
    double d_re = m_mat[re(r1, r1)];
    double d_im = m_mat[im(r1, r1)];

    // get the last eigen value
    if (Neigen == 1) {
      eigen_val[re(0)] = d_re;
      eigen_val[im(0)] = d_im;
      break;
    }

    int     r2   = Neigen - 2;
    double  a_re = m_mat[re(r2, r2)];
    double  a_im = m_mat[im(r2, r2)];
    double& c_re = m_mat[re(r1, r2)];
    double& c_im = m_mat[im(r1, r2)];

    double norm_ad = a_re * a_re + a_im * a_im + d_re * d_re + d_im * d_im;
    double norm_c  = c_re * c_re + c_im * c_im;

    if (norm_c < norm_ad * Stop_cond * Stop_cond) {
      // When c is very small, d become eigen value.
      eigen_val[re(r1)] = d_re;
      eigen_val[im(r1)] = d_im;
      c_re    = c_im = 0;
      d_re    = a_re;
      d_im    = a_im;
      Neigen -= 1;
    }

    // shift
    for (int r = 0; r < Neigen; ++r) {
      m_mat[re(r, r)] -= d_re;
      m_mat[im(r, r)] -= d_im;
    }

    // qr step
    qr_step(Neigen);

    // shift back
    for (int r = 0; r < Neigen; ++r) {
      m_mat[re(r, r)] += d_re;
      m_mat[im(r, r)] += d_im;
    }
  }

  return eigen_val;
}


//====================================================================
void Eigen_QR_Cmplx::qr_step(const int rank)
{
  std::valarray<double> tau(rank * 2);

  for (int r = 0; r < rank - 1; ++r) {
    // diagnal part
    double& alpha_re = m_mat[re(r, r)];
    double& alpha_im = m_mat[im(r, r)];
    // sub diagnal part
    double& v1_re = m_mat[re(r + 1, r)];
    double& v1_im = m_mat[im(r + 1, r)];

    double nu = alpha_re * alpha_re + alpha_im * alpha_im + v1_re * v1_re + v1_im * v1_im;

    nu         = std::sqrt(nu) * ((alpha_re < 0) ? -1 : 1);
    tau[re(r)] = alpha_re + nu;
    tau[im(r)] = alpha_im;

    double norm2_sigma  = tau[re(r)] * tau[re(r)] + tau[im(r)] * tau[im(r)];
    double inv_sigma_re = tau[re(r)] / norm2_sigma;
    double inv_sigma_im = -tau[im(r)] / norm2_sigma;

    // householder factor
    tau[re(r)] /= nu;
    tau[im(r)] /= nu;

    alpha_re = -nu;
    alpha_im = 0;

    // householder vector
    double v_re = v1_re;
    double v_im = v1_im;
    v1_re = v_re * inv_sigma_re - v_im * inv_sigma_im;
    v1_im = v_im * inv_sigma_re + v_re * inv_sigma_im;

    // v = (1, v_1)^T
    for (int j = r + 1; j < N; ++j) {
      // va = <v|a_j>
      double va_re = m_mat[re(r, j)];
      double va_im = m_mat[im(r, j)];
      va_re += v1_re * m_mat[re(r + 1, j)] + v1_im * m_mat[im(r + 1, j)];
      va_im += v1_re * m_mat[im(r + 1, j)] - v1_im * m_mat[re(r + 1, j)];

      double beta_re = tau[re(r)] * va_re + tau[im(r)] * va_im;
      double beta_im = tau[re(r)] * va_im - tau[im(r)] * va_re;
      // i = r
      m_mat[re(r, j)] -= beta_re;
      m_mat[im(r, j)] -= beta_im;
      // i = r+1
      m_mat[re(r + 1, j)] -= beta_re * v1_re - beta_im * v1_im;
      m_mat[im(r + 1, j)] -= beta_re * v1_im + beta_im * v1_re;
    }

    // Q -> Q Q_r
    for (int i = 0; i < N; ++i) {
      // va = <v|Q_i^*>
      double va_re = m_q[re(i, r)];
      double va_im = -m_q[im(i, r)];
      va_re += v1_re * m_q[re(i, r + 1)] - v1_im * m_q[im(i, r + 1)];
      va_im -= v1_re * m_q[im(i, r + 1)] + v1_im * m_q[re(i, r + 1)];

      double beta_re = tau[re(r)] * va_re + tau[im(r)] * va_im;
      double beta_im = tau[re(r)] * va_im - tau[im(r)] * va_re;
      // j = r
      m_q[re(i, r)] -= beta_re;
      m_q[im(i, r)] += beta_im;
      // j = r+1
      m_q[re(i, r + 1)] -= beta_re * v1_re - beta_im * v1_im;
      m_q[im(i, r + 1)] += beta_re * v1_im + beta_im * v1_re;
    }
  }

  //A -> RQ
  for (int r = 0; r < rank - 1; ++r) {
    double v1_re = m_mat[re(r + 1, r)];
    double v1_im = m_mat[im(r + 1, r)];
    m_mat[re(r + 1, r)] = m_mat[im(r + 1, r)] = 0;
    for (int i = 0; i <= r + 1; ++i) {
      double va_re = m_mat[re(i, r)];
      double va_im = -m_mat[im(i, r)];
      va_re += v1_re * m_mat[re(i, r + 1)] - v1_im * m_mat[im(i, r + 1)];
      va_im -= v1_re * m_mat[im(i, r + 1)] + v1_im * m_mat[re(i, r + 1)];

      // beta = tau^* * va
      // a_i^* -= bata * v
      double beta_re = tau[re(r)] * va_re + tau[im(r)] * va_im;
      double beta_im = tau[re(r)] * va_im - tau[im(r)] * va_re;
      // j = r
      m_mat[re(i, r)] -= beta_re;
      m_mat[im(i, r)] += beta_im;
      // j = r+1
      m_mat[re(i, r + 1)] -= beta_re * v1_re - beta_im * v1_im;
      m_mat[im(i, r + 1)] += beta_re * v1_im + beta_im * v1_re;
    }
  }
}


//====================================================================
void Eigen_QR_Cmplx::get_Q(double *q)
{
  for (int i = 0; i < size; ++i) {
    q[i] = m_q[i];
  }
}


//====================================================================
void Eigen_QR_Cmplx::get_R(double *r)
{
  for (int i = 0; i < size; ++i) {
    r[i] = m_mat[i];
  }

  /*
  for (int i = 0; i < N; ++i) {
    for (int j = i; j < N; ++j) {
      r[re(i,j)] = m_mat[re(i,j)];
      r[im(i,j)] = m_mat[im(i,j)];
    }
  }
  */
}


//====================================================================
