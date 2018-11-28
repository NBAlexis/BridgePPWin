#define  NC      3
#define  NCOL    3

//====================================================================
namespace {
  void check_Nc()
  {
    vout.general(CommonParameters::Vlevel(),
                 "Fopr_Wilson_impl: implementation for SU(3).\n");
  }


  inline double mult_uv_r(const double *g, const double *w, int Nc)
  {
    return g[0] * w[0] - g[1] * w[1]
           + g[2] * w[2] - g[3] * w[3]
           + g[4] * w[4] - g[5] * w[5];
  }


  inline double mult_uv_i(const double *g, const double *w, int Nc)
  {
    return g[0] * w[1] + g[1] * w[0]
           + g[2] * w[3] + g[3] * w[2]
           + g[4] * w[5] + g[5] * w[4];
  }


  inline double mult_udagv_r(const double *g, const double *w, int Nc)
  {
    return g[0] * w[0] + g[1] * w[1]
           + g[6] * w[2] + g[7] * w[3]
           + g[12] * w[4] + g[13] * w[5];
  }


  inline double mult_udagv_i(const double *g, const double *w, int Nc)
  {
    return g[0] * w[1] - g[1] * w[0]
           + g[6] * w[3] - g[7] * w[2]
           + g[12] * w[5] - g[13] * w[4];
  }


  inline dcomplex mult_uv(const dcomplex *u, const dcomplex *v, int Nc)
  {
    return u[0] * v[0] + u[1] * v[1] + u[2] * v[2];
  }


  inline dcomplex mult_udagv(const dcomplex *u, const dcomplex *v, int Nc)
  {
    return conj(u[0]) * v[0] + conj(u[3]) * v[1] + conj(u[6]) * v[2];
  }
} // end of nameless namespace
//====================================================================
//============================================================END=====
