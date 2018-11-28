#define  NC      2
#define  NCOL    2

//====================================================================
namespace {
  void check_Nc()
  {
    vout.general(CommonParameters::Vlevel(),
                 "Fopr_Wilson_impl: implementation for SU(2).\n");
  }


  inline double mult_uv_r(const double *u, const double *v, int Nc)
  {
    return u[0] * v[0] - u[1] * v[1]
           + u[2] * v[2] - u[3] * v[3];
  }


  inline double mult_uv_i(const double *u, const double *v, int Nc)
  {
    return u[0] * v[1] + u[1] * v[0]
           + u[2] * v[3] + u[3] * v[2];
  }


  inline double mult_udagv_r(const double *u, const double *v, int Nc)
  {
    return u[0] * v[0] + u[1] * v[1]
           + u[4] * v[2] + u[5] * v[3];
  }


  inline double mult_udagv_i(const double *u, const double *v, int Nc)
  {
    return u[0] * v[1] - u[1] * v[0]
           + u[4] * v[3] - u[5] * v[2];
  }


  inline dcomplex mult_uv(const dcomplex *u, const dcomplex *v, int Nc)
  {
    return u[0] * v[0] + u[1] * v[1];
  }


  inline dcomplex mult_udagv(const dcomplex *u, const dcomplex *v, int Nc)
  {
    return conj(u[0]) * v[0] + conj(u[2]) * v[1];
  }
} // end of nameless namespace
//====================================================================
//============================================================END=====
