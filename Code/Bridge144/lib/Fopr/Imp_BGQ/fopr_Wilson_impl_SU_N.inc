#define  NC      m_Nc
#define  NCOL    m_Nc

//====================================================================
namespace {
  void check_Nc()
  {
    vout.general(CommonParameters::Vlevel(),
                 "Fopr_Wilson_impl: implementation for general SU(N).\n");
  }


  inline double mult_uv_r(const double *g, const double *w, int Nc)
  {
    double a = 0.0;

    for (int i = 0; i < Nc; ++i) {
      a += g[2 * i] * w[2 * i] - g[2 * i + 1] * w[2 * i + 1];
    }
    return a;
  }


  inline double mult_uv_i(const double *g, const double *w, int Nc)
  {
    double a = 0.0;

    for (int i = 0; i < Nc; ++i) {
      a += g[2 * i] * w[2 * i + 1] + g[2 * i + 1] * w[2 * i];
    }
    return a;
  }


  inline double mult_udagv_r(const double *g, const double *w, int Nc)
  {
    double a = 0.0;

    for (int i = 0; i < Nc; ++i) {
      a += g[2 * i * Nc] * w[2 * i] + g[2 * i * Nc + 1] * w[2 * i + 1];
    }
    return a;
  }


  inline double mult_udagv_i(const double *g, const double *w, int Nc)
  {
    double a = 0.0;

    for (int i = 0; i < Nc; ++i) {
      a += g[2 * i * Nc] * w[2 * i + 1] - g[2 * i * Nc + 1] * w[2 * i];
    }
    return a;
  }
} // end of nameless namespace
//====================================================================
//============================================================END=====
