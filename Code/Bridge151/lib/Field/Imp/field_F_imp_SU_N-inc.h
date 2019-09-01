// This implementation only applies to SU(3) group and Nd=4 case.
#define NC      m_Nc
#define NC2     (2 * m_Nc)
#define NDF     (2 * m_Nc * m_Nc)
#define ND      4
#define NCD     (4 * m_Nc)
#define NCD2    (8 * m_Nc)

//====================================================================
namespace {
  void check_Nc()
  {
    vout.paranoiac(CommonParameters::Vlevel(),
                   "Field_F: implementation for general SU(N).\n");
  }


  double mult_Gn_r(const double *g, const double *w, int Nc)
  {
    double a = 0.0;

    for (int i = 0; i < Nc; ++i) {
      a += g[2 * i] * w[2 * i] - g[2 * i + 1] * w[2 * i + 1];
    }
    return a;
  }


  double mult_Gn_i(const double *g, const double *w, int Nc)
  {
    double a = 0.0;

    for (int i = 0; i < Nc; ++i) {
      a += g[2 * i] * w[2 * i + 1] + g[2 * i + 1] * w[2 * i];
    }
    return a;
  }


  double mult_Gd_r(const double *g, const double *w, int Nc)
  {
    double a = 0.0;

    for (int i = 0; i < Nc; ++i) {
      a += g[2 * i * Nc] * w[2 * i] + g[2 * i * Nc + 1] * w[2 * i + 1];
    }
    return a;
  }


  double mult_Gd_i(const double *g, const double *w, int Nc)
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
