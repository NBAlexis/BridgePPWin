#define  NC      m_Nc
#define  NCOL    m_Nc

//====================================================================
namespace {
  void check_Nc()
  {
    vout.paranoiac(CommonParameters::Vlevel(),
                   "Field_G: implementation for general SU(N).\n");
  }


  double mult_Gnn_r(const double *g1, const double *g2, int Nc)
  {
    double a = 0.0;

    for (int i = 0; i < Nc; ++i) {
      a += g1[2 * i] * g2[2 * i * Nc] - g1[2 * i + 1] * g2[2 * i * Nc + 1];
    }
    return a;
  }


  double mult_Gnn_i(const double *g1, const double *g2, int Nc)
  {
    double a = 0.0;

    for (int i = 0; i < Nc; ++i) {
      a += g1[2 * i] * g2[2 * i * Nc + 1] + g1[2 * i + 1] * g2[2 * i * Nc];
    }
    return a;
  }


  double mult_Gdn_r(const double *g1, const double *g2, int Nc)
  {
    double a = 0.0;

    for (int i = 0; i < Nc; ++i) {
      a += g1[2 * i * Nc] * g2[2 * i * Nc] + g1[2 * i * Nc + 1] * g2[2 * i * Nc + 1];
    }
    return a;
  }


  double mult_Gdn_i(const double *g1, const double *g2, int Nc)
  {
    double a = 0.0;

    for (int i = 0; i < Nc; ++i) {
      a += g1[2 * i * Nc] * g2[2 * i * Nc + 1] - g1[2 * i * Nc + 1] * g2[2 * i * Nc];
    }
    return a;
  }


  double mult_Gnd_r(const double *g1, const double *g2, int Nc)
  {
    double a = 0.0;

    for (int i = 0; i < Nc; ++i) {
      a += g1[2 * i] * g2[2 * i] + g1[2 * i + 1] * g2[2 * i + 1];
    }
    return a;
  }


  double mult_Gnd_i(const double *g1, const double *g2, int Nc)
  {
    double a = 0.0;

    for (int i = 0; i < Nc; ++i) {
      a += -g1[2 * i] * g2[2 * i + 1] + g1[2 * i + 1] * g2[2 * i];
    }
    return a;
  }


  double mult_Gdd_r(const double *g1, const double *g2, int Nc)
  {
    double a = 0.0;

    for (int i = 0; i < Nc; ++i) {
      a += g1[2 * i * Nc] * g2[2 * i] - g1[2 * i * Nc + 1] * g2[2 * i + 1];
    }
    return a;
  }


  double mult_Gdd_i(const double *g1, const double *g2, int Nc)
  {
    double a = 0.0;

    for (int i = 0; i < Nc; ++i) {
      a += -g1[2 * i * Nc] * g2[2 * i + 1] - g1[2 * i * Nc + 1] * g2[2 * i];
    }
    return a;
  }
} // end of nameless namespace
//====================================================================
//============================================================END=====
