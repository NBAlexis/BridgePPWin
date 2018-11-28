// This implementation only applies to SU(3) group and Nd=4 case.
#define NC      3
#define NC2     6
#define NDF     18
#define ND      4
#define NCD     12
#define NCD2    24

//====================================================================
namespace {
  void check_Nc()
  {
    vout.paranoiac(CommonParameters::Vlevel(),
                   "Field_F: implementation for SU(3).\n");
  }


  double mult_Gn_r(const double *g, const double *w, int)
  {
    return g[0] * w[0] - g[1] * w[1]
           + g[2] * w[2] - g[3] * w[3]
           + g[4] * w[4] - g[5] * w[5];
  }


  double mult_Gn_i(const double *g, const double *w, int)
  {
    return g[0] * w[1] + g[1] * w[0]
           + g[2] * w[3] + g[3] * w[2]
           + g[4] * w[5] + g[5] * w[4];
  }


  double mult_Gd_r(const double *g, const double *w, int)
  {
    return g[0] * w[0] + g[1] * w[1]
           + g[6] * w[2] + g[7] * w[3]
           + g[12] * w[4] + g[13] * w[5];
  }


  double mult_Gd_i(const double *g, const double *w, int)
  {
    return g[0] * w[1] - g[1] * w[0]
           + g[6] * w[3] - g[7] * w[2]
           + g[12] * w[5] - g[13] * w[4];
  }
} // end of nameless namespace
//====================================================================
//============================================================END=====
