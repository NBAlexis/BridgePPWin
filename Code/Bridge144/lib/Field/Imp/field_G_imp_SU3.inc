#define  NC      3
#define  NCOL    3

//====================================================================
namespace {
  void check_Nc()
  {
    vout.paranoiac(CommonParameters::Vlevel(),
                   "Field_G: implementation for SU(3).\n");
  }


  double mult_Gnn_r(const double *g1, const double *g2, int )
  {
    return g1[0] * g2[0] - g1[1] * g2[1]
           + g1[2] * g2[6] - g1[3] * g2[7]
           + g1[4] * g2[12] - g1[5] * g2[13];
  }


  double mult_Gnn_i(const double *g1, const double *g2, int )
  {
    return g1[0] * g2[1] + g1[1] * g2[0]
           + g1[2] * g2[7] + g1[3] * g2[6]
           + g1[4] * g2[13] + g1[5] * g2[12];
  }


  double mult_Gdn_r(const double *g1, const double *g2, int )
  {
    return g1[0] * g2[0] + g1[1] * g2[1]
           + g1[6] * g2[6] + g1[7] * g2[7]
           + g1[12] * g2[12] + g1[13] * g2[13];
  }


  double mult_Gdn_i(const double *g1, const double *g2, int )
  {
    return g1[0] * g2[1] - g1[1] * g2[0]
           + g1[6] * g2[7] - g1[7] * g2[6]
           + g1[12] * g2[13] - g1[13] * g2[12];
  }


  double mult_Gnd_r(const double *g1, const double *g2, int )
  {
    return g1[0] * g2[0] + g1[1] * g2[1]
           + g1[2] * g2[2] + g1[3] * g2[3]
           + g1[4] * g2[4] + g1[5] * g2[5];
  }


  double mult_Gnd_i(const double *g1, const double *g2, int )
  {
    return -g1[0] * g2[1] + g1[1] * g2[0]
           - g1[2] * g2[3] + g1[3] * g2[2]
           - g1[4] * g2[5] + g1[5] * g2[4];
  }


  double mult_Gdd_r(const double *g1, const double *g2, int )
  {
    return g1[0] * g2[0] - g1[1] * g2[1]
           + g1[6] * g2[2] - g1[7] * g2[3]
           + g1[12] * g2[4] - g1[13] * g2[5];
  }


  double mult_Gdd_i(const double *g1, const double *g2, int )
  {
    return -g1[0] * g2[1] - g1[1] * g2[0]
           - g1[6] * g2[3] - g1[7] * g2[2]
           - g1[12] * g2[5] - g1[13] * g2[4];
  }
} // end of nameless namespace
//====================================================================
//============================================================END=====
