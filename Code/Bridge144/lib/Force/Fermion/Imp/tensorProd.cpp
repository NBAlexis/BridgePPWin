#include "BridgeLib_Private.h"
#if USE_IMP

/*!
        @file    $Id: tensorProd.cpp #$
        @brief
        @author
                 $LastChangedBy: aoym $
        @date    $LastChangedDate: 2013-03-21 15:28:34 #$
        @version $LastChangedRevision: 1561 $
*/

#include "Force/Fermion/tensorProd.h"
#include <cassert>

// This implementation only applies to SU(3) group and Nd=4 case.
#if defined USE_GROUP_SU3
#define NC      3
#define NC2     6
#define NDF     18
#define ND      4
#define NCD2    24
#define C1      0
#define C2      2
#define C3      4
#elif defined USE_GROUP_SU2
#define NC      2
#define NC2     4
#define NDF     8
#define ND      4
#define NCD2    16
#endif

//====================================================================
void tensorProd_Field_F(Field_G& u, const Field_F& v1, const Field_F& v2)
{
  return tensorProd_Field_F(u, 0, v1, v2);
}


//====================================================================
void tensorProd_Field_F(Field_G& u, const int ex, const Field_F& v1, const Field_F& v2)
{
  int Nvol = u.nvol();

  assert(Nvol == v1.nvol());
  assert(Nvol == v2.nvol());
  assert(ex < u.nex());
  assert(v1.nex() == 1);
  assert(v2.nex() == 1);

#if defined USE_GROUP_SU_N
  int NC   = CommonParameters::Nc();
  int ND   = CommonParameters::Nd();
  int NC2  = 2 * NC;
  int NDF  = 2 * NC * NC;
  int NCD2 = NC2 * ND;
#endif

  const double *w1 = v1.ptr(0);
  const double *w2 = v2.ptr(0);
  double       *g  = u.ptr(0, 0, ex);

  for (int site = 0; site < Nvol; ++site) {
    int iw = NCD2 * site;
    int ig = NDF * site;
    for (int c1 = 0; c1 < NC; ++c1) {
      for (int c2 = 0; c2 < NC; ++c2) {
        int ig2 = c2 * 2 + c1 * NC2 + ig;
        g[ig2]     = 0.0;
        g[ig2 + 1] = 0.0;
        for (int s = 0; s < ND; ++s) {
          g[ig2] +=
            w1[2 * c2 + s * NC2 + iw] * w2[2 * c1 + s * NC2 + iw]
            + w1[2 * c2 + 1 + s * NC2 + iw] * w2[2 * c1 + 1 + s * NC2 + iw];
          g[ig2 + 1] +=
            w1[2 * c2 + s * NC2 + iw] * w2[2 * c1 + 1 + s * NC2 + iw]
            - w1[2 * c2 + 1 + s * NC2 + iw] * w2[2 * c1 + s * NC2 + iw];
        }
      }
    }
  }
}


//================================================================
#endif
