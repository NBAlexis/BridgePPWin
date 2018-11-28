#include "BridgeLib_Private.h"
#if USE_ORG

/*!
        @file    $Id:: field_F.cpp #$

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 00:49:46 #$

        @version $LastChangedRevision: 1561 $
*/

#include "Field/field_F.h"
#include "Field/index_lex.h"

using std::valarray;

//====================================================================
void Field_F::check()
{
  // do nothing.
}


//====================================================================
void mult_Field_Gn(Field_F& y, const int ex,
                   const Field_G& U, int ex1,
                   const Field_F& x, int ex2)
{
  assert(ex < y.nex());
  assert(ex1 < U.nex());
  assert(ex2 < x.nex());
  assert(U.nvol() == y.nvol());
  assert(x.nvol() == y.nvol());

  const int Nvol = x.nvol();
  const int Nd   = x.nd();
  const int Nc   = x.nc();

  for (int site = 0; site < Nvol; ++site) {
    for (int s = 0; s < Nd; ++s) {
      y.set_vec(s, site, ex, U.mat(site, ex1) * x.vec(s, site, ex2));
    }
  }
}


//====================================================================
void mult_Field_Gd(Field_F& y, const int ex,
                   const Field_G& U, int ex1,
                   const Field_F& x, int ex2)
{
  assert(ex < y.nex());
  assert(ex1 < U.nex());
  assert(ex2 < x.nex());
  assert(U.nvol() == y.nvol());
  assert(x.nvol() == y.nvol());

  const int Nvol = y.nvol();
  const int Nd   = y.nd();
  const int Nc   = y.nc();

  for (int site = 0; site < Nvol; ++site) {
    for (int s = 0; s < Nd; ++s) {
      y.set_vec(s, site, ex, U.mat_dag(site, ex1) * x.vec(s, site, ex2));
    }
  }
}


//====================================================================
void multadd_Field_Gn(Field_F& y, const int ex,
                      const Field_G& U, int ex1,
                      const Field_F& x, int ex2,
                      const double a)
{
  assert(ex < y.nex());
  assert(ex1 < U.nex());
  assert(ex2 < x.nex());
  assert(U.nvol() == y.nvol());
  assert(x.nvol() == y.nvol());

  const int Nvol = y.nvol();
  const int Nd   = y.nd();
  const int Nc   = y.nc();

  Vec_SU_N vec(Nc);

  for (int site = 0; site < Nvol; ++site) {
    for (int s = 0; s < Nd; ++s) {
      y.add_vec(s, site, ex, U.mat(site, ex1) * x.vec(s, site, ex2) * a);
    }
  }
}


//====================================================================
void multadd_Field_Gd(Field_F& y, const int ex,
                      const Field_G& U, int ex1,
                      const Field_F& x, int ex2,
                      const double a)
{
  assert(ex < y.nex());
  assert(ex1 < U.nex());
  assert(ex2 < x.nex());
  assert(U.nvol() == y.nvol());
  assert(x.nvol() == y.nvol());

  const int Nvol = y.nvol();
  const int Nd   = y.nd();
  const int Nc   = y.nc();

  for (int site = 0; site < Nvol; ++site) {
    for (int s = 0; s < Nd; ++s) {
      y.add_vec(s, site, ex, U.mat_dag(site, ex1) * x.vec(s, site, ex2) * a);
    }
  }
}


//====================================================================
void mult_GM(Field_F& v, const GammaMatrix& gm, const Field_F& w)
{
  assert(w.nex() == v.nex());
  assert(w.nvol() == v.nvol());

  const int Nvol = v.nvol();
  const int Nex  = v.nex();
  const int Nd   = v.nd();
  const int Nc   = v.nc();

  valarray<int>    id(Nd);
  valarray<int>    idc_r(Nd);
  valarray<int>    idc_i(Nd);
  valarray<double> gv_r(Nd);
  valarray<double> gv_i(Nd);

  for (int s = 0; s < Nd; ++s) {
    id[s]    = gm.index(s);
    gv_r[s]  = gm.value_r(s);
    gv_i[s]  = gm.value_i(s);
    idc_r[s] = gm.index_c(s);
    idc_i[s] = 1 - idc_r[s];
  }

  for (int ex = 0; ex < Nex; ++ex) {
    for (int site = 0; site < Nvol; ++site) {
      for (int s = 0; s < Nd; ++s) {
        for (int c = 0; c < Nc; ++c) {
          double ww[2];
          ww[0] = w.cmp_r(c, id[s], site, ex);
          ww[1] = w.cmp_i(c, id[s], site, ex);

          v.set_ri(c, s, site, ex,
                   gv_r[s] * ww[idc_r[s]],
                   gv_i[s] * ww[idc_i[s]]);
        }
      }
    }
  }
}


//====================================================================
void mult_iGM(Field_F& v, const GammaMatrix& gm, const Field_F& w)
{
  assert(w.nex() == v.nex());
  assert(w.nvol() == v.nvol());

  const int Nvol = v.nvol();
  const int Nex  = v.nex();
  const int Nd   = v.nd();
  const int Nc   = v.nc();

  valarray<int>    id(Nd);
  valarray<int>    idc_r(Nd);
  valarray<int>    idc_i(Nd);
  valarray<double> gv_r(Nd);
  valarray<double> gv_i(Nd);

  for (int s = 0; s < Nd; ++s) {
    id[s]    = gm.index(s);
    gv_r[s]  = -gm.value_i(s);
    gv_i[s]  = gm.value_r(s);
    idc_r[s] = 1 - gm.index_c(s);
    idc_i[s] = 1 - idc_r[s];
  }

  for (int ex = 0; ex < Nex; ++ex) {
    for (int site = 0; site < Nvol; ++site) {
      for (int s = 0; s < Nd; ++s) {
        for (int c = 0; c < Nc; ++c) {
          double ww[2];
          ww[0] = w.cmp_r(c, id[s], site, ex);
          ww[1] = w.cmp_i(c, id[s], site, ex);

          v.set_ri(c, s, site, ex,
                   gv_r[s] * ww[idc_r[s]],
                   gv_i[s] * ww[idc_i[s]]);
        }
      }
    }
  }
}


//====================================================================
void mult_GMproj(Field_F& v,
                 const int pm, const GammaMatrix& gm,
                 const Field_F& w)
{
  assert(w.nvol() == v.nvol());
  assert(w.nex() == v.nex());
  assert(pm == 1 || pm == -1);

  mult_GM(v, gm, w);

  if (pm == 1) {
    axpy(v, 1.0, w);  // v += w;
    scal(v, 0.5);     // v *= 0.5;
  } else if (pm == -1) {
    axpy(v, -1.0, w); // v -= w     i.e. v = (gamma - 1) * w
    scal(v, -0.5);    // v *= -0.5  i.e. v = (1 - gamma)/2 * w
  } else {
    vout.crucial("Error at %s: wrong pm = %d\n", __func__, pm);
    exit(EXIT_FAILURE);
  }
}


//====================================================================
void mult_GMproj2(Field_F& v,
                  const int pm, const GammaMatrix& gm,
                  const Field_F& w)
{
  assert(w.nvol() == v.nvol());
  assert(w.nex() == v.nex());
  assert(pm == 1 || pm == -1);

  mult_GM(v, gm, w);

  if (pm == 1) {
    axpy(v, 1.0, w);   // v += w;
  } else if (pm == -1) {
    axpy(v, -1.0, w);  // v -= w;
    scal(v, -1.0);     // v *= -1.0;
  } else {
    vout.crucial("Error at %s: wrong pm = %d\n", __func__, pm);
    exit(EXIT_FAILURE);
  }
}


//====================================================================
void mult_GMproj2(Field_F& v,
                  const double nu_s,
                  const int pm,
                  const double r_s,
                  const GammaMatrix& gm,
                  const Field_F& w)
{
  assert(w.nvol() == v.nvol());
  assert(w.nex() == v.nex());
  assert(pm == 1 || pm == -1);

  mult_GM(v, gm, w);
  scal(v, nu_s);  // v *= nu_s;

  if (pm == 1) {
    axpy(v, r_s, w);   // v += r_s * w;
  } else if (pm == -1) {
    axpy(v, -r_s, w);  // v -= r_s * w  i.e. v = (nu_s * gamma - r_s) * w
    scal(v, -1.0);     // v *= -1.0     i.e. v = (r_s - nu_s * gamma) * w
  } else {
    vout.crucial("Error at %s: wrong pm = %d\n", __func__, pm);
    exit(EXIT_FAILURE);
  }
}


//====================================================================
//============================================================END=====
#endif 
