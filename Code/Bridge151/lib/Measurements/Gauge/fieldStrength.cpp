/*!
        @file    fieldStrength.cpp

        @brief

        @author  Yusuke Namekawa  (namekawa)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/
#include "BridgeLib_Private.h"
#include "fieldStrength.h"

const std::string FieldStrength::class_name = "FieldStrengh";

//====================================================================
void FieldStrength::construct_Fmunu_1x1(Field_G& Fmunu_1x1,
                                        const int mu, const int nu, const Field_G& U)
{
  const int Nvol = CommonParameters::Nvol();

  //- building blocks
  // (1)  mu (2)
  //    +->-+
  // nu |   |
  //   i+   +
  Field_G Cup;

  m_staple.upper(Cup, U, mu, nu);

  // (1)  mu (2)
  //   i+   +
  // nu |   |
  //    +->-+
  Field_G Cdn;
  m_staple.lower(Cdn, U, mu, nu);

  Field_G Umu;
  Umu.setpart_ex(0, U, mu);

  //---------------------------------
  //- right clover leaf
  // (2)  +-<-+
  //  nu  |   |
  //     i+->-+
  //      mu (1)
  mult_Field_Gnd(Fmunu_1x1, 0, Umu, 0, Cup, 0);

  // (1) i+-<-+
  //  nu  |   |
  //      +->-+
  //      mu (2)
  multadd_Field_Gnd(Fmunu_1x1, 0, Cdn, 0, Umu, 0, 1.0);
  //---------------------------------

  //---------------------------------
  //- left clover leaf
  //   mu (2)
  //   +-<-+ (1)
  //   |   | nu    mu (1)
  //   +->-+i   +  +-<-+i (2)
  //               |   | nu
  //               +->-+
  Field_G v;
  mult_Field_Gdn(v, 0, Cup, 0, Umu, 0);
  multadd_Field_Gdn(v, 0, Umu, 0, Cdn, 0, 1.0);

  // NB. shift.forward(mu)
  //   mu (2)        mu (2)
  //   +-<-+ (1)     +-<-+ (1)
  //   |   | nu  ->  |   | nu
  //  i+->-+         +->-+i
  Field_G v2;
  m_shift.forward(v2, v, mu);

  axpy(Fmunu_1x1, 1.0, v2);
  //---------------------------------

  ah_Field_G(Fmunu_1x1, 0);

  scal(Fmunu_1x1, -0.25);
  Fmunu_1x1.xI();
}


//====================================================================
void FieldStrength::construct_Fmunu_1x2(Field_G& Fmunu_1x2,
                                        const int mu, const int nu, const Field_G& U)
{
  const int Nvol = CommonParameters::Nvol();


  //-- building blocks
  //      mu (2)
  // (1)  +->-+
  //  nu  |   |
  //     i+   +
  Field_G Cup1;

  m_staple.upper(Cup1, U, mu, nu);

  //      +-<-+ (2)
  //          | nu
  //     i+->-+
  //      mu (1)
  Field_G Cup2;
  m_staple.upper(Cup2, U, nu, mu);

  // (1) i+   +
  //  nu  |   |
  //      +->-+
  //      mu (2)
  Field_G Cdn1;
  m_staple.lower(Cdn1, U, mu, nu);

  // (2)  +->-+
  //  nu  |
  //      +-<-+i
  //      mu (1)
  Field_G Cdn2;
  m_staple.lower(Cdn2, U, nu, mu);

  Field_G Umu;
  Umu.setpart_ex(0, U, mu);

  Field_G Unu;
  Unu.setpart_ex(0, U, nu);

  //  +->-+
  //
  // i+
  Field_G Umu_nu;
  m_shift.backward(Umu_nu, Umu, nu);

  //      +
  //      |
  // i+   +
  Field_G Unu_mu;
  m_shift.backward(Unu_mu, Unu, mu);

  //---------------------------------
  //-- 2x1 part

  //- upper-right(2x1)
  //    +---<---+      +      +-<-+         <---+
  // nu |       |  =   |  +          +          |  +
  //   i+--->---+     i+     i+         i+  >---+     i+->-+
  Field_G c;
  m_shift.backward(c, Cup2, mu);

  Field_G v;
  mult_Field_Gdd(v, 0, Umu_nu, 0, Unu, 0);

  Field_G w;
  mult_Field_Gnn(w, 0, c, 0, v, 0);

  Field_G rect;
  mult_Field_Gnn(rect, 0, Umu, 0, w, 0);

  Fmunu_1x2 = rect;

  //- lower-right(2x1)
  // (1) +---<---+     (1) i+---<---+
  //  nu |       |  ->  nu  |       |
  //    i+--->---+          +--->---+
  //      mu (2)             mu (2)
  mult_Field_Gnd(v, 0, c, 0, Umu_nu, 0);
  mult_Field_Gnn(w, 0, Umu, 0, v, 0);

  mult_Field_Gdn(rect, 0, Unu, 0, w, 0);

  m_shift.forward(v, rect, nu);

  axpy(Fmunu_1x2, 1.0, v);

  //- upper-left(2x1)
  //      mu (2)
  //  +---<---+ (1)              +-<-+      +-<-+           +
  //  |       | nu  =         +  |       +         +        |
  //  +---i---+       i+->-+     +->-+i     +i         +i   +
  //
  // NB. shift.forward(mu)
  //      mu (2)             mu (2)
  //  +---<---+ (1)      +---<---+ (1)
  //  |       | nu  ->   |       | nu
  //  +---i---+          +--->---+i
  mult_Field_Gdn(v, 0, Cdn2, 0, Umu, 0);
  mult_Field_Gdn(w, 0, Umu_nu, 0, v, 0);

  mult_Field_Gnn(rect, 0, Unu_mu, 0, w, 0);

  m_shift.forward(v, rect, mu);

  axpy(Fmunu_1x2, 1.0, v);

  //- lower-left(2x1)
  //      mu (1)
  //  +---<---+ (2)        +                +-<-+      +-<-+
  //  |       | nu  =      |  +          +  |       +
  //  +---i---+       +i   +     i+->-+     +->-+i     +i
  //
  // NB. shift.forward(mu+nu)
  //      mu (1)             mu (1)            mu (1)
  //  +---<---+ (2)      +---<---+ (2)     +---<---+i (2)
  //  |       | nu  ->   |       | nu  ->  |       |  nu
  //  +---i---+          +--->---+i        +--->---+
  mult_Field_Gnn(v, 0, Umu, 0, Unu_mu, 0);
  mult_Field_Gdn(w, 0, Cdn2, 0, v, 0);

  mult_Field_Gdn(rect, 0, Umu_nu, 0, w, 0);

  m_shift.forward(w, rect, mu);
  m_shift.forward(v, w, nu);

  axpy(Fmunu_1x2, 1.0, v);
  //---------------------------------


  //---------------------------------
  //-- 1x2 part

  //- upper-right(1x2)
  //     +-<-+             +-<-+
  //     |   |             |   |
  // (2) +   +  =  +   +   +   +   +      +   +
  //  nu |   |     |                      |
  //    i+->-+    i+      i+         i+   +      i+->-+
  //     mu (1)
  m_shift.backward(c, Cup1, nu);

  mult_Field_Gdd(v, 0, c, 0, Unu, 0);
  mult_Field_Gnn(w, 0, Unu_mu, 0, v, 0);

  mult_Field_Gnn(rect, 0, Umu, 0, w, 0);

  axpy(Fmunu_1x2, 1.0, rect);

  //- lower-right(1x2)
  //      mu (2)
  // (1)  +-<-+      +-<-+           +                     +
  //  nu  |   |                      |                     |
  //     i+   +  =  i+      +   i+   +   +   i+   +   +   i+
  //      |   |                               |   |
  //      +->-+                               +->-+
  //
  // NB. shift.forward(nu)
  //      mu (2)        mu (2)
  // (1)  +-<-+    (1) i+-<-+
  //  nu  |   |     nu  |   |
  //     i+   +  ->     +   +
  //      |   |         |   |
  //      +->-+         +->-+
  mult_Field_Gnd(v, 0, Unu_mu, 0, Umu_nu, 0);
  mult_Field_Gnn(w, 0, Cdn1, 0, v, 0);

  mult_Field_Gdn(rect, 0, Unu, 0, w, 0);

  m_shift.forward(v, rect, nu);

  axpy(Fmunu_1x2, 1.0, v);

  //- upper-left(1x2)
  //  +-<-+                            +-<-+
  //  |   |                            |   |
  //  +   + (1)  =         +   +   +   +   +   +       +
  //  |   | nu                 |                       |
  // i+->-+        i+->-+     i+      i+          i+   +
  //  mu (2)
  //
  // NB. shift.forward(mu)
  //  +-<-+         +-<-+
  //  |   |         |   |
  //  +   + (1) ->  +   + (1)
  //  |   | nu      |   | nu
  // i+->-+         +->-+i
  //  mu (2)        mu (2)
  mult_Field_Gdn(v, 0, Unu, 0, Umu, 0);
  mult_Field_Gdn(w, 0, c, 0, v, 0);

  mult_Field_Gnn(rect, 0, Unu_mu, 0, w, 0);

  m_shift.forward(v, rect, mu);

  axpy(Fmunu_1x2, 1.0, v);

  //- lower-left(1x2)
  //      mu (1)
  // (2)  +-<-+          +                     +        +-<-+
  //  nu  |   |          |                     |
  //     i+   +  =  i+   +   +   i+   +   +   i+   +   i+
  //      |   |                   |   |
  //      +->-+                   +->-+
  //
  // NB. shift.forward(mu+nu)
  //      mu (1)        mu (1)
  // (2)  +-<-+    (2)  +-<-+i
  //  nu  |   |     nu  |   |
  //     i+   +  ->     +   +
  //      |   |         |   |
  //      +->-+         +->-+
  mult_Field_Gnn(v, 0, Cdn1, 0, Unu_mu, 0);
  mult_Field_Gdn(w, 0, Unu, 0, v, 0);

  mult_Field_Gdn(rect, 0, Umu_nu, 0, w, 0);

  m_shift.forward(w, rect, mu);
  m_shift.forward(v, w, nu);

  axpy(Fmunu_1x2, 1.0, v);
  //---------------------------------

  ah_Field_G(Fmunu_1x2, 0);

  //- normalization = 1/8
  //  i.e. overall factor = 1/4, average of 1x2 and 2x1 = 1/2
  scal(Fmunu_1x2, -0.125);
  Fmunu_1x2.xI();
}


//====================================================================
void FieldStrength::construct_Fmunu_1x1_traceless(Field_G& Fmunu_1x1,
                                                  const int mu, const int nu, const Field_G& U)
{
  const int Nvol = CommonParameters::Nvol();

  //- building blocks
  // (1)  mu (2)
  //    +->-+
  // nu |   |
  //   i+   +
  Field_G Cup;

  m_staple.upper(Cup, U, mu, nu);

  // (1)  mu (2)
  //   i+   +
  // nu |   |
  //    +->-+
  Field_G Cdn;
  m_staple.lower(Cdn, U, mu, nu);

  Field_G Umu;
  Umu.setpart_ex(0, U, mu);

  //---------------------------------
  //- right clover leaf
  // (2)  +-<-+
  //  nu  |   |
  //     i+->-+
  //      mu (1)
  mult_Field_Gnd(Fmunu_1x1, 0, Umu, 0, Cup, 0);

  // (1) i+-<-+
  //  nu  |   |
  //      +->-+
  //      mu (2)
  multadd_Field_Gnd(Fmunu_1x1, 0, Cdn, 0, Umu, 0, 1.0);
  //---------------------------------

  //---------------------------------
  //- left clover leaf
  //   mu (2)
  //   +-<-+ (1)
  //   |   | nu    mu (1)
  //   +->-+i   +  +-<-+i (2)
  //               |   | nu
  //               +->-+
  Field_G v;
  mult_Field_Gdn(v, 0, Cup, 0, Umu, 0);
  multadd_Field_Gdn(v, 0, Umu, 0, Cdn, 0, 1.0);

  // NB. shift.forward(mu)
  //   mu (2)        mu (2)
  //   +-<-+ (1)     +-<-+ (1)
  //   |   | nu  ->  |   | nu
  //  i+->-+         +->-+i
  Field_G v2;
  m_shift.forward(v2, v, mu);

  axpy(Fmunu_1x1, 1.0, v2);
  //---------------------------------

  // ah_Field_G(Fmunu_1x1, 0);
  at_Field_G(Fmunu_1x1, 0);

  scal(Fmunu_1x1, -0.25);
  Fmunu_1x1.xI();
}


//====================================================================
void FieldStrength::construct_Fmunu_1x2_traceless(Field_G& Fmunu_1x2,
                                                  const int mu, const int nu, const Field_G& U)
{
  const int Nvol = CommonParameters::Nvol();

  //-- building blocks
  //      mu (2)
  // (1)  +->-+
  //  nu  |   |
  //     i+   +
  Field_G Cup1;

  m_staple.upper(Cup1, U, mu, nu);

  //      +-<-+ (2)
  //          | nu
  //     i+->-+
  //      mu (1)
  Field_G Cup2;
  m_staple.upper(Cup2, U, nu, mu);

  // (1) i+   +
  //  nu  |   |
  //      +->-+
  //      mu (2)
  Field_G Cdn1;
  m_staple.lower(Cdn1, U, mu, nu);

  // (2)  +->-+
  //  nu  |
  //      +-<-+i
  //      mu (1)
  Field_G Cdn2;
  m_staple.lower(Cdn2, U, nu, mu);

  Field_G Umu;
  Umu.setpart_ex(0, U, mu);

  Field_G Unu;
  Unu.setpart_ex(0, U, nu);

  //  +->-+
  //
  // i+
  Field_G Umu_nu;
  m_shift.backward(Umu_nu, Umu, nu);

  //      +
  //      |
  // i+   +
  Field_G Unu_mu;
  m_shift.backward(Unu_mu, Unu, mu);

  //---------------------------------
  //-- 2x1 part

  //- upper-right(2x1)
  //    +---<---+      +      +-<-+         <---+
  // nu |       |  =   |  +          +          |  +
  //   i+--->---+     i+     i+         i+  >---+     i+->-+
  Field_G c;
  m_shift.backward(c, Cup2, mu);

  Field_G v;
  mult_Field_Gdd(v, 0, Umu_nu, 0, Unu, 0);

  Field_G w;
  mult_Field_Gnn(w, 0, c, 0, v, 0);

  Field_G rect;
  mult_Field_Gnn(rect, 0, Umu, 0, w, 0);

  Fmunu_1x2 = rect;

  //- lower-right(2x1)
  // (1) +---<---+     (1) i+---<---+
  //  nu |       |  ->  nu  |       |
  //    i+--->---+          +--->---+
  //      mu (2)             mu (2)
  mult_Field_Gnd(v, 0, c, 0, Umu_nu, 0);
  mult_Field_Gnn(w, 0, Umu, 0, v, 0);

  mult_Field_Gdn(rect, 0, Unu, 0, w, 0);

  m_shift.forward(v, rect, nu);

  axpy(Fmunu_1x2, 1.0, v);

  //- upper-left(2x1)
  //      mu (2)
  //  +---<---+ (1)              +-<-+      +-<-+           +
  //  |       | nu  =         +  |       +         +        |
  //  +---i---+       i+->-+     +->-+i     +i         +i   +
  //
  // NB. shift.forward(mu)
  //      mu (2)             mu (2)
  //  +---<---+ (1)      +---<---+ (1)
  //  |       | nu  ->   |       | nu
  //  +---i---+          +--->---+i
  mult_Field_Gdn(v, 0, Cdn2, 0, Umu, 0);
  mult_Field_Gdn(w, 0, Umu_nu, 0, v, 0);

  mult_Field_Gnn(rect, 0, Unu_mu, 0, w, 0);

  m_shift.forward(v, rect, mu);

  axpy(Fmunu_1x2, 1.0, v);

  //- lower-left(2x1)
  //      mu (1)
  //  +---<---+ (2)        +                +-<-+      +-<-+
  //  |       | nu  =      |  +          +  |       +
  //  +---i---+       +i   +     i+->-+     +->-+i     +i
  //
  // NB. shift.forward(mu+nu)
  //      mu (1)             mu (1)            mu (1)
  //  +---<---+ (2)      +---<---+ (2)     +---<---+i (2)
  //  |       | nu  ->   |       | nu  ->  |       |  nu
  //  +---i---+          +--->---+i        +--->---+
  mult_Field_Gnn(v, 0, Umu, 0, Unu_mu, 0);
  mult_Field_Gdn(w, 0, Cdn2, 0, v, 0);

  mult_Field_Gdn(rect, 0, Umu_nu, 0, w, 0);

  m_shift.forward(w, rect, mu);
  m_shift.forward(v, w, nu);

  axpy(Fmunu_1x2, 1.0, v);
  //---------------------------------


  //---------------------------------
  //-- 1x2 part

  //- upper-right(1x2)
  //     +-<-+             +-<-+
  //     |   |             |   |
  // (2) +   +  =  +   +   +   +   +      +   +
  //  nu |   |     |                      |
  //    i+->-+    i+      i+         i+   +      i+->-+
  //     mu (1)
  m_shift.backward(c, Cup1, nu);

  mult_Field_Gdd(v, 0, c, 0, Unu, 0);
  mult_Field_Gnn(w, 0, Unu_mu, 0, v, 0);

  mult_Field_Gnn(rect, 0, Umu, 0, w, 0);

  axpy(Fmunu_1x2, 1.0, rect);

  //- lower-right(1x2)
  //      mu (2)
  // (1)  +-<-+      +-<-+           +                     +
  //  nu  |   |                      |                     |
  //     i+   +  =  i+      +   i+   +   +   i+   +   +   i+
  //      |   |                               |   |
  //      +->-+                               +->-+
  //
  // NB. shift.forward(nu)
  //      mu (2)        mu (2)
  // (1)  +-<-+    (1) i+-<-+
  //  nu  |   |     nu  |   |
  //     i+   +  ->     +   +
  //      |   |         |   |
  //      +->-+         +->-+
  mult_Field_Gnd(v, 0, Unu_mu, 0, Umu_nu, 0);
  mult_Field_Gnn(w, 0, Cdn1, 0, v, 0);

  mult_Field_Gdn(rect, 0, Unu, 0, w, 0);

  m_shift.forward(v, rect, nu);

  axpy(Fmunu_1x2, 1.0, v);

  //- upper-left(1x2)
  //  +-<-+                            +-<-+
  //  |   |                            |   |
  //  +   + (1)  =         +   +   +   +   +   +       +
  //  |   | nu                 |                       |
  // i+->-+        i+->-+     i+      i+          i+   +
  //  mu (2)
  //
  // NB. shift.forward(mu)
  //  +-<-+         +-<-+
  //  |   |         |   |
  //  +   + (1) ->  +   + (1)
  //  |   | nu      |   | nu
  // i+->-+         +->-+i
  //  mu (2)        mu (2)
  mult_Field_Gdn(v, 0, Unu, 0, Umu, 0);
  mult_Field_Gdn(w, 0, c, 0, v, 0);

  mult_Field_Gnn(rect, 0, Unu_mu, 0, w, 0);

  m_shift.forward(v, rect, mu);

  axpy(Fmunu_1x2, 1.0, v);

  //- lower-left(1x2)
  //      mu (1)
  // (2)  +-<-+          +                     +        +-<-+
  //  nu  |   |          |                     |
  //     i+   +  =  i+   +   +   i+   +   +   i+   +   i+
  //      |   |                   |   |
  //      +->-+                   +->-+
  //
  // NB. shift.forward(mu+nu)
  //      mu (1)        mu (1)
  // (2)  +-<-+    (2)  +-<-+i
  //  nu  |   |     nu  |   |
  //     i+   +  ->     +   +
  //      |   |         |   |
  //      +->-+         +->-+
  mult_Field_Gnn(v, 0, Cdn1, 0, Unu_mu, 0);
  mult_Field_Gdn(w, 0, Unu, 0, v, 0);

  mult_Field_Gdn(rect, 0, Umu_nu, 0, w, 0);

  m_shift.forward(w, rect, mu);
  m_shift.forward(v, w, nu);

  axpy(Fmunu_1x2, 1.0, v);
  //---------------------------------

  // ah_Field_G(Fmunu_1x2, 0);
  at_Field_G(Fmunu_1x2, 0);

  //- normalization = 1/8
  //  i.e. overall factor = 1/4, average of 1x2 and 2x1 = 1/2
  scal(Fmunu_1x2, -0.125);
  Fmunu_1x2.xI();
}


//====================================================================
void FieldStrength::construct_Fmunu_plaq_traceless(Field_G& Fmunu_plaq,
                                                   const int mu, const int nu, const Field_G& U)
{
  const int Nvol = CommonParameters::Nvol();

  //- building blocks
  // (1)  mu (2)
  //    +->-+
  // nu |   |
  //   i+   +
  Field_G Cup;

  m_staple.upper(Cup, U, mu, nu);

  Field_G Umu;
  Umu.setpart_ex(0, U, mu);

  //- right clover leaf
  // (2)  +-<-+
  //  nu  |   |
  //     i+->-+
  //      mu (1)
  mult_Field_Gnd(Fmunu_plaq, 0, Umu, 0, Cup, 0);

  // ah_Field_G(Fmunu_plaq, 0);
  at_Field_G(Fmunu_plaq, 0);

  Fmunu_plaq.xI();
}


//====================================================================
double FieldStrength::contract(const Field_G& Fmunu_1, const Field_G& Fmunu_2)
{
  const int Nvol = CommonParameters::Nvol();

  double FF = 0.0;

  for (int site = 0; site < Nvol; ++site) {
    FF += ReTr(Fmunu_1.mat(site) * Fmunu_2.mat(site));
  }
  const double result = Communicator::reduce_sum(FF);

  return result;
}


//====================================================================
void FieldStrength::contract_at_t(std::vector<double>& corr_global,
                                  const Field_G& Fmunu_1, const Field_G& Fmunu_2)
{
  const int Nx = CommonParameters::Nx();
  const int Ny = CommonParameters::Ny();
  const int Nz = CommonParameters::Nz();
  const int Nt = CommonParameters::Nt();
  const int Lt = CommonParameters::Lt();

  assert(corr_global.size() == Lt);

  std::vector<double> corr_local(Nt, 0.0);
  for (int t = 0; t < Nt; ++t) {
    for (int z = 0; z < Nz; ++z) {
      for (int y = 0; y < Ny; ++y) {
        for (int x = 0; x < Nx; ++x) {
          int site = x + Nx * (y + Ny * (z + Nz * t));

          corr_local[t] += ReTr(Fmunu_1.mat(site) * Fmunu_2.mat(site));
        }
      }
    }
  }
  global_corr_t(corr_global, corr_local);
}


//====================================================================
void FieldStrength::contract_at_t(std::vector<double>& corr_global,
                                  const Field_G& Fmunu_1,
                                  const Field_G& Fmunu_2,
                                  const std::vector<int>& momentum_sink,
                                  const std::vector<int>& source_position)
{
  int       Nx = CommonParameters::Nx();
  int       Ny = CommonParameters::Ny();
  int       Nz = CommonParameters::Nz();
  int       Nt = CommonParameters::Nt();
  const int Lx = CommonParameters::Lx();
  const int Ly = CommonParameters::Ly();
  const int Lz = CommonParameters::Lz();
  const int Lt = CommonParameters::Lt();
  const int ND = CommonParameters::Nd();

  assert(corr_global.size() == Lt);
  assert(momentum_sink.size() == ND - 1);
  assert(source_position.size() == ND);

  static const double PI = 4.0 * atan(1.0);
  std::vector<double> p_unit(ND - 1);
  p_unit[0] = (2.0 * PI / Lx) * momentum_sink[0];
  p_unit[1] = (2.0 * PI / Ly) * momentum_sink[1];
  p_unit[2] = (2.0 * PI / Lz) * momentum_sink[2];

  std::vector<int> ipe(ND - 1);
  ipe[0] = Communicator::ipe(0);
  ipe[1] = Communicator::ipe(1);
  ipe[2] = Communicator::ipe(2);

  std::vector<double> corr_local(Nt, 0);
  for (int t = 0; t < Nt; ++t) {
    for (int z = 0; z < Nz; ++z) {
      for (int y = 0; y < Ny; ++y) {
        for (int x = 0; x < Nx; ++x) {
          int    site      = x + Nx * (y + Ny * (z + Nz * t));
          int    x_global  = x + ipe[0] * Nx;
          int    y_global  = y + ipe[1] * Ny;
          int    z_global  = z + ipe[2] * Nz;
          double p_x       = p_unit[0] * (x_global - source_position[0]);
          double p_y       = p_unit[1] * (y_global - source_position[1]);
          double p_z       = p_unit[2] * (z_global - source_position[2]);
          double cos_p_xyz = cos(p_x + p_y + p_z);

          double FF = ReTr(Fmunu_1.mat(site) * Fmunu_2.mat(site));
          corr_local[t] += FF * cos_p_xyz;
        }
      }
    }
  }
  global_corr_t(corr_global, corr_local);
}


//====================================================================
void FieldStrength::contract_at_x(std::vector<double>& corr_global,
                                  const Field_G& Fmunu_1,
                                  const Field_G& Fmunu_2)
{
  const int Nx = CommonParameters::Nx();
  const int Ny = CommonParameters::Ny();
  const int Nz = CommonParameters::Nz();
  const int Nt = CommonParameters::Nt();
  const int Lx = CommonParameters::Lx();

  assert(corr_global.size() == Lx);

  std::vector<double> corr_local(Nx, 0.0);
  for (int x = 0; x < Nx; ++x) {
    for (int t = 0; t < Nt; ++t) {
      for (int z = 0; z < Nz; ++z) {
        for (int y = 0; y < Ny; ++y) {
          int site = x + Nx * (y + Ny * (z + Nz * t));

          corr_local[x] += ReTr(Fmunu_1.mat(site) * Fmunu_2.mat(site));
        }
      }
    }
  }
  global_corr_x(corr_global, corr_local);
}


//====================================================================
void FieldStrength::contract_at_x(std::vector<double>& corr_global,
                                  const Field_G& Fmunu_1,
                                  const Field_G& Fmunu_2,
                                  const std::vector<int>& momentum_sink,
                                  const std::vector<int>& source_position)
{
  const int Nx = CommonParameters::Nx();
  const int Ny = CommonParameters::Ny();
  const int Nz = CommonParameters::Nz();
  const int Nt = CommonParameters::Nt();
  const int Lx = CommonParameters::Lx();
  const int Ly = CommonParameters::Ly();
  const int Lz = CommonParameters::Lz();
  const int Lt = CommonParameters::Lt();
  const int ND = CommonParameters::Nd();

  assert(corr_global.size() == Lx);
  assert(momentum_sink.size() == ND - 1);
  assert(source_position.size() == ND);

  static const double PI = 4.0 * atan(1.0);
  std::vector<double> p_unit(ND - 1);
  p_unit[0] = (2.0 * PI / Ly) * momentum_sink[0];
  p_unit[1] = (2.0 * PI / Lz) * momentum_sink[1];
  p_unit[2] = (2.0 * PI / Lt) * momentum_sink[2];

  std::vector<int> ipe(ND);
  ipe[0] = Communicator::ipe(0);
  ipe[1] = Communicator::ipe(1);
  ipe[2] = Communicator::ipe(2);
  ipe[3] = Communicator::ipe(3);

  std::vector<double> corr_local(Nx, 0);
  for (int x = 0; x < Nx; ++x) {
    for (int t = 0; t < Nt; ++t) {
      for (int z = 0; z < Nz; ++z) {
        for (int y = 0; y < Ny; ++y) {
          int    site      = x + Nx * (y + Ny * (z + Nz * t));
          int    y_global  = y + ipe[1] * Ny;
          int    z_global  = z + ipe[2] * Nz;
          int    t_global  = t + ipe[3] * Nt;
          double p_y       = p_unit[0] * (y_global - source_position[1]);
          double p_z       = p_unit[1] * (z_global - source_position[2]);
          double p_t       = p_unit[2] * (t_global - source_position[3]);
          double cos_p_yzt = cos(p_t + p_y + p_z);
          double FF        = ReTr(Fmunu_1.mat(site) * Fmunu_2.mat(site));
          corr_local[x] += FF * cos_p_yzt;
        }
      }
    }
  }
  global_corr_x(corr_global, corr_local);
}


//====================================================================
void FieldStrength::contract_at_y(std::vector<double>& corr_global,
                                  const Field_G& Fmunu_1,
                                  const Field_G& Fmunu_2)
{
  const int Nx = CommonParameters::Nx();
  const int Ny = CommonParameters::Ny();
  const int Nz = CommonParameters::Nz();
  const int Nt = CommonParameters::Nt();
  const int Ly = CommonParameters::Ly();

  assert(corr_global.size() == Ly);

  std::vector<double> corr_local(Ny, 0.0);
  for (int y = 0; y < Ny; ++y) {
    for (int t = 0; t < Nt; ++t) {
      for (int z = 0; z < Nz; ++z) {
        for (int x = 0; x < Nx; ++x) {
          int site = x + Nx * (y + Ny * (z + Nz * t));

          corr_local[y] += ReTr(Fmunu_1.mat(site) * Fmunu_2.mat(site));
        }
      }
    }
  }
  global_corr_y(corr_global, corr_local);
}


//====================================================================
void FieldStrength::contract_at_y(std::vector<double>& corr_global,
                                  const Field_G& Fmunu_1,
                                  const Field_G& Fmunu_2,
                                  const std::vector<int>& momentum_sink,
                                  const std::vector<int>& source_position)
{
  const int Nx = CommonParameters::Nx();
  const int Ny = CommonParameters::Ny();
  const int Nz = CommonParameters::Nz();
  const int Nt = CommonParameters::Nt();
  const int Lx = CommonParameters::Lx();
  const int Ly = CommonParameters::Ly();
  const int Lz = CommonParameters::Lz();
  const int Lt = CommonParameters::Lt();
  const int ND = CommonParameters::Nd();

  assert(corr_global.size() == Ly);
  assert(momentum_sink.size() == ND - 1);
  assert(source_position.size() == ND);

  static const double PI = 4.0 * atan(1.0);
  std::vector<double> p_unit(ND - 1);
  p_unit[0] = (2.0 * PI / Ly) * momentum_sink[0];
  p_unit[1] = (2.0 * PI / Lz) * momentum_sink[1];
  p_unit[2] = (2.0 * PI / Lt) * momentum_sink[2];

  std::vector<int> ipe(ND);
  ipe[0] = Communicator::ipe(0);
  ipe[1] = Communicator::ipe(1);
  ipe[2] = Communicator::ipe(2);
  ipe[3] = Communicator::ipe(3);

  std::vector<double> corr_local(Ny, 0);
  for (int y = 0; y < Ny; ++y) {
    for (int t = 0; t < Nt; ++t) {
      for (int z = 0; z < Nz; ++z) {
        for (int x = 0; x < Nx; ++x) {
          int    site      = x + Nx * (y + Ny * (z + Nz * t));
          int    x_global  = x + ipe[0] * Nx;
          int    z_global  = z + ipe[2] * Nz;
          int    t_global  = t + ipe[3] * Nt;
          double p_x       = p_unit[0] * (x_global - source_position[0]);
          double p_z       = p_unit[1] * (z_global - source_position[2]);
          double p_t       = p_unit[2] * (t_global - source_position[3]);
          double cos_p_xzt = cos(p_t + p_x + p_z);
          double FF        = ReTr(Fmunu_1.mat(site) * Fmunu_2.mat(site));
          corr_local[y] += FF * cos_p_xzt;
        }
      }
    }
  }
  global_corr_y(corr_global, corr_local);
}


//====================================================================
void FieldStrength::contract_at_z(std::vector<double>& corr_global,
                                  const Field_G& Fmunu_1,
                                  const Field_G& Fmunu_2)
{
  const int Nx = CommonParameters::Nx();
  const int Ny = CommonParameters::Ny();
  const int Nz = CommonParameters::Nz();
  const int Nt = CommonParameters::Nt();
  const int Lz = CommonParameters::Lz();

  assert(corr_global.size() == Lz);

  std::vector<double> corr_local(Nz, 0.0);
  for (int z = 0; z < Nz; ++z) {
    for (int t = 0; t < Nt; ++t) {
      for (int y = 0; y < Ny; ++y) {
        for (int x = 0; x < Nx; ++x) {
          int site = x + Nx * (y + Ny * (z + Nz * t));

          corr_local[z] += ReTr(Fmunu_1.mat(site) * Fmunu_2.mat(site));
        }
      }
    }
  }
  global_corr_z(corr_global, corr_local);
}


//====================================================================
void FieldStrength::contract_at_z(std::vector<double>& corr_global,
                                  const Field_G& Fmunu_1,
                                  const Field_G& Fmunu_2,
                                  const std::vector<int>& momentum_sink,
                                  const std::vector<int>& source_position)
{
  const int Nx = CommonParameters::Nx();
  const int Ny = CommonParameters::Ny();
  const int Nz = CommonParameters::Nz();
  const int Nt = CommonParameters::Nt();
  const int Lx = CommonParameters::Lx();
  const int Ly = CommonParameters::Ly();
  const int Lz = CommonParameters::Lz();
  const int Lt = CommonParameters::Lt();
  const int ND = CommonParameters::Nd();

  assert(corr_global.size() == Lz);
  assert(momentum_sink.size() == ND - 1);
  assert(source_position.size() == ND);

  static const double PI = 4.0 * atan(1.0);
  std::vector<double> p_unit(ND - 1);
  p_unit[0] = (2.0 * PI / Ly) * momentum_sink[0];
  p_unit[1] = (2.0 * PI / Lz) * momentum_sink[1];
  p_unit[2] = (2.0 * PI / Lt) * momentum_sink[2];

  std::vector<int> ipe(ND);
  ipe[0] = Communicator::ipe(0);
  ipe[1] = Communicator::ipe(1);
  ipe[2] = Communicator::ipe(2);
  ipe[3] = Communicator::ipe(3);

  std::vector<double> corr_local(Nz, 0);
  for (int z = 0; z < Nz; ++z) {
    for (int t = 0; t < Nt; ++t) {
      for (int y = 0; y < Ny; ++y) {
        for (int x = 0; x < Nx; ++x) {
          int    site      = x + Nx * (y + Ny * (z + Nz * t));
          int    x_global  = x + ipe[0] * Nx;
          int    y_global  = y + ipe[1] * Ny;
          int    t_global  = t + ipe[3] * Nt;
          double p_x       = p_unit[0] * (x_global - source_position[0]);
          double p_y       = p_unit[1] * (y_global - source_position[1]);
          double p_t       = p_unit[2] * (t_global - source_position[3]);
          double cos_p_xyt = cos(p_t + p_x + p_y);
          double FF        = ReTr(Fmunu_1.mat(site) * Fmunu_2.mat(site));
          corr_local[z] += FF * cos_p_xyt;
        }
      }
    }
  }
  global_corr_z(corr_global, corr_local);
}


//====================================================================
void FieldStrength::global_corr_x(std::vector<double>& corr_global,
                                  const std::vector<double>& corr_local)
{
  const int Lx = CommonParameters::Lx();
  const int Nx = CommonParameters::Nx();

  assert(corr_global.size() == Lx);
  assert(corr_local.size() == Nx);

  const int ipe_x = Communicator::ipe(0);

  std::vector<double> corr_tmp(Lx, 0.0);

  for (int x = 0; x < Nx; ++x) {
    int x_global = x + ipe_x * Nx;
    corr_tmp[x_global] = corr_local[x];
  }

  for (int x_global = 0; x_global < Lx; ++x_global) {
    double cr_r = Communicator::reduce_sum(corr_tmp[x_global]);
    corr_global[x_global] = cr_r;
  }
}


//====================================================================
void FieldStrength::global_corr_y(std::vector<double>& corr_global,
                                  const std::vector<double>& corr_local)
{
  const int Ly = CommonParameters::Ly();
  const int Ny = CommonParameters::Ny();

  assert(corr_global.size() == Ly);
  assert(corr_local.size() == Ny);

  const int ipe_y = Communicator::ipe(1);

  std::vector<double> corr_tmp(Ly, 0.0);

  for (int y = 0; y < Ny; ++y) {
    int y_global = y + ipe_y * Ny;
    corr_tmp[y_global] = corr_local[y];
  }

  for (int y_global = 0; y_global < Ly; ++y_global) {
    double cr_r = Communicator::reduce_sum(corr_tmp[y_global]);
    corr_global[y_global] = cr_r;
  }
}


//====================================================================
void FieldStrength::global_corr_z(std::vector<double>& corr_global,
                                  const std::vector<double>& corr_local)
{
  const int Lz = CommonParameters::Lz();
  const int Nz = CommonParameters::Nz();

  assert(corr_global.size() == Lz);
  assert(corr_local.size() == Nz);

  const int ipe_z = Communicator::ipe(2);

  std::vector<double> corr_tmp(Lz, 0.0);

  for (int z = 0; z < Nz; ++z) {
    int z_global = z + ipe_z * Nz;
    corr_tmp[z_global] = corr_local[z];
  }

  for (int z_global = 0; z_global < Lz; ++z_global) {
    double cr_r = Communicator::reduce_sum(corr_tmp[z_global]);
    corr_global[z_global] = cr_r;
  }
}


//====================================================================
void FieldStrength::global_corr_t(std::vector<double>& corr_global,
                                  const std::vector<double>& corr_local)
{
  const int Lt = CommonParameters::Lt();
  const int Nt = CommonParameters::Nt();

  assert(corr_global.size() == Lt);
  assert(corr_local.size() == Nt);

  const int ipe_t = Communicator::ipe(3);

  std::vector<double> corr_tmp(Lt, 0.0);

  for (int t = 0; t < Nt; ++t) {
    int t_global = t + ipe_t * Nt;
    corr_tmp[t_global] = corr_local[t];
  }

  for (int t_global = 0; t_global < Lt; ++t_global) {
    double cr_r = Communicator::reduce_sum(corr_tmp[t_global]);
    corr_global[t_global] = cr_r;
  }
}


//====================================================================
//============================================================END=====
