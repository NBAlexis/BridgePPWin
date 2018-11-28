#include "BridgeLib_Private.h"

/*!
        @file    $Id:: fieldStrength.cpp #$

        @brief

        @author  Yusuke Namekawa  (namekawa)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 00:49:46 #$

        @version $LastChangedRevision: 1561 $
*/

#include "fieldStrength.h"

const std::string FieldStrength::class_name = "FieldStrengh";

//====================================================================
void FieldStrength::construct_Fmunu_1x1(Field_G& Fmunu_1x1,
                                        const int mu, const int nu, const Field_G& U)
{
  int Nvol = CommonParameters::Nvol();

  Field_G Cup(Nvol, 1), Cdn(Nvol, 1);
  Field_G Umu(Nvol, 1);
  Field_G v(Nvol, 1), v2(Nvol, 1);


  //- building blocks
  // (1)  mu (2)
  //    +->-+
  // nu |   |
  //   i+   +
  m_staple.upper(Cup, U, mu, nu);

  // (1)  mu (2)
  //   i+   +
  // nu |   |
  //    +->-+
  m_staple.lower(Cdn, U, mu, nu);
  Umu.setpart_ex(0, U, mu);


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


  //- left clover leaf
  //   mu (2)
  //   +-<-+ (1)
  //   |   | nu    mu (1)
  //   +->-+i   +  +-<-+i (2)
  //               |   | nu
  //               +->-+
  mult_Field_Gdn(v, 0, Cup, 0, Umu, 0);
  multadd_Field_Gdn(v, 0, Umu, 0, Cdn, 0, 1.0);

  // NB. shift.forward(mu)
  //   mu (2)        mu (2)
  //   +-<-+ (1)     +-<-+ (1)
  //   |   | nu  ->  |   | nu
  //  i+->-+         +->-+i
  m_shift.forward(v2, v, mu);

  axpy(Fmunu_1x1, 1.0, v2);

  ah_Field_G(Fmunu_1x1, 0);

  scal(Fmunu_1x1, -0.25);
  Fmunu_1x1.xI();
}


//====================================================================
void FieldStrength::construct_Fmunu_1x2(Field_G& Fmunu_1x2,
                                        const int mu, const int nu, const Field_G& U)
{
  int Nvol = CommonParameters::Nvol();

  Field_G Cup1(Nvol, 1), Cup2(Nvol, 1);
  Field_G Cdn1(Nvol, 1), Cdn2(Nvol, 1);
  Field_G Umu(Nvol, 1), Unu(Nvol, 1);
  Field_G Umu_nu(Nvol, 1), Unu_mu(Nvol, 1);
  Field_G rect(Nvol, 1);
  Field_G v(Nvol, 1), w(Nvol, 1), c(Nvol, 1);

  //-- building blocks
  //      mu (2)
  // (1)  +->-+
  //  nu  |   |
  //     i+   +
  m_staple.upper(Cup1, U, mu, nu);

  //      +-<-+ (2)
  //          | nu
  //     i+->-+
  //      mu (1)
  m_staple.upper(Cup2, U, nu, mu);

  // (1) i+   +
  //  nu  |   |
  //      +->-+
  //      mu (2)
  m_staple.lower(Cdn1, U, mu, nu);

  // (2)  +->-+
  //  nu  |
  //      +-<-+i
  //      mu (1)
  m_staple.lower(Cdn2, U, nu, mu);

  Umu.setpart_ex(0, U, mu);
  Unu.setpart_ex(0, U, nu);

  //  +->-+
  //
  // i+
  m_shift.backward(Umu_nu, Umu, nu);

  //      +
  //      |
  // i+   +
  m_shift.backward(Unu_mu, Unu, mu);

  //============================================================
  //-- 2x1 part

  //- upper-right(2x1)
  //    +---<---+      +      +-<-+         <---+
  // nu |       |  =   |  +          +          |  +
  //   i+--->---+     i+     i+         i+  >---+     i+->-+
  m_shift.backward(c, Cup2, mu);

  mult_Field_Gdd(v, 0, Umu_nu, 0, Unu, 0);
  mult_Field_Gnn(w, 0, c, 0, v, 0);

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
  //============================================================


  //============================================================
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
  //============================================================

  ah_Field_G(Fmunu_1x2, 0);

  //- normalization = 1/8
  //  i.e. overall factor = 1/4, average of 1x2 and 2x1 = 1/2
  scal(Fmunu_1x2, -0.125);
  Fmunu_1x2.xI();
}


//====================================================================
//============================================================END=====
