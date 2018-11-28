#include "BridgeLib_Private.h"

/*!
        @file    $Id:: langevin_Momentum.cpp #$

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 00:49:46 #$

        @version $LastChangedRevision: 1561 $
*/

#include "langevin_Momentum.h"

//====================================================================
double Langevin_Momentum::set_iP(Field_G& iP)
{
  int Nc = CommonParameters::Nc();

  double iP2;

  if (Nc == 3) {
    iP2 = set_iP_SU3(iP);
    //  iP2 = set_iP_general_SU_N(iP);
    //  iP2 = set_iP_SU3_alt(iP);
    // Comment:
    //  if you want to check that set_iP_general_SU_N() works
    //  correctly for Nc=3 case, compare the result with
    //  set_iP_SU3_alt().          [13 Feb 2013 H.Matsufuru]
  } else {
    iP2 = set_iP_general_SU_N(iP);
  }

  return iP2;
}


//====================================================================
double Langevin_Momentum::set_iP_general_SU_N(Field_G& iP)
{
  int Nin  = iP.nin();
  int Nvol = iP.nvol(); // local volume (SA)
  int Nex  = iP.nex();

  int Nc  = CommonParameters::Nc();
  int NcA = Nc * Nc - 1;

  vout.general(m_vl, "  Conjugate momenta for SU(%d) gauge field.\n", Nc);

  GeneratorSet_Mat_SU_N   gen_set(Nc);
  std::vector<Mat_SU_N *> iT(NcA);
  for (int i = 0; i < NcA; ++i) {
    iT[i]  = new Mat_SU_N(Nc);
    *iT[i] = gen_set.get_generator(i);
    iT[i]->xI();
  }
  Mat_SU_N u_tmp(Nc), u_tmp2(Nc);

  Field Hrand(NcA, Nvol, Nex);     // Local random number (SA)
  m_rand->gauss_lex_global(Hrand); // generate gaussian random number Hrand. (SA)

  for (int ex = 0; ex < Nex; ++ex) {
    for (int site = 0; site < Nvol; ++site) {
      u_tmp = 0.0;
      for (int i = 0; i < NcA; ++i) {
        u_tmp2  = *iT[i];
        u_tmp2 *= 2.0 * Hrand.cmp(i, site, ex);
        u_tmp  += u_tmp2;
      }

      iP.set_mat(site, ex, u_tmp);
    }
  }

  for (int i = 0; i < NcA; ++i) {
    delete iT[i];
  }

  double iP2 = iP.norm(); //calculate global norm sqrt(|iP2|^2) (SA)
  iP2 = 0.5 * iP2 * iP2;  // |iP|^2/2 (SA)

  return iP2;
}


//====================================================================
double Langevin_Momentum::set_iP_SU3_alt(Field_G& iP)
{
  // alternative implementation to set_iP_SU3():
  // this gives the same result as set_iP_general_SU_N() for Nc=3 case.

  int Nin  = iP.nin();
  int Nvol = iP.nvol(); // local volume (SA)
  int Nex  = iP.nex();

  int Nc  = CommonParameters::Nc();
  int NcA = Nc * Nc - 1;

  // confirm that now gauge group is SU(3)
  assert(NcA == 8);

  Field Hrand(NcA, Nvol, Nex);     // Local random number (SA)

  m_rand->gauss_lex_global(Hrand); // generate gaussian random number Hrand. (SA)

  double sq3r  = 1.0 / sqrt(3.0);
  double sq3r2 = 2.0 * sq3r;

  for (int ex = 0; ex < Nex; ++ex) {
    for (int site = 0; site < Nvol; ++site) {
      //here SU(3) is explicitly assumed. need to generalize.
      // return (jin,site,ex) component(double) of Hrand. (SA)
      double hc1 = Hrand.cmp(0, site, ex);
      double hc2 = Hrand.cmp(1, site, ex);
      double hc4 = Hrand.cmp(2, site, ex);
      double hc5 = Hrand.cmp(3, site, ex);
      double hc6 = Hrand.cmp(4, site, ex);
      double hc7 = Hrand.cmp(5, site, ex);
      double hc3 = Hrand.cmp(6, site, ex);
      double hc8 = Hrand.cmp(7, site, ex);

      iP.set(0, site, ex, 0.0);
      iP.set(1, site, ex, hc3 + hc8 * sq3r);
      iP.set(2, site, ex, hc2);
      iP.set(3, site, ex, hc1);
      iP.set(4, site, ex, hc5);
      iP.set(5, site, ex, hc4);
      iP.set(6, site, ex, -hc2);
      iP.set(7, site, ex, hc1);
      iP.set(8, site, ex, 0.0);
      iP.set(9, site, ex, -hc3 + hc8 * sq3r);
      iP.set(10, site, ex, hc7);
      iP.set(11, site, ex, hc6);
      iP.set(12, site, ex, -hc5);
      iP.set(13, site, ex, hc4);
      iP.set(14, site, ex, -hc7);
      iP.set(15, site, ex, hc6);
      iP.set(16, site, ex, 0.0);
      iP.set(17, site, ex, -hc8 * sq3r2);
    }
  }

  double iP2 = iP.norm(); //calculate global norm sqrt(|iP2|^2) (SA)
  iP2 = 0.5 * iP2 * iP2;  // |iP|^2/2 (SA)

  return iP2;
}


//====================================================================
double Langevin_Momentum::set_iP_SU3(Field_G& iP)
{
  // implementation for SU(3) case: original version.

  int Nin  = iP.nin();
  int Nvol = iP.nvol(); // local volume (SA)
  int Nex  = iP.nex();

  int Nc  = CommonParameters::Nc();
  int NcA = Nc * Nc - 1;

  // confirm that now gauge group is SU(3)
  assert(NcA == 8);

  Field Hrand(NcA, Nvol, Nex);     // Local random number (SA)

  m_rand->gauss_lex_global(Hrand); // generate gaussian random number Hrand. (SA)

  double sq3r  = 1.0 / sqrt(3.0);
  double sq3r2 = 2.0 * sq3r;

  for (int ex = 0; ex < Nex; ++ex) {
    for (int site = 0; site < Nvol; ++site) {
      //here SU(3) is explicitly assumed. need to generalize.
      // return (jin,site,ex) component(double) of Hrand. (SA)
      double hc1 = Hrand.cmp(0, site, ex);
      double hc2 = Hrand.cmp(1, site, ex);
      double hc3 = Hrand.cmp(2, site, ex);
      double hc4 = Hrand.cmp(3, site, ex);
      double hc5 = Hrand.cmp(4, site, ex);
      double hc6 = Hrand.cmp(5, site, ex);
      double hc7 = Hrand.cmp(6, site, ex);
      double hc8 = Hrand.cmp(7, site, ex);

      /*
      iP = i P^a T^a: T^a: Gellmann matrix
       P=
       | hc3+hc8/sqrt(3), hc1-i*hc2,       hc4-i*hc5     |
       | hc1+i*hc2,      -hc3+hc8/sqrt(3), hc6-i*hc7     |
       | hc4+i*hc5,       hc6+i*hc7,      -2*hc8/sqrt(3) |
      */
      iP.set(0, site, ex, 0.0);
      iP.set(1, site, ex, hc3 + hc8 * sq3r);
      iP.set(2, site, ex, hc2);
      iP.set(3, site, ex, hc1);
      iP.set(4, site, ex, hc5);
      iP.set(5, site, ex, hc4);
      iP.set(6, site, ex, -hc2);
      iP.set(7, site, ex, hc1);
      iP.set(8, site, ex, 0.0);
      iP.set(9, site, ex, -hc3 + hc8 * sq3r);
      iP.set(10, site, ex, hc7);
      iP.set(11, site, ex, hc6);
      iP.set(12, site, ex, -hc5);
      iP.set(13, site, ex, hc4);
      iP.set(14, site, ex, -hc7);
      iP.set(15, site, ex, hc6);
      iP.set(16, site, ex, 0.0);
      iP.set(17, site, ex, -hc8 * sq3r2);
    }
  }

  double iP2 = iP.norm(); //calculate global norm sqrt(|iP2|^2) (SA)
  iP2 = 0.5 * iP2 * iP2;  // |iP|^2/2 (SA)

  return iP2;
}


//====================================================================
//============================================================END=====
