#include "BridgeLib_Private.h"

/*!
        @file    $Id:: projection_Stout_SU3.cpp #$

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-03-14 06:41:34 #$

        @version $LastChangedRevision: 1593 $
*/

#include "projection_Stout_SU3.h"


#ifdef USE_FACTORY
namespace {
  Projection *create_object()
  {
    return new Projection_Stout_SU3();
  }


  bool init = Projection::Factory::Register("Stout_SU3", create_object);
}
#endif



const std::string Projection_Stout_SU3::class_name = "Projection_Stout_SU3";

//====================================================================
void Projection_Stout_SU3::set_parameters(const Parameters& params)
{
  const string str_vlevel = params.get_string("verbose_level");

  m_vl = vout.set_verbose_level(str_vlevel);
}


//====================================================================
void Projection_Stout_SU3::setup()
{
  m_flop = 0;
  m_time = 0.0;

  assert(CommonParameters::Nc() == NC);
}


//====================================================================
void Projection_Stout_SU3::print_stat()
{
  double gflops = 1.e-9 * double(m_flop) / m_time;

  vout.general(m_vl, "%s:\n", class_name.c_str());
  vout.general(m_vl, "  total time: %f\n", m_time);
  vout.general(m_vl, "  total flop: %d\n", m_flop);
  vout.general(m_vl, "  GFlops    : %f\n", gflops);
}


//====================================================================
void Projection_Stout_SU3::project(Field_G& U,
                                   double ,
                                   const Field_G& Cst, const Field_G& Uorg)
{
  //  int id = 31;
  //  KEK_FopCountStart(id);
  double time0 = Communicator::get_time();

  // in stout projection, parameter alpha is dummy.

  int Nex  = Uorg.nex();
  int Nvol = Uorg.nvol();

  assert(Cst.nex() == Nex);
  assert(Cst.nvol() == Nvol);
  assert(U.nex() == Nex);
  assert(U.nvol() == Nvol);

  //int NinG = Uorg.nin();

  Mat_SU_N iQ0(NC), iQ1(NC), iQ2(NC), iQ3(NC);
  Mat_SU_N ct(NC), ut(NC), ut2(NC), e_iQ(NC);
  iQ0.unit();

  double   u, w;
  dcomplex f0, f1, f2;

  for (int mu = 0; mu < Nex; ++mu) {
    for (int site = 0; site < Nvol; ++site) {
      Uorg.mat(ut, site, mu);
      Cst.mat(ct, site, mu);
      iQ1.mult_nd(ct, ut);
      iQ1.at();
      iQ2.mult_nn(iQ1, iQ1);
      iQ3.mult_nn(iQ1, iQ2);

      double norm = iQ1.norm2();
      if (norm > 1.0e-10) {
        set_uw(u, w, iQ2, iQ3);
        set_fj(f0, f1, f2, u, w);

        for (int cc = 0; cc < NC * NC; ++cc) {
          dcomplex qt = f0 * cmplx(iQ0.r(cc), iQ0.i(cc))
                        + f1 * cmplx(iQ1.i(cc), -iQ1.r(cc))
                        - f2 * cmplx(iQ2.r(cc), iQ2.i(cc));
          e_iQ.set_r(cc, real(qt));
          e_iQ.set_i(cc, imag(qt));
        }
      } else {
        //  vout.general(m_vl,"project: |iQ1|^2 too small: %lf. Set e_iQ=1.\n",norm);
        e_iQ.unit();
      }

      ut2.mult_nn(e_iQ, ut);
      U.set_mat(site, mu, ut2);
    }
  }

  /*
  unsigned long count;
  double time;
  KEK_FopCountFinish(id,&count,&time);
  m_time += time;
  m_flop += count;
  */
  double time1 = Communicator::get_time();
  m_time += time1 - time0;
}


//====================================================================
void Projection_Stout_SU3::exp_iQ(Field_G& e_iQ, const Field_G& iQ)
{
  int Nvol = iQ.nvol();
  int Nex  = iQ.nex();

  Mat_SU_N iQ0(NC), iQ1(NC), iQ2(NC), iQ3(NC);

  iQ0.unit();

  double   u, w;
  dcomplex f0, f1, f2;

  for (int mu = 0; mu < Nex; ++mu) {
    for (int site = 0; site < Nvol; ++site) {
      iQ1 = iQ.mat(site, mu);
      iQ2 = iQ1 * iQ1;
      iQ3 = iQ1 * iQ2;

      double norm = iQ1.norm2();
      if (norm > 1.0e-10) {
        set_uw(u, w, iQ2, iQ3);
        set_fj(f0, f1, f2, u, w);

        for (int cc = 0; cc < NC * NC; ++cc) {
          dcomplex qt = f0 * cmplx(iQ0.r(cc), iQ0.i(cc))
                        + f1 * cmplx(iQ1.i(cc), -iQ1.r(cc))
                        - f2 * cmplx(iQ2.r(cc), iQ2.i(cc));
          e_iQ.set_ri(cc, site, mu, real(qt), imag(qt));
        }
      } else {
        //      vout.general(m_vl,"exp_iQ: |iQ1|^2 too small: %lf. Set e_iQ=1.\n",norm);
        e_iQ.set_mat(site, mu, iQ0);
      }
    }
  }
}


//====================================================================

/*!
<ul>
<li>See the implementation note "note_cloverHMC.pdf" (21 Mar 2012) by H.Matsufuru.
<li>Evaluate \f$\Sigma_\mu'(x)\exp(iQ_\mu(x))+C_\mu^\dagger i\Lambda_\mu(x)\f$ and \f$i\Lambda_\mu(x)U_\mu(x)\f$ in eq.(93)
<li>argument Sigmap \f$=\Sigma_\mu'(x)\f$ in eq.(93) in terms of (k)-th smearing.
<li>argument Cst \f$=C_\mu(x)\f$ of eq.(43) in terms of (k-1)-th smearing.
<li>argument Uorg \f$=U_\mu(x)\f$ in (k-1)-th smearing.
<li>argument Xi is a resultant \f$\Sigma_\mu'(x)\exp(iQ_\mu(x))+C_\mu^\dagger i\Lambda_\mu(x)\f$.
<li>argument iTheta is a resultant \f$i\Lambda_\mu(x)U_\mu(x)\f$.
<li>Comment by [Y.Taniguchi 2012.04.16]
</ul>
*/
//====================================================================
void Projection_Stout_SU3::force_recursive(Field_G& Xi, Field_G& iTheta,
                                           double , const Field_G& Sigmap,
                                           const Field_G& Cst, const Field_G& Uorg)
{
  // in stout projection, parameter alpha is dummy.

  //  int id = 31;
  //  KEK_FopCountStart(id);
  double time0 = Communicator::get_time();

  int Nvol = CommonParameters::Nvol();
  // int NinG = 2 * NC * NC;
  int Nex = Xi.nex();

  assert(Xi.nvol() == Nvol);
  assert(iTheta.nvol() == Nvol);
  assert(Sigmap.nvol() == Nvol);
  assert(Cst.nvol() == Nvol);
  assert(Uorg.nvol() == Nvol);
  assert(iTheta.nex() == Nex);
  assert(Sigmap.nex() == Nex);
  assert(Cst.nex() == Nex);
  assert(Uorg.nex() == Nex);

  Mat_SU_N C_tmp(NC), U_tmp(NC), Sigmap_tmp(NC);
  Mat_SU_N iQ0(NC), iQ1(NC), iQ2(NC), iQ3(NC), e_iQ(NC);
  Mat_SU_N B1(NC), B2(NC);
  Mat_SU_N USigmap(NC), iQUS(NC), iUSQ(NC), iGamma(NC);
  Mat_SU_N Xi_tmp(NC), iTheta_tmp(NC);
  Mat_SU_N tmp1(NC), tmp2(NC);
  iQ0.unit();

  double   u, w;
  dcomplex f0, f1, f2;

  for (int mu = 0; mu < Nex; ++mu) {
    for (int site = 0; site < Nvol; ++site) {
      //! C_tmp \f$=C_\mu(x)\f$
      Cst.mat(C_tmp, site, mu);
      //! U_tmp \f$=U_\mu(x)\f$
      Uorg.mat(U_tmp, site, mu);
      // Sigmap_tmp \f$=\Sigma_\mu'(x)\f$
      Sigmap.mat(Sigmap_tmp, site, mu);

      //! iQ1 \f$=iQ_\mu\f$
      iQ1.mult_nd(C_tmp, U_tmp);
      iQ1.at();
      iQ2.mult_nn(iQ1, iQ1);
      iQ3.mult_nn(iQ1, iQ2);

      // In order to aviod 1Q1=0
      double norm = iQ1.norm2();
      if (norm > 1.0e-10) {
        set_uw(u, w, iQ2, iQ3);
        set_fj(f0, f1, f2, u, w);

        for (int cc = 0; cc < NC * NC; ++cc) {
          dcomplex qt = f0 * cmplx(iQ0.r(cc), iQ0.i(cc))
                        + f1 * cmplx(iQ1.i(cc), -iQ1.r(cc))
                        - f2 * cmplx(iQ2.r(cc), iQ2.i(cc));
          e_iQ.set(cc, real(qt), imag(qt));
        }

        double xi0   = func_xi0(w);
        double xi1   = func_xi1(w);
        double u2    = u * u;
        double w2    = w * w;
        double cos_w = cos(w);

        dcomplex emiu = cmplx(cos(u), -sin(u));
        dcomplex e2iu = cmplx(cos(2.0 * u), sin(2.0 * u));

        dcomplex r01 = cmplx(2.0 * u, 2.0 * (u2 - w2)) * e2iu
                       + emiu * cmplx(16.0 * u * cos_w + 2.0 * u * (3.0 * u2 + w2) * xi0,
                                      -8.0 * u2 * cos_w + 2.0 * (9.0 * u2 + w2) * xi0);

        dcomplex r11 = cmplx(2.0, 4.0 * u) * e2iu
                       + emiu * cmplx(-2.0 * cos_w + (3.0 * u2 - w2) * xi0,
                                      2.0 * u * cos_w + 6.0 * u * xi0);

        dcomplex r21 = cmplx(0.0, 2.0) * e2iu
                       + emiu * cmplx(-3.0 * u * xi0, cos_w - 3.0 * xi0);

        dcomplex r02 = cmplx(-2.0, 0.0) * e2iu
                       + emiu * cmplx(-8.0 * u2 * xi0,
                                      2.0 * u * (cos_w + xi0 + 3.0 * u2 * xi1));

        dcomplex r12 = emiu * cmplx(2.0 * u * xi0,
                                    -cos_w - xi0 + 3.0 * u2 * xi1);

        dcomplex r22 = emiu * cmplx(xi0, -3.0 * u * xi1);

        double fden = 1.0 / (2 * (9.0 * u2 - w2) * (9.0 * u2 - w2));

        dcomplex b10 = cmplx(2.0 * u, 0.0) * r01 + cmplx(3.0 * u2 - w2, 0.0) * r02
                       - cmplx(30.0 * u2 + 2.0 * w2, 0.0) * f0;
        dcomplex b11 = cmplx(2.0 * u, 0.0) * r11 + cmplx(3.0 * u2 - w2, 0.0) * r12
                       - cmplx(30.0 * u2 + 2.0 * w2, 0.0) * f1;
        dcomplex b12 = cmplx(2.0 * u, 0.0) * r21 + cmplx(3.0 * u2 - w2, 0.0) * r22
                       - cmplx(30.0 * u2 + 2.0 * w2, 0.0) * f2;

        dcomplex b20 = r01 - cmplx(3.0 * u, 0.0) * r02 - cmplx(24.0 * u, 0.0) * f0;
        dcomplex b21 = r11 - cmplx(3.0 * u, 0.0) * r12 - cmplx(24.0 * u, 0.0) * f1;
        dcomplex b22 = r21 - cmplx(3.0 * u, 0.0) * r22 - cmplx(24.0 * u, 0.0) * f2;

        b10 *= cmplx(fden, 0.0);
        b11 *= cmplx(fden, 0.0);
        b12 *= cmplx(fden, 0.0);
        b20 *= cmplx(fden, 0.0);
        b21 *= cmplx(fden, 0.0);
        b22 *= cmplx(fden, 0.0);

        for (int cc = 0; cc < NC * NC; ++cc) {
          dcomplex qt1 = b10 * cmplx(iQ0.r(cc), iQ0.i(cc))
                         + b11 * cmplx(iQ1.i(cc), -iQ1.r(cc))
                         - b12 * cmplx(iQ2.r(cc), iQ2.i(cc));
          B1.set(cc, real(qt1), imag(qt1));

          dcomplex qt2 = b20 * cmplx(iQ0.r(cc), iQ0.i(cc))
                         + b21 * cmplx(iQ1.i(cc), -iQ1.r(cc))
                         - b22 * cmplx(iQ2.r(cc), iQ2.i(cc));
          B2.set(cc, real(qt2), imag(qt2));
        }

        USigmap.mult_nn(U_tmp, Sigmap_tmp);

        tmp1.mult_nn(USigmap, B1);
        tmp2.mult_nn(USigmap, B2);

        dcomplex tr1 = cmplx(tmp1.r(0) + tmp1.r(4) + tmp1.r(8),
                             tmp1.i(0) + tmp1.i(4) + tmp1.i(8));
        dcomplex tr2 = cmplx(tmp2.r(0) + tmp2.r(4) + tmp2.r(8),
                             tmp2.i(0) + tmp2.i(4) + tmp2.i(8));

        iQUS.mult_nn(iQ1, USigmap);
        iUSQ.mult_nn(USigmap, iQ1);

        for (int cc = 0; cc < NC * NC; ++cc) {
          dcomplex qt = tr1 * cmplx(iQ1.i(cc), -iQ1.r(cc))
                        - tr2 * cmplx(iQ2.r(cc), iQ2.i(cc))
                        + f1 * cmplx(USigmap.r(cc), USigmap.i(cc))
                        + f2 * cmplx(iQUS.i(cc), -iQUS.r(cc))
                        + f2 * cmplx(iUSQ.i(cc), -iUSQ.r(cc));
          iGamma.set(cc, -imag(qt), real(qt));
        }
      } else {
        // vout.general(m_vl,"force_recursive: |iQ1|^2 too small: %lf. Set e_iQ=1.\n",norm);
        iGamma.zero();
        e_iQ.unit();
      }

      //! iGamma \f$=i\Lambda\f$
      iGamma.at();
      iTheta_tmp.mult_nn(iGamma, U_tmp);
      //! iTheta \f$=i\Lambda U_\mu(x)\f$
      iTheta.set_mat(site, mu, iTheta_tmp);

      Xi_tmp.mult_nn(Sigmap_tmp, e_iQ);
      Xi_tmp.multadd_dn(C_tmp, iGamma);
      //! Xi \f$=\Sigma_\mu'(x)\exp(iQ_\mu(x))+C_\mu^\dagger i\Lambda_\mu(x)\f$.
      Xi.set_mat(site, mu, Xi_tmp);
    }
  }

  /*
  unsigned long count;
  double time;
  KEK_FopCountFinish(id,&count,&time);
  m_time += time;
  m_flop += count;
  */
  double time1 = Communicator::get_time();
  m_time += time1 - time0;
}


//====================================================================
void Projection_Stout_SU3::set_fj(dcomplex& f0, dcomplex& f1, dcomplex& f2,
                                  const double& u, const double& w)
{
  double xi0   = func_xi0(w);
  double u2    = u * u;
  double w2    = w * w;
  double cos_w = cos(w);

  double cos_u = cos(u);
  double sin_u = sin(u);

  dcomplex emiu = cmplx(cos_u, -sin_u);
  dcomplex e2iu = cmplx(cos_u * cos_u - sin_u * sin_u, 2.0 * sin_u * cos_u);

  dcomplex h0 = e2iu * cmplx(u2 - w2, 0.0)
                + emiu * cmplx(8.0 * u2 * cos_w, 2.0 * u * (3.0 * u2 + w2) * xi0);
  dcomplex h1 = cmplx(2 * u, 0.0) * e2iu
                - emiu * cmplx(2.0 * u * cos_w, -(3.0 * u2 - w2) * xi0);
  dcomplex h2 = e2iu - emiu * cmplx(cos_w, 3.0 * u * xi0);

  double fden = 1.0 / (9.0 * u2 - w2);

  //- output
  f0 = h0 * fden;
  f1 = h1 * fden;
  f2 = h2 * fden;
}


//====================================================================
void Projection_Stout_SU3::set_uw(double& u, double& w,
                                  const Mat_SU_N& iQ2, const Mat_SU_N& iQ3)
{
  double c0    = -(iQ3.i(0, 0) + iQ3.i(1, 1) + iQ3.i(2, 2)) / 3.0;
  double c1    = -0.5 * (iQ2.r(0, 0) + iQ2.r(1, 1) + iQ2.r(2, 2));
  double c13r  = sqrt(c1 / 3.0);
  double c0max = 2.0 * c13r * c13r * c13r;

  double theta = acos(c0 / c0max);

  //- output
  u = c13r * cos(theta / 3.0);
  w = sqrt(c1) * sin(theta / 3.0);
}


//====================================================================
double Projection_Stout_SU3::func_xi0(double w)
{
  if (w == 0.0) {
    return 1.0;
  } else {
    return sin(w) / w;
  }
}


//====================================================================
double Projection_Stout_SU3::func_xi1(double w)
{
  if (w < 0.25) {
    double        w2 = w * w;
    static double c0 = -1.0 / 3.0;
    static double c1 = 1.0 / 30.0;
    static double c2 = -1.0 / 840.0;
    static double c3 = 1.0 / 45360.0;
    static double c4 = -1.0 / 3991680.0;

    return c0 + w2 * (c1 + w2 * (c2 + w2 * (c3 + w2 * c4)));
  } else {
    return (w * cos(w) - sin(w)) / (w * w * w);
  }
}


//====================================================================
void Projection_Stout_SU3::exp_iQ_bf(Field_G& e_iQ, const Field_G& iQ)
{
  // brute force version of exponentiation: for check

  int Nprec = 32;

  int Nvol = iQ.nvol();
  int Nex  = iQ.nex();

  Mat_SU_N u0(NC), u1(NC), u2(NC);
  Mat_SU_N h1(NC);

  for (int ex = 0; ex < Nex; ++ex) {
    for (int site = 0; site < Nvol; ++site) {
      u0.unit();
      u1.unit();
      h1 = iQ.mat(site, ex);

      for (int iprec = 0; iprec < Nprec; ++iprec) {
        double exf = 1.0 / (Nprec - iprec);
        u2  = h1 * u1;
        u2 *= exf;
        u1  = u2;
        u1 += u0;
      }

      u1.reunit();
      e_iQ.set_mat(site, ex, u1);
    }
  }
}


//====================================================================
//============================================================END=====
