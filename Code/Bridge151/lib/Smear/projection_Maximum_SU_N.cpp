/*!
        @file    projection_Maximum_SU_N.cpp

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/
#include "BridgeLib_Private.h"
#include "projection_Maximum_SU_N.h"

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init = Projection_Maximum_SU_N::register_factory();
}
#endif

const std::string Projection_Maximum_SU_N::class_name = "Projection_Maximum_SU_N";

//====================================================================
void Projection_Maximum_SU_N::set_parameters(const Parameters& params)
{
  const string str_vlevel = params.get_string("verbose_level");

  m_vl = vout.set_verbose_level(str_vlevel);

  //- fetch and check input parameters
  int    Niter;
  double Enorm;

  int err = 0;
  err += params.fetch_int("maximum_number_of_iteration", Niter);
  err += params.fetch_double("convergence_criterion", Enorm);

  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameter not found.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }


  set_parameters(Niter, Enorm);
}


//====================================================================
void Projection_Maximum_SU_N::set_parameters(const int Niter, const double Enorm)
{
  //- print input parameters
  vout.general(m_vl, "%s:\n", class_name.c_str());
  vout.general(m_vl, "  Niter  = %d\n", Niter);
  vout.general(m_vl, "  Enorm  = %12.4e\n", Enorm);

  //- range check
  int err = 0;
  err += ParameterCheck::non_negative(Niter);
  err += ParameterCheck::non_zero(Enorm);

  if (err) {
    vout.crucial(m_vl, "Error at %s: parameter range check failed.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }

  //- store values
  m_Niter = Niter;
  m_Enorm = Enorm;
}


//====================================================================
void Projection_Maximum_SU_N::project(Field_G& U,
                                      const double alpha,
                                      const Field_G& Cst, const Field_G& Uorg)
{
  const int Nex  = Uorg.nex();
  const int Nvol = Uorg.nvol();
  const int Nc   = CommonParameters::Nc();

  assert(Cst.nex() == Nex);
  assert(Cst.nvol() == Nvol);
  assert(U.nex() == Nex);
  assert(U.nvol() == Nvol);

  Field_G u_tmp(Nvol, Nex);
  for (int ex = 0; ex < Nex; ++ex) {
    u_tmp.setpart_ex(ex, Cst, ex);
    u_tmp.addpart_ex(ex, Uorg, ex, 1.0 - alpha);
  }

  maxTr(U, u_tmp);
}


//====================================================================
void Projection_Maximum_SU_N::force_recursive(Field_G& Xi, Field_G& iTheta,
                                              const double alpha, const Field_G& Sigmap,
                                              const Field_G& Cst, const Field_G& Uorg)
{
  vout.crucial(m_vl, "Error at %s: force_recursive() is not available.\n", class_name.c_str());
  exit(EXIT_FAILURE);
}


//====================================================================
void Projection_Maximum_SU_N::maxTr(Field_G& G0, const Field_G& Cst)
{
  const int Nc   = CommonParameters::Nc();
  const int Nvol = Cst.nvol();
  const int Nex  = Cst.nex();

  assert(Nvol == G0.nvol());
  assert(Nex == G0.nex());

  const static int Nmt = 1;  // number of subgroup maximization loop:
                             // seems not important because of outer iter-loop.

  Mat_SU_N unity(Nc);
  unity.unit();

  Field_G A(Nvol, Nex);
  A = Cst;

  vout.detailed(m_vl, "Maximum projection start.\n");

  for (int ex = 0; ex < Nex; ++ex) {
    for (int site = 0; site < Nvol; ++site) {
      G0.set_mat(site, ex, unity);
    }
  }

  for (int iter = 0; iter < m_Niter; ++iter) {
    Field_G Udelta(Nvol, Nex);

    for (int ex = 0; ex < Nex; ++ex) {
      for (int site = 0; site < Nvol; ++site) {
        Udelta.set_mat(site, ex, unity);
      }
    }

    for (int imt = 0; imt < Nmt; ++imt) {
      for (int i1 = 0; i1 < Nc; ++i1) {
        int i2 = (i1 + 1) % Nc;
        maxTr_SU2(i1, i2, G0, A, Udelta);
      }
    }
    // for exact check with Fortran version (passed).

    /*
    for(int imt = 0; imt < Nmt; ++imt){
     for(int i = Nc; i > 0; --i){
       int i1 = i % Nc;
       int i2 = (i1 + 1) % Nc;
       maxTr_SU2(i1, i2, G0, A, Udelta);
     }
    }
    */

    //- convergence test
    double retr1 = 0.0;
    for (int ex = 0; ex < Nex; ++ex) {
      for (int site = 0; site < Nvol; ++site) {
        for (int cc = 0; cc < Nc; ++cc) {
          retr1 += Udelta.cmp_r(cc * (1 + Nc), site, ex);
        }
      }
    }
    double retr = Communicator::reduce_sum(retr1);

    const int Npe    = Communicator::size();
    double    deltaV = 1.0 - retr / (Nc * Nvol * Nex * Npe);
    vout.detailed(m_vl, "  iter = %d  deltaV = %12.4e\n", iter, deltaV);

    if (deltaV < m_Enorm) {
      for (int ex = 0; ex < Nex; ++ex) {
        for (int site = 0; site < Nvol; ++site) {
          Mat_SU_N ut(Nc);
          G0.mat_dag(ut, site, ex);
          G0.set_mat(site, ex, ut);
        }
      }

      vout.detailed(m_vl, "Maximum projection converged.\n");

      return;
    }
  }


  vout.crucial(m_vl, "Error at %s: Maximum projection not converged.\n", class_name.c_str());
  exit(EXIT_FAILURE);
}


//====================================================================
void Projection_Maximum_SU_N::maxTr_SU2(const int i1, const int i2, Field_G& Gmax,
                                        Field_G& A, Field_G& Udelta)
{
  const int Nc   = CommonParameters::Nc();
  const int Nvol = A.nvol();
  const int Nex  = A.nex();

  assert(i1 < Nc);
  assert(i2 < Nc);

  const int j1 = mindex(i1, i1, Nc);
  const int j2 = mindex(i2, i2, Nc);
  const int k1 = mindex(i1, i2, Nc);
  const int k2 = mindex(i2, i1, Nc);

  //----------[     | # # 0 | <i1  ]--------------------------
  //----------[ V = | # # 0 | <i2  ]--------------------------
  //----------[     | 0 0 1 |      ]--------------------------

  Field_G v(Nvol, Nex);
  for (int ex = 0; ex < Nex; ++ex) {
    for (int site = 0; site < Nvol; ++site) {
      Mat_SU_N at(Nc);
      at = A.mat(site, ex);

      double xlamd =
        at.r(j1) * at.r(j1) + at.i(j1) * at.i(j1) + 2.0 * at.r(j1) * at.r(j2)
        + at.r(k1) * at.r(k1) + at.i(k1) * at.i(k1) - 2.0 * at.i(j1) * at.i(j2)
        + at.r(k2) * at.r(k2) + at.i(k2) * at.i(k2) - 2.0 * at.r(k1) * at.r(k2)
        + at.r(j2) * at.r(j2) + at.i(j2) * at.i(j2) + 2.0 * at.i(k1) * at.i(k2);
      xlamd = 1.0 / sqrt(xlamd);

      Mat_SU_N vt(Nc);
      vt.unit();
      vt.set(j1, (at.r(j1) + at.r(j2)) * xlamd, (-at.i(j1) + at.i(j2)) * xlamd);
      vt.set(k1, (at.r(k2) - at.r(k1)) * xlamd, (-at.i(k2) - at.i(k1)) * xlamd);
      vt.set(k2, (at.r(k1) - at.r(k2)) * xlamd, (-at.i(k1) - at.i(k2)) * xlamd);
      vt.set(j2, (at.r(j1) + at.r(j2)) * xlamd, (at.i(j1) - at.i(j2)) * xlamd);

      v.set_mat(site, ex, vt);
    }
  }

  Field_G w(Nvol, Nex);
  for (int ex = 0; ex < Nex; ++ex) {
    mult_Field_Gnn(w, ex, A, ex, v, ex);
  }
  A = w;

  for (int ex = 0; ex < Nex; ++ex) {
    mult_Field_Gnn(w, ex, Gmax, ex, v, ex);
  }
  Gmax = w;

  for (int ex = 0; ex < Nex; ++ex) {
    mult_Field_Gnn(w, ex, Udelta, ex, v, ex);
  }
  Udelta = w;
}


//====================================================================
//============================================================END=====
