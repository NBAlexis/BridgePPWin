/*!
        @file    $Id:: projection_Maximum_SU_N.cpp #$

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 00:49:46 #$

        @version $LastChangedRevision: 1561 $
*/

#include "projection_Maximum_SU_N.h"


#ifdef USE_FACTORY
namespace {
  Projection *create_object()
  {
    return new Projection_Maximum_SU_N();
  }


  bool init = Projection::Factory::Register("Maximum_SU_N", create_object);
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
                                      double alpha,
                                      const Field_G& Cst, const Field_G& Uorg)
{
  int Nex  = Uorg.nex();
  int Nvol = Uorg.nvol();
  int Nc   = CommonParameters::Nc();

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
                                              double alpha, const Field_G& Sigmap,
                                              const Field_G& Cst, const Field_G& Uorg)
{
  vout.crucial(m_vl, "Error at %s: force_recursive() is not available.\n", class_name.c_str());
  exit(EXIT_FAILURE);
}


//====================================================================
void Projection_Maximum_SU_N::maxTr(Field_G& G0, Field_G& Cst)
{
  int Nc   = CommonParameters::Nc();
  int Nvol = Cst.nvol();
  int Nex  = Cst.nex();

  assert(Nvol == G0.nvol());
  assert(Nex == G0.nex());

  int Nmt = 1;  // number of subgroup maximization loop:
                // seems not important because of outer iter-loop.

  Mat_SU_N unity(Nc);
  unity.unit();

  Field_G Udelta(Nvol, Nex), A(Nvol, Nex);
  A = Cst;

  vout.detailed(m_vl, "Maximum projection start.\n");

  for (int ex = 0; ex < Nex; ++ex) {
    for (int site = 0; site < Nvol; ++site) {
      G0.set_mat(site, ex, unity);
    }
  }

  for (int iter = 0; iter < m_Niter; ++iter) {
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

    double retr   = Communicator::reduce_sum(retr1);
    int    Npe    = Communicator::size();
    double deltaV = 1.0 - retr / (Nc * Nvol * Nex * Npe);
    vout.detailed(m_vl, "  iter = %d  deltaV = %12.4e\n", iter, deltaV);

    if (deltaV < m_Enorm) {
      Mat_SU_N ut(Nc);
      for (int ex = 0; ex < Nex; ++ex) {
        for (int site = 0; site < Nvol; ++site) {
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
void Projection_Maximum_SU_N::maxTr_SU2(int i1, int i2, Field_G& Gmax,
                                        Field_G& A, Field_G& Udelta)
{
  int Nc   = CommonParameters::Nc();
  int Nvol = A.nvol();
  int Nex  = A.nex();

  assert(i1 < Nc);
  assert(i2 < Nc);

  int j1 = mindex(i1, i1, Nc);
  int j2 = mindex(i2, i2, Nc);
  int k1 = mindex(i1, i2, Nc);
  int k2 = mindex(i2, i1, Nc);

  Mat_SU_N at(Nc), vt(Nc);
  Field_G  v(Nvol, Nex), w(Nvol, Nex);

  //----------[     | # # 0 | <i1  ]--------------------------
  //----------[ V = | # # 0 | <i2  ]--------------------------
  //----------[     | 0 0 1 |      ]--------------------------

  for (int ex = 0; ex < Nex; ++ex) {
    for (int site = 0; site < Nvol; ++site) {
      at = A.mat(site, ex);

      double xlamd =
        at.r(j1) * at.r(j1) + at.i(j1) * at.i(j1) + 2.0 * at.r(j1) * at.r(j2)
        + at.r(k1) * at.r(k1) + at.i(k1) * at.i(k1) - 2.0 * at.i(j1) * at.i(j2)
        + at.r(k2) * at.r(k2) + at.i(k2) * at.i(k2) - 2.0 * at.r(k1) * at.r(k2)
        + at.r(j2) * at.r(j2) + at.i(j2) * at.i(j2) + 2.0 * at.i(k1) * at.i(k2);
      xlamd = 1.0 / sqrt(xlamd);

      vt.unit();
      vt.set(j1, (at.r(j1) + at.r(j2)) * xlamd, (-at.i(j1) + at.i(j2)) * xlamd);
      vt.set(k1, (at.r(k2) - at.r(k1)) * xlamd, (-at.i(k2) - at.i(k1)) * xlamd);
      vt.set(k2, (at.r(k1) - at.r(k2)) * xlamd, (-at.i(k1) - at.i(k2)) * xlamd);
      vt.set(j2, (at.r(j1) + at.r(j2)) * xlamd, (at.i(j1) - at.i(j2)) * xlamd);

      v.set_mat(site, ex, vt);
    }
  }

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
