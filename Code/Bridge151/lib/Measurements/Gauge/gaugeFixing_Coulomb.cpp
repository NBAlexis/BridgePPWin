/*!
        @file    gaugeFixing_Coulomb.cpp

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/
#include "BridgeLib_Private.h"
#include "gaugeFixing_Coulomb.h"

#include "staple_lex.h"

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init = GaugeFixing_Coulomb::register_factory();
}
#endif

const std::string GaugeFixing_Coulomb::class_name = "GaugeFixing_Coulomb";

//====================================================================
void GaugeFixing_Coulomb::set_parameters(const Parameters& params)
{
  const string str_vlevel = params.get_string("verbose_level");

  m_vl = vout.set_verbose_level(str_vlevel);

  //- fetch and check input parameters
  int    Niter, Nnaive, Nmeas, Nreset;
  double Enorm, wp;

  int err = 0;
  err += params.fetch_int("maximum_number_of_iteration", Niter);
  err += params.fetch_int("number_of_naive_iteration", Nnaive);
  err += params.fetch_int("interval_of_measurement", Nmeas);
  err += params.fetch_int("iteration_to_reset", Nreset);
  err += params.fetch_double("convergence_criterion_squared", Enorm);
  err += params.fetch_double("overrelaxation_parameter", wp);

  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameter not found.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }


  set_parameters(Niter, Nnaive, Nmeas, Nreset, Enorm, wp);
}


//====================================================================
void GaugeFixing_Coulomb::set_parameters(const int Niter, const int Nnaive,
                                         const int Nmeas, const int Nreset,
                                         const double Enorm, const double wp)
{
  //- print input parameters
  vout.general(m_vl, "%s:\n", class_name.c_str());
  vout.general(m_vl, "  Niter  = %d\n", Niter);
  vout.general(m_vl, "  Nnaive = %d\n", Nnaive);
  vout.general(m_vl, "  Nmeas  = %d\n", Nmeas);
  vout.general(m_vl, "  Nreset = %d\n", Nreset);
  vout.general(m_vl, "  Enorm  = %12.4e\n", Enorm);
  vout.general(m_vl, "  wp     = %8.4f\n", wp);

  //- range check
  int err = 0;
  err += ParameterCheck::non_negative(Niter);
  err += ParameterCheck::non_negative(Nnaive);
  err += ParameterCheck::non_negative(Nmeas);
  err += ParameterCheck::non_negative(Nreset);
  err += ParameterCheck::square_non_zero(Enorm);
  err += ParameterCheck::non_zero(wp);

  if (err) {
    vout.crucial(m_vl, "Error at %s: parameter range check failed.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }

  //- store values
  m_Niter  = Niter;
  m_Nnaive = Nnaive;
  m_Nmeas  = Nmeas;
  m_Nreset = Nreset;
  m_Enorm  = Enorm;
  m_wp     = wp;
}


//====================================================================
void GaugeFixing_Coulomb::fix(Field_G& Ufix, const Field_G& Uorg)
{
  const int Nvol = Uorg.nvol();
  const int Nex  = Uorg.nex();
  const int Lt   = CommonParameters::Lt();

  const int Nvol2 = Nvol / 2;

  Field_G Ue(Nvol2, Nex);

  m_index.convertField(Ue, Uorg, 0);

  Field_G Uo(Nvol2, Nex);
  m_index.convertField(Uo, Uorg, 1);

  int Nconv = -1;

  Staple_lex   staple;
  const double plaq = staple.plaquette(Uorg);

  vout.general(m_vl, "plaq(original) = %18.14f\n", plaq);


  //- gauge fixing iteration
  for (int iter = 0; iter < m_Niter; ++iter) {
    std::valarray<double> sg(Lt);
    std::valarray<double> Fval(Lt);

    if ((iter % m_Nmeas) == 0) {
      calc_SG(sg, Fval, Ue, Uo);

      double sg_max = 0.0;
      for (int t = 0; t < Lt; ++t) {
        if (sg[t] > sg_max) sg_max = sg[t];
      }

      vout.detailed(m_vl, "  iter = %6d: sg_max = %16.8e\n", iter, sg_max);

      for (int t = 0; t < Lt; ++t) {
        vout.paranoiac(m_vl, "    t = %4d  sg = %16.8e  Fval = %16.8e\n",
                       t, sg[t], Fval[t]);
      }

      if (sg_max < m_Enorm) {
        Nconv = iter;
        vout.general(m_vl, "converged at iter = %d\n", Nconv);
        break;
      }
    }

    double wp2 = m_wp;
    if ((iter % m_Nreset) < m_Nnaive) wp2 = 1.0;
    gfix_step(Ue, Uo, wp2);

    if (((iter % m_Nreset) == 0) && (iter > 0)) {
      vout.detailed(m_vl, "  random gauge transformation performed.\n");

      Field_G Ge(Nvol2, 1);
      set_randomGaugeTrans(sg, Ge);
      gauge_trans_eo(Ue, Uo, Ge, 0);

      Field_G Go(Nvol2, 1);
      set_randomGaugeTrans(sg, Go);
      gauge_trans_eo(Ue, Uo, Go, 1);
    }
  }

  if (Nconv < 0) {
    vout.crucial(m_vl, "Error at %s: not converged.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }


  m_index.reverseField(Ufix, Ue, 0);
  m_index.reverseField(Ufix, Uo, 1);


  const double plaq2     = staple.plaquette(Ufix);
  const double plaq_diff = std::abs(plaq - plaq2);

  vout.general(m_vl, "plaq(fixed)    = %18.14f\n", plaq2);
  vout.general(m_vl, "plaq(diff)     = %18.10e\n", plaq_diff);

  if (plaq_diff > sqrt(m_Enorm)) {
    vout.crucial(m_vl, "Error at %s: too large plaq(diff) = %20.14e\n",
                 class_name.c_str(), plaq_diff);
    exit(EXIT_FAILURE);
  }
}


//====================================================================
void GaugeFixing_Coulomb::gfix_step(Field_G& Ue, Field_G& Uo, const double wp)
{
  const int Nc    = CommonParameters::Nc();
  const int Nvol2 = Ue.nvol();

  Mat_SU_N uwp(Nc);

  uwp.unit();
  uwp *= (1.0 - wp);

  for (int ieo = 0; ieo < 2; ++ieo) {
    Field_G Weo(Nvol2, 1);
    calc_W(Weo, Ue, Uo, ieo);

    Field_G Geo(Nvol2, 1);
    maxTr(Geo, Weo);

    scal(Geo, wp);

    for (int site = 0; site < Nvol2; ++site) {
      Mat_SU_N ut(Nc);
      ut  = Geo.mat(site, 0);
      ut += uwp;
      ut.reunit();
      Geo.set_mat(site, 0, ut);
    }

    gauge_trans_eo(Ue, Uo, Geo, ieo);
  }
}


//====================================================================
void GaugeFixing_Coulomb::set_randomGaugeTrans(const std::valarray<double>& sg,
                                               Field_G& Geo)
{
  const int Lt    = CommonParameters::Lt();
  const int Nt    = CommonParameters::Nt();
  const int Nz    = CommonParameters::Nz();
  const int Ny    = CommonParameters::Ny();
  const int Nx2   = CommonParameters::Nx() / 2;
  const int Nc    = CommonParameters::Nc();
  const int ipe_t = Communicator::ipe(3);

  assert(Geo.nex() == 1);
  assert(sg.size() == Lt);

  for (int t = 0; t < Nt; ++t) {
    int tg = t + ipe_t * Nt;

    if (sg[tg] > m_Enorm) {
      for (int z = 0; z < Nz; ++z) {
        for (int y = 0; y < Ny; ++y) {
          for (int x = 0; x < Nx2; ++x) {
            int site = m_index.siteh(x, y, z, t);

            Mat_SU_N gt(Nc);
            gt.set_random(m_rand);
            Geo.set_mat(site, 0, gt);
          }
        }
      }
    } else {
      for (int z = 0; z < Nz; ++z) {
        for (int y = 0; y < Ny; ++y) {
          for (int x = 0; x < Nx2; ++x) {
            int site = m_index.siteh(x, y, z, t);

            Mat_SU_N gt(Nc);
            gt.unit();
            Geo.set_mat(site, 0, gt);
          }
        }
      }
    }
  }
}


//====================================================================
void GaugeFixing_Coulomb::gauge_trans_eo(Field_G& Ue, Field_G& Uo,
                                         const Field_G& Geo, const int Ieo)
{
  //  Ieo = 0: gauge transformation on even sites.
  //  Ieo = 1:                      on odd sites.

  const int Nvol2 = Geo.nvol();
  const int Ndim  = CommonParameters::Ndim();

  ShiftField_eo shift;

  if (Ieo == 0) {
    for (int mu = 0; mu < Ndim; ++mu) {
      Field_G Ut(Nvol2, 1);
      mult_Field_Gnn(Ut, 0, Geo, 0, Ue, mu);
      Ue.setpart_ex(mu, Ut, 0);

      Field_G Gt(Nvol2, 1);
      shift.backward_h(Gt, Geo, mu, 1);
      mult_Field_Gnd(Ut, 0, Uo, mu, Gt, 0);
      Uo.setpart_ex(mu, Ut, 0);
    }
  } else {
    for (int mu = 0; mu < Ndim; ++mu) {
      Field_G Ut(Nvol2, 1);
      mult_Field_Gnn(Ut, 0, Geo, 0, Uo, mu);
      Uo.setpart_ex(mu, Ut, 0);

      Field_G Gt(Nvol2, 1);
      shift.backward_h(Gt, Geo, mu, 0);
      mult_Field_Gnd(Ut, 0, Ue, mu, Gt, 0);
      Ue.setpart_ex(mu, Ut, 0);
    }
  }
}


//====================================================================
void GaugeFixing_Coulomb::calc_SG(std::valarray<double>& sg,
                                  std::valarray<double>& Fval,
                                  const Field_G& Ue, const Field_G& Uo)
{
  const int Nc    = CommonParameters::Nc();
  const int Lt    = CommonParameters::Lt();
  const int Nt    = CommonParameters::Nt();
  const int Nz    = CommonParameters::Nz();
  const int Ny    = CommonParameters::Ny();
  const int Nx2   = CommonParameters::Nx() / 2;
  const int Nvol2 = Nx2 * Ny * Nz * Nt;
  const int Ndim  = CommonParameters::Ndim();
  const int NPE   = CommonParameters::NPE();

  assert(Ue.nex() == Ndim);
  assert(Ue.nvol() == Nvol2);
  assert(Uo.nex() == Ndim);
  assert(Uo.nvol() == Nvol2);

  sg   = 0.0;
  Fval = 0.0;

  std::valarray<double> sg_local(Nt);
  sg_local = 0.0;

  for (int ieo = 0; ieo < 2; ++ieo) {
    Field_G DLT(Nvol2, 1);
    calc_DLT(DLT, Ue, Uo, ieo);

    for (int t = 0; t < Nt; ++t) {
      for (int z = 0; z < Nz; ++z) {
        for (int y = 0; y < Ny; ++y) {
          for (int x = 0; x < Nx2; ++x) {
            int site = m_index.siteh(x, y, z, t);

            Mat_SU_N ut(Nc);
            ut           = DLT.mat(site, 0);
            sg_local[t] += ut.norm2();
          }
        }
      }
    }
  }

  sum_global_t(sg, sg_local);
  for (int t = 0; t < Lt; ++t) {
    sg[t] = sg[t] / (Ndim * Nc * (2 * Nvol2 * NPE) / Lt);
  }


  std::valarray<double> Fval_local(Nt);
  Fval_local = 0.0;

  for (int mu = 0; mu < Ndim - 1; ++mu) {
    for (int t = 0; t < Nt; ++t) {
      for (int z = 0; z < Nz; ++z) {
        for (int y = 0; y < Ny; ++y) {
          for (int x = 0; x < Nx2; ++x) {
            int site = m_index.siteh(x, y, z, t);

            Mat_SU_N ut(Nc);
            ut             = Ue.mat(site, mu);
            Fval_local[t] += ReTr(ut);
            ut             = Uo.mat(site, mu);
            Fval_local[t] += ReTr(ut);
          }
        }
      }
    }
  }

  sum_global_t(Fval, Fval_local);
  for (int t = 0; t < Lt; ++t) {
    Fval[t] = Fval[t] / (Ndim * (2 * Nvol2 * NPE) / Lt);
  }
}


//====================================================================
void GaugeFixing_Coulomb::sum_global_t(std::valarray<double>& val_global,
                                       const std::valarray<double>& val_local)
{
  const int Lt    = CommonParameters::Lt();
  const int Nt    = CommonParameters::Nt();
  const int ipe_t = Communicator::ipe(3);

  assert(val_global.size() == Lt);
  assert(val_local.size() == Nt);

  for (int t_global = 0; t_global < Lt; ++t_global) {
    val_global[t_global] = 0.0;
  }

  for (int t = 0; t < Nt; ++t) {
    int t_global = t + ipe_t * Nt;
    val_global[t_global] = val_local[t];
  }

  for (int t_global = 0; t_global < Lt; ++t_global) {
    double val = val_global[t_global];
    val_global[t_global] = Communicator::reduce_sum(val);
  }
}


//====================================================================
void GaugeFixing_Coulomb::calc_DLT(Field_G& DLT,
                                   const Field_G& Ue, const Field_G& Uo,
                                   const int Ieo)
{
  const int Nvol2 = Ue.nvol();
  const int Nc    = CommonParameters::Nc();
  const int Ndim  = CommonParameters::Ndim();

  ShiftField_eo shift;

  DLT.set(0.0);

  if (Ieo == 0) { // on even sites
    for (int mu = 0; mu < Ndim - 1; ++mu) {
      DLT.addpart_ex(0, Ue, mu, -1.0);

      Field_G Ut1(Nvol2, 1);
      Ut1.setpart_ex(0, Uo, mu);

      Field_G Ut2(Nvol2, 1);
      shift.forward_h(Ut2, Ut1, mu, 0);
      DLT.addpart_ex(0, Ut2, 0);
    }
  } else {        // on odd sites
    for (int mu = 0; mu < Ndim - 1; ++mu) {
      DLT.addpart_ex(0, Uo, mu, -1.0);

      Field_G Ut1(Nvol2, 1);
      Ut1.setpart_ex(0, Ue, mu);

      Field_G Ut2(Nvol2, 1);
      shift.forward_h(Ut2, Ut1, mu, 1);
      DLT.addpart_ex(0, Ut2, 0);
    }
  }

  for (int site = 0; site < Nvol2; ++site) {
    Mat_SU_N u_tmp(Nc);

    u_tmp = DLT.mat(site, 0);
    u_tmp.at();
    u_tmp *= 2.0;
    DLT.set_mat(site, 0, u_tmp);
  }
}


//====================================================================
void GaugeFixing_Coulomb::calc_W(Field_G& Weo,
                                 const Field_G& Ue, const Field_G& Uo,
                                 const int Ieo)
{
  const int Nvol2 = Ue.nvol();
  const int Nc    = CommonParameters::Nc();
  const int Ndim  = CommonParameters::Ndim();

  assert(Weo.nex() == 1);

  ShiftField_eo shift;

  Weo.set(0.0);

  if (Ieo == 0) {       // on even sites
    for (int mu = 0; mu < Ndim - 1; ++mu) {
      Weo.addpart_ex(0, Ue, mu);

      Field_G Ut1(Nvol2, 1);
      Ut1.setpart_ex(0, Uo, mu);

      Field_G Ut2(Nvol2, 1);
      shift.forward_h(Ut2, Ut1, mu, 0);
      for (int site = 0; site < Nvol2; ++site) {
        Mat_SU_N u_tmp(Nc);
        u_tmp = Ut2.mat_dag(site, 0);
        Weo.add_mat(site, 0, u_tmp);
      }
    }
  } else if (Ieo == 1) { // on odd sites
    for (int mu = 0; mu < Ndim - 1; ++mu) {
      Weo.addpart_ex(0, Uo, mu);

      Field_G Ut1(Nvol2, 1);
      Ut1.setpart_ex(0, Ue, mu);

      Field_G Ut2(Nvol2, 1);
      shift.forward_h(Ut2, Ut1, mu, 1);
      for (int site = 0; site < Nvol2; ++site) {
        Mat_SU_N u_tmp(Nc);
        u_tmp = Ut2.mat_dag(site, 0);
        Weo.add_mat(site, 0, u_tmp);
      }
    }
  } else {
    vout.crucial(m_vl, "Error at %s: Wrong ieo=%d\n", class_name.c_str(), Ieo);
    exit(EXIT_FAILURE);
  }
}


//====================================================================
void GaugeFixing_Coulomb::maxTr(Field_G& G0, Field_G& W)
{
  // Present implementation only applys to SU(3) case.
  const int Nc    = CommonParameters::Nc();
  const int Nvol2 = G0.nvol();

  const int Nmt = 1;

  Mat_SU_N unity(Nc);

  unity.unit();

  for (int site = 0; site < Nvol2; ++site) {
    G0.set_mat(site, 0, unity);
  }

  for (int imt = 0; imt < Nmt; ++imt) {
    maxTr1(G0, W);
    maxTr2(G0, W);
    maxTr3(G0, W);
  }
}


//====================================================================
void GaugeFixing_Coulomb::maxTr1(Field_G& G, Field_G& W)
{
  const int Nc    = CommonParameters::Nc();
  const int Nvol2 = W.nvol();

  for (int site = 0; site < Nvol2; ++site) {
    Mat_SU_N wt(Nc);
    wt = W.mat(site, 0);

    Mat_SU_N gt(Nc);
    gt.set(2, 0.0, 0.0);
    gt.set(5, 0.0, 0.0);
    gt.set(6, 0.0, 0.0);
    gt.set(7, 0.0, 0.0);
    gt.set(8, 1.0, 0.0);

    double fn1 = (wt.r(0) + wt.r(4)) * (wt.r(0) + wt.r(4))
                 + (wt.i(0) - wt.i(4)) * (wt.i(0) - wt.i(4));
    double fn2 = (wt.r(1) - wt.r(3)) * (wt.r(1) - wt.r(3))
                 + (wt.i(1) + wt.i(3)) * (wt.i(1) + wt.i(3));
    double fn = 1.0 / sqrt(fn1 + fn2);

    gt.set(0, fn * (wt.r(0) + wt.r(4)), fn * (-wt.i(0) + wt.i(4)));
    gt.set(1, fn * (-wt.r(1) + wt.r(3)), fn * (-wt.i(1) - wt.i(3)));
    gt.set(3, fn * (wt.r(1) - wt.r(3)), fn * (-wt.i(1) - wt.i(3)));
    gt.set(4, fn * (wt.r(0) + wt.r(4)), fn * (wt.i(0) - wt.i(4)));

    Mat_SU_N wt2(Nc);
    wt2 = gt * wt;
    W.set_mat(site, 0, wt2);

    Mat_SU_N gt2(Nc);
    gt2 = G.mat(site, 0);
    wt2 = gt * gt2;
    G.set_mat(site, 0, wt2);
  }
}


//====================================================================
void GaugeFixing_Coulomb::maxTr2(Field_G& G, Field_G& W)
{
  const int Nc    = CommonParameters::Nc();
  const int Nvol2 = W.nvol();

  Mat_SU_N gt2(Nc), wt2(Nc);

  for (int site = 0; site < Nvol2; ++site) {
    Mat_SU_N wt(Nc);
    wt = W.mat(site, 0);

    Mat_SU_N gt(Nc);
    gt.set(1, 0.0, 0.0);
    gt.set(3, 0.0, 0.0);
    gt.set(4, 1.0, 0.0);
    gt.set(5, 0.0, 0.0);
    gt.set(7, 0.0, 0.0);

    double fn1 = (wt.r(8) + wt.r(0)) * (wt.r(8) + wt.r(0))
                 + (wt.i(8) - wt.i(0)) * (wt.i(8) - wt.i(0));
    double fn2 = (wt.r(2) - wt.r(6)) * (wt.r(2) - wt.r(6))
                 + (wt.i(2) + wt.i(6)) * (wt.i(2) + wt.i(6));
    double fn = 1.0 / sqrt(fn1 + fn2);

    gt.set(0, fn * (wt.r(8) + wt.r(0)), fn * (wt.i(8) - wt.i(0)));
    gt.set(2, fn * (wt.r(6) - wt.r(2)), fn * (-wt.i(6) - wt.i(2)));
    gt.set(6, fn * (-wt.r(6) + wt.r(2)), fn * (-wt.i(6) - wt.i(2)));
    gt.set(8, fn * (wt.r(8) + wt.r(0)), fn * (-wt.i(8) + wt.i(0)));

    Mat_SU_N wt2_l(Nc);
    wt2_l = gt * wt;
    W.set_mat(site, 0, wt2_l);

    Mat_SU_N gt2_l(Nc);
    gt2_l = G.mat(site, 0);
    wt2_l = gt * gt2_l;
    G.set_mat(site, 0, wt2_l);
  }
}


//====================================================================
void GaugeFixing_Coulomb::maxTr3(Field_G& G, Field_G& W)
{
  const int Nc    = CommonParameters::Nc();
  const int Nvol2 = W.nvol();

  Mat_SU_N gt2(Nc), wt2(Nc);

  for (int site = 0; site < Nvol2; ++site) {
    Mat_SU_N wt(Nc);
    wt = W.mat(site, 0);

    Mat_SU_N gt(Nc);
    gt.set(0, 1.0, 0.0);
    gt.set(1, 0.0, 0.0);
    gt.set(2, 0.0, 0.0);
    gt.set(3, 0.0, 0.0);
    gt.set(6, 0.0, 0.0);

    double fn1 = (wt.r(4) + wt.r(8)) * (wt.r(4) + wt.r(8))
                 + (wt.i(4) - wt.i(8)) * (wt.i(4) - wt.i(8));
    double fn2 = (wt.r(7) - wt.r(5)) * (wt.r(7) - wt.r(5))
                 + (wt.i(7) + wt.i(5)) * (wt.i(7) + wt.i(5));
    double fn = 1.0 / sqrt(fn1 + fn2);

    gt.set(4, fn * (wt.r(4) + wt.r(8)), fn * (-wt.i(4) + wt.i(8)));
    gt.set(5, fn * (-wt.r(5) + wt.r(7)), fn * (-wt.i(5) - wt.i(7)));
    gt.set(7, fn * (wt.r(5) - wt.r(7)), fn * (-wt.i(5) - wt.i(7)));
    gt.set(8, fn * (wt.r(4) + wt.r(8)), fn * (wt.i(4) - wt.i(8)));

    Mat_SU_N wt2_l(Nc);
    wt2_l = gt * wt;
    W.set_mat(site, 0, wt2_l);

    Mat_SU_N gt2_l(Nc);
    gt2_l = G.mat(site, 0);
    wt2_l = gt * gt2_l;
    G.set_mat(site, 0, wt2_l);
  }
}


//====================================================================
//============================================================END=====
