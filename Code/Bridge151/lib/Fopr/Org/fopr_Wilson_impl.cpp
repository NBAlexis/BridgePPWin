/*!
        @file    fopr_Wilson_impl.cpp

        @brief

        @author  Hideo Matsufuru  (matsufuru)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/
#include "BridgeLib_Private.h"
#if !USE_IMP
#include "fopr_Wilson_impl.h"

namespace Org {
#ifdef USE_FACTORY_AUTOREGISTER
  namespace {
    bool init = Fopr_Wilson::register_factory();
  }
#endif

  const std::string Fopr_Wilson::class_name = "Org::Fopr_Wilson";

//====================================================================
  void Fopr_Wilson::init(const std::string repr)
  {
    m_vl = CommonParameters::Vlevel();

    m_Nc = CommonParameters::Nc();
    m_Nd = CommonParameters::Nd();

    m_Nvol = CommonParameters::Nvol();
    m_Ndim = CommonParameters::Ndim();
    m_boundary.resize(m_Ndim);

    m_U = 0;

    m_repr = repr;

    m_GM.resize(m_Ndim + 1);

    unique_ptr<GammaMatrixSet> gmset(GammaMatrixSet::New(m_repr));

    m_GM[0] = gmset->get_GM(GammaMatrixSet::GAMMA1);
    m_GM[1] = gmset->get_GM(GammaMatrixSet::GAMMA2);
    m_GM[2] = gmset->get_GM(GammaMatrixSet::GAMMA3);
    m_GM[3] = gmset->get_GM(GammaMatrixSet::GAMMA4);
    m_GM[4] = gmset->get_GM(GammaMatrixSet::GAMMA5);

    m_mult     = &Fopr_Wilson::mult_undef;
    m_mult_dag = &Fopr_Wilson::mult_undef;
  }


//====================================================================
  void Fopr_Wilson::set_mode(const std::string mode)
  {
    m_mode = mode;

    if (m_mode == "D") {
      m_mult     = &Fopr_Wilson::D;
      m_mult_dag = &Fopr_Wilson::Ddag;
    } else if (m_mode == "Ddag") {
      m_mult     = &Fopr_Wilson::Ddag;
      m_mult_dag = &Fopr_Wilson::D;
    } else if (m_mode == "DdagD") {
      m_mult     = &Fopr_Wilson::DdagD;
      m_mult_dag = &Fopr_Wilson::DdagD;
    } else if (m_mode == "DDdag") {
      m_mult     = &Fopr_Wilson::DDdag;
      m_mult_dag = &Fopr_Wilson::DDdag;
    } else if (m_mode == "H") {
      m_mult     = &Fopr_Wilson::H;
      m_mult_dag = &Fopr_Wilson::H;
    } else {
      vout.crucial(m_vl, "Error at %s: input mode is undefined.\n", class_name.c_str());
      exit(EXIT_FAILURE);
    }
  }


//====================================================================
  std::string Fopr_Wilson::get_mode() const
  {
    return m_mode;
  }


//====================================================================
  void Fopr_Wilson::set_parameters(const Parameters& params)
  {
    const std::string str_vlevel = params.get_string("verbose_level");

    m_vl = vout.set_verbose_level(str_vlevel);

    //- fetch and check input parameters
    double           kappa;
    std::vector<int> bc;

    int err = 0;
    err += params.fetch_double("hopping_parameter", kappa);
    err += params.fetch_int_vector("boundary_condition", bc);

    if (err) {
      vout.crucial(m_vl, "Error at %s: input parameter not found.\n", class_name.c_str());
      exit(EXIT_FAILURE);
    }

    set_parameters(kappa, bc);
  }


//====================================================================
  void Fopr_Wilson::set_parameters(const double kappa, const std::vector<int> bc)
  {
    //- print input parameters
    vout.general(m_vl, "%s:\n", class_name.c_str());
    vout.general(m_vl, "  kappa = %12.8f\n", kappa);
    for (int mu = 0; mu < m_Ndim; ++mu) {
      vout.general(m_vl, "  boundary[%d] = %2d\n", mu, bc[mu]);
    }

    //- range check
    // NB. kappa = 0 is allowed.
    assert(bc.size() == m_Ndim);

    //- store values
    m_kappa = kappa;

    // m_boundary.resize(m_Ndim);  // NB. already resized in init.
    m_boundary = bc;
  }


//====================================================================
  void Fopr_Wilson::D(Field& v, const Field& f)
  {
    v = f;

    Field_F w(f.nvol(), f.nex());
    w.set(0.0);

    for (int mu = 0; mu < m_Ndim; ++mu) {
      mult_up(mu, w, f);
      mult_dn(mu, w, f);
    }

    axpy(v, -m_kappa, w);
  }


//====================================================================
  void Fopr_Wilson::D_ex(Field& v, const int ex1, const Field& f, const int ex2)
  {
    Field ff(f.nin(), f.nvol(), 1);

    ff.setpart_ex(0, f, ex2);

    Field w(f.nin(), f.nvol(), 1);

    for (int mu = 0; mu < m_Ndim; ++mu) {
      mult_up(mu, w, ff);
      mult_dn(mu, w, ff);
    }

    scal(w, -m_kappa);

    v.addpart_ex(ex1, w, 0);
  }


//====================================================================
  void Fopr_Wilson::mult_gm5(Field& v, const Field& f)
  {
    assert(v.nvol() == f.nvol());
    assert(v.nex() == f.nex());
    assert(v.nin() == f.nin());

    Field_F vt(f.nvol(), f.nex());
    mult_GM(vt, m_GM[4], (Field_F)f);
    v = (Field)vt;
  }


//====================================================================
  void Fopr_Wilson::proj_chiral(Field& w,
                                const int ex1, const Field& v, const int ex2, const int ipm)
  {
    assert(ipm == 1 || ipm == -1);

    Field vv(v.nin(), v.nvol(), 1);
    vv.setpart_ex(0, v, ex2);

    Field ww(v.nin(), v.nvol(), 1);
    mult_gm5(ww, vv);

    if (ipm == 1) {
    } else if (ipm == -1) {
      scal(ww, -1.0);
    }

    w.addpart_ex(ex1, ww, 0);
  }


//====================================================================

/*
const Field_F Fopr_Wilson::mult_gm5p(int mu, const Field_F& w)
{
  Field_F vt, v;

  vt.set(0.0);

  assert(mu >= 0);
  assert(mu < m_Ndim);

  mult_up(mu, vt, w);

  mult_gm5(v, vt);

  return v;
}
*/

//====================================================================
  void Fopr_Wilson::mult_gm5p(const int mu, Field_F& v, const Field_F& w)
  {
    assert(mu >= 0);
    assert(mu < m_Ndim);

    Field_F vt;
    mult_up(mu, vt, w);
    mult_gm5(v, vt);
  }


//====================================================================
  void Fopr_Wilson::mult_up(const int mu, Field& w, const Field& f)
  {
    for (int ex = 0; ex < f.nex(); ++ex) {
      Field_F vt(f.nvol(), 1);
      vt.setpart_ex(0, f, ex);

      m_shift.backward(m_trf, f, m_boundary[mu], mu);

      mult_Field_Gn(m_trf2, 0, *m_U, mu, m_trf, 0);
      mult_GMproj2(vt, -1, m_GM[mu], m_trf2);

      w.addpart_ex(ex, vt, 0);
    }
  }


//====================================================================
  void Fopr_Wilson::mult_dn(const int mu, Field& w, const Field& f)
  {
    for (int ex = 0; ex < f.nex(); ++ex) {
      mult_Field_Gd(m_trf, 0, *m_U, mu, (Field_F)f, ex);
      m_shift.forward(m_trf2, m_trf, m_boundary[mu], mu);

      Field_F vt(f.nvol(), 1);
      mult_GMproj2(vt, 1, m_GM[mu], m_trf2);

      w.addpart_ex(ex, vt, 0);
    }
  }


//====================================================================
  double Fopr_Wilson::flop_count()
  {
    // This counting is based on the Org-implementation of ver.1.2.0.
    // Flop count of mult_GMproj2 is different for upward and downward
    // directions due to the implemetation in Field_F.cpp.
    // The present counting is based on rev.1130. [10 Sep 2014 H.Matsufuru]

    const int Nvol = CommonParameters::Nvol();
    const int NPE  = CommonParameters::NPE();

    const int Nc = m_Nc;
    const int Nd = m_Nd;

    int flop_per_site = Nc * Nd * 2 * 8 * (4 * Nc - 1); // #(mult_Field_Gn/d)

    flop_per_site += Nc * Nd * 2 * (4 * 3 + 4 * 2);     // #(mult_GMproj2)
    flop_per_site += Nc * Nd * 2 * 8;                   // #(addpart_ex)
    flop_per_site += Nc * Nd * 2 * 2;                   // #(aypx(kappa))

    double gflop = flop_per_site * (Nvol * (NPE / 1.0e+9));

    if ((m_mode == "DdagD") || (m_mode == "DDdag")) {
      gflop *= 2;
    }

    return gflop;
  }


//====================================================================
}
//============================================================END=====
#endif
