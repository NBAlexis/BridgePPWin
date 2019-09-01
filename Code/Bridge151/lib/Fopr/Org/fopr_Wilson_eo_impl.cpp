/*!
        @file    fopr_Wilson_eo_impl.cpp

        @brief

        @author  UEDA, Satoru  (sueda)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/
#include "BridgeLib_Private.h"
#if !USE_IMP
#include "fopr_Wilson_eo_impl.h"

namespace Org {
  const std::string Fopr_Wilson_eo::class_name = "Org::Fopr_Wilson_eo";

#ifdef USE_FACTORY_AUTOREGISTER
  namespace {
    bool init = Fopr_Wilson_eo::register_factory();
  }
#endif

//===================================================================
  namespace {
    //! U: even-odd preconditioned gauge field (full size)
    //! f: even-odd preconditioned fermion field (half size)
    //! w = U_ieo * f_ieo
    void mult_Field_Gn_eo(Field_F& w, int ex,
                          const Field_G& U, const int ex1,
                          const Field_F& x, const int ex2,
                          const int ieo)
    {
      assert(ex < w.nex());
      assert(ex1 < U.nex());
      assert(ex2 < x.nex());
      assert(U.nvol() == w.nvol() * 2);
      assert(x.nvol() == w.nvol());

      for (int site = 0, nvol = w.nvol(); site < nvol; ++site) {
        for (int s = 0, nd = w.nd(); s < nd; ++s) {
          Vec_SU_N vec(w.nc());
          vec = U.mat(site + ieo * nvol, ex1) * x.vec(s, site, ex2);
          w.set_vec(s, site, ex, vec);
        }
      }
    }


    //====================================================================
    //! U: even-odd preconditioned gauge field (full size)
    //! f: even-odd preconditioned fermion field (half size)
    //! w = U^dag_ieo * f_ieo
    void mult_Field_Gd_eo(Field_F& w, int ex,
                          const Field_G& U, int ex1,
                          const Field_F& x, int ex2,
                          const int ieo)
    {
      assert(ex < w.nex());
      assert(ex1 < U.nex());
      assert(ex2 < x.nex());
      assert(U.nvol() == w.nvol() * 2);
      assert(x.nvol() == w.nvol());

      for (int site = 0, nvol = w.nvol(); site < nvol; ++site) {
        for (int s = 0, nd = w.nd(); s < nd; ++s) {
          Vec_SU_N vec(w.nc());
          vec = U.mat_dag(site + ieo * nvol, ex1) * x.vec(s, site, ex2);
          w.set_vec(s, site, ex, vec);
        }
      }
    }
  }

//====================================================================
  void Fopr_Wilson_eo::set_parameters(const Parameters& params)
  {
    const string str_vlevel = params.get_string("verbose_level");

    m_vl = vout.set_verbose_level(str_vlevel);

    //- fetch and check input parameters
    double           kappa;
    std::vector<int> bc;

    int err = 0;
    err += params.fetch_double("hopping_parameter", kappa);
    err += params.fetch_int_vector("boundary_condition", bc);

    if (err) {
      vout.crucial(m_vl, "Error at %s: fetch error, input parameter not found.\n", class_name.c_str());
      exit(EXIT_FAILURE);
    }

    set_parameters(kappa, bc);
  }


//====================================================================
  void Fopr_Wilson_eo::set_parameters(const double kappa, const std::vector<int> bc)
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
  void Fopr_Wilson_eo::set_config(Field *U)
  {
    const int Ndim = CommonParameters::Ndim();
    const int Nvol = CommonParameters::Nvol();

    m_index.convertField(*m_Ueo, *U);
  }


//====================================================================
  void Fopr_Wilson_eo::set_mode(const std::string mode)
  {
    m_mode = mode;

    if (m_mode == "D") {
      m_mult     = &Fopr_Wilson_eo::D;
      m_mult_dag = &Fopr_Wilson_eo::Ddag;
      m_preProp  = &Fopr_Wilson_eo::prePropD;
      m_postProp = &Fopr_Wilson_eo::postPropD;
    } else if (m_mode == "Ddag") {
      m_mult     = &Fopr_Wilson_eo::Ddag;
      m_mult_dag = &Fopr_Wilson_eo::D;
      m_preProp  = &Fopr_Wilson_eo::prePropDag;
      m_postProp = &Fopr_Wilson_eo::postPropDag;
    } else if (m_mode == "DdagD") {
      m_mult     = &Fopr_Wilson_eo::DdagD;
      m_mult_dag = &Fopr_Wilson_eo::DdagD;
    } else if (m_mode == "DDdag") {
      m_mult     = &Fopr_Wilson_eo::DDdag;
      m_mult_dag = &Fopr_Wilson_eo::DDdag;
    } else if (m_mode == "H") {
      m_mult     = &Fopr_Wilson_eo::H;
      m_mult_dag = &Fopr_Wilson_eo::H;
    } else {
      vout.crucial("Error at %s: input mode is undefined.\n", class_name.c_str());
      exit(EXIT_FAILURE);
    }
  }


//====================================================================
  std::string Fopr_Wilson_eo::get_mode() const
  {
    return m_mode;
  }


//====================================================================
  void Fopr_Wilson_eo::init(const std::string repr)
  {
    m_repr = repr;

    m_boundary.resize(m_Ndim);
    m_GM.resize(m_Ndim + 1);

    unique_ptr<GammaMatrixSet> gmset(GammaMatrixSet::New(m_repr));

    m_GM[0] = gmset->get_GM(gmset->GAMMA1);
    m_GM[1] = gmset->get_GM(gmset->GAMMA2);
    m_GM[2] = gmset->get_GM(gmset->GAMMA3);
    m_GM[3] = gmset->get_GM(gmset->GAMMA4);
    m_GM[4] = gmset->get_GM(gmset->GAMMA5);

    m_Ueo = new Field_G(m_Nvol, m_Ndim);
  }


//====================================================================
  void Fopr_Wilson_eo::prePropD(Field& Be, Field& bo, const Field& b)
  {
#pragma omp master
    {
      const int Nin = b.nin();
      const int Nex = b.nex();

      Field be(Nin, m_Nvol2, Nex);
      m_index.convertField(be, b, 0);

      m_index.convertField(bo, b, 1);

      Be.reset(Nin, m_Nvol2, Nex);

      //  Be = be - (Field)Meo(bo, 0);
      copy(Be, be);

      Field_F v1(m_Nvol2, Nex);
      Meo(v1, bo, 0);

      axpy(Be, -1.0, v1);
    }
#pragma omp barrier
  }


//====================================================================
  void Fopr_Wilson_eo::postPropD(Field& x, const Field& xe, const Field& bo)
  {
#pragma omp master
    {
      const int Nin = xe.nin();
      const int Nex = xe.nex();

      //  xo = bo - (Field)Meo(xe, 1);
      Field xo(Nin, m_Nvol2, Nex);
      copy(xo, bo);

      Field_F v1(m_Nvol2, Nex);
      Meo(v1, xe, 1);

      axpy(xo, -1.0, v1);

      x.reset(Nin, m_Nvol2 * 2, Nex);
      m_index.reverseField(x, xe, 0);
      m_index.reverseField(x, xo, 1);
    }
#pragma omp barrier
  }


//====================================================================
  void Fopr_Wilson_eo::prePropDag(Field& Be, Field& bo, const Field& b)
  {
#pragma omp master
    {
      const int Nin = b.nin();
      const int Nex = b.nex();

      Field be(Nin, m_Nvol2, Nex);
      m_index.convertField(be, b, 0);

      m_index.convertField(bo, b, 1);

      Be.reset(Nin, m_Nvol2, Nex);

      //  Be = be - (Field)Mdageo(bo, 0);
      copy(Be, be);

      Field_F v1(m_Nvol2, Nex);
      Mdageo(v1, bo, 0);

      axpy(Be, -1.0, v1);
    }
#pragma omp barrier
  }


//====================================================================
  void Fopr_Wilson_eo::postPropDag(Field& x, const Field& xe, const Field& bo)
  {
#pragma omp master
    {
      //Field bo();
      //index.convertField(bo, b, 1);
      const int Nin = xe.nin();
      const int Nex = xe.nex();

      //  xo = bo - (Field)Mdageo(xe, 1);
      Field xo(Nin, m_Nvol2, Nex);
      copy(xo, bo);

      Field_F v1(m_Nvol2, Nex);
      Mdageo(v1, xe, 1);
      axpy(xo, -1.0, v1);

      x.reset(Nin, m_Nvol2 * 2, Nex);
      m_index.reverseField(x, xe, 0);
      m_index.reverseField(x, xo, 1);
    }
#pragma omp barrier
  }


//====================================================================
  void Fopr_Wilson_eo::D(Field& v, const Field& f)
  {
#pragma omp master
    {
      assert(f.nex() == 1);

      copy(v, f); //  v  = f;

      //  v -= (Field)Meo(Meo(f, 1), 0);
      Field_F v1(m_Nvol2, 1);
      Meo(v1, f, 1);

      Field_F v2(m_Nvol2, 1);
      Meo(v2, v1, 0);

      axpy(v, -1.0, v2);
    }
#pragma omp barrier
  }


//====================================================================
  void Fopr_Wilson_eo::Ddag(Field& v, const Field& f)
  {
#pragma omp master
    {
      copy(v, f); //  v  = f;

      //  v -= (Field)Mdageo(Mdageo(f, 1), 0);
      Field_F v1(m_Nvol2, 1);
      Mdageo(v1, f, 1);

      Field_F v2(m_Nvol2, 1);
      Mdageo(v2, v1, 0);

      axpy(v, -1.0, v2);
    }
#pragma omp barrier
  }


//====================================================================
  void Fopr_Wilson_eo::DdagD(Field& v, const Field& f)
  {
#pragma omp master
    {
      Field w(v);
      D(w, f);
      Ddag(v, w);
    }
#pragma omp barrier
  }


//====================================================================
  void Fopr_Wilson_eo::DDdag(Field& v, const Field& f)
  {
#pragma omp master
    {
      Field w(v);
      Ddag(w, f);
      D(v, w);
    }
#pragma omp barrier
  }


//====================================================================
  void Fopr_Wilson_eo::H(Field& v, const Field& f)
  {
#pragma omp master
    {
      Field w(v);
      D(w, f);
      mult_gm5(v, w);
    }
#pragma omp barrier
  }


//====================================================================
  void Fopr_Wilson_eo::MeoMoe(Field& v, const Field& f)
  {
#pragma omp master
    {
      //  v -= (Field)Meo(Meo(f, 1), 0);
      Field_F v1(m_Nvol2, 1);
      Meo(v1, f, 1);

      Field_F v2(m_Nvol2, 1);
      Meo(v2, v1, 0);

      axpy(v, -1.0, v2);
    }
#pragma omp barrier
  }


//====================================================================
  void Fopr_Wilson_eo::Meo(Field& w,
                           const Field& f, const int ieo)
  {
#pragma omp master
    {
      Field_F f2(m_Nvol2, 1);
      copy(f2, f);

      Field_F w2(m_Nvol2, 1);
      w2.set(0.0); // w2 = 0.0

      for (int mu = 0; mu < m_Ndim; ++mu) {
        mult_p(mu, w2, f2, ieo);
        mult_m(mu, w2, f2, ieo);
      }

      scal(w2, -m_kappa); //  w *= -m_kappa;

      copy(w, w2);
    }
  }


//====================================================================
  void Fopr_Wilson_eo::Mdageo(Field& w,
                              const Field& f, const int ieo)
  {
#pragma omp master
    {
      Field_F v(m_Nvol2, f.nex());
      mult_gm5(v, f);
      Meo_gm5(w, v, ieo);
    }
  }


//====================================================================
  void Fopr_Wilson_eo::Meo_gm5(Field& v,
                               const Field& f, const int ieo)
  {
    Field_F w(m_Nvol2, f.nex());

    Meo(w, f, ieo);
    mult_gm5(v, w);
  }


//====================================================================
  void Fopr_Wilson_eo::mult_gm5(Field& w, const Field& f)
  {
#pragma omp master
    {
      Field_F w2 = (Field_F)f;
      mult_GM(w2, m_GM[4], (Field_F)f);
      w = w2;
    }
#pragma omp barrier
  }


//====================================================================
  void Fopr_Wilson_eo::gm5p(const int mu,
                            Field& w, const Field& f)
  {
    for (int ex = 0; ex < f.nex(); ++ex) {
      Field_F vt1(f.nvol(), 1);
      vt1.setpart_ex(0, f, ex);

      Field_F vt2(f.nvol(), 1);
      shift.backward(vt2, vt1, m_boundary[mu], mu); // vt2 = vt1(x+mu)

      mult_GMproj2(vt1, -1, m_GM[mu], vt2);         // vt1 = (1 - gamma_mu) vt2
      mult_GM(vt2, m_GM[4], vt1);                   // vt2 = gamma_5 vt1

      w.addpart_ex(ex, vt2, 0);
    }
  }


//====================================================================
  void Fopr_Wilson_eo::mult_p(const int mu,
                              Field_F& w, const Field_F& f,
                              const int ieo)
  {
    for (int ex = 0; ex < f.nex(); ++ex) {
      Field_F vt(f.nvol(), 1);
      vt.setpart_ex(0, f, ex);

      shift.backward_h(trf, vt, m_boundary[mu], mu, ieo);

      mult_Field_Gn_eo(trf2, 0, *m_Ueo, mu, trf, 0, ieo);
      mult_GMproj2(vt, -1, m_GM[mu], trf2);

      w.addpart_ex(ex, vt, 0);
    }
  }


//=====================================================================
  void Fopr_Wilson_eo::mult_m(const int mu,
                              Field_F& w, const Field_F& f,
                              const int ieo)
  {
    for (int ex = 0; ex < f.nex(); ++ex) {
      mult_Field_Gd_eo(trf, 0, *m_Ueo, mu, f, ex, 1 - ieo);
      shift.forward_h(trf2, trf, m_boundary[mu], mu, ieo);

      Field_F vt(f.nvol(), 1);
      mult_GMproj2(vt, 1, m_GM[mu], trf2);

      w.addpart_ex(ex, vt, 0);
    }
  }


//====================================================================
  double Fopr_Wilson_eo::flop_count()
  {
    // Counting of floating point operations in giga unit.
    // this counting is based on the Org-implementation of ver.1.2.0.
    // Flop count of mult_GMproj2 is different for upward and downward
    // directions due to the implemetation in Field_F.cpp.
    // Present counting is based on rev.1130. [10 Sep 2014 H.Matsufuru]

    const int Nvol = CommonParameters::Nvol();
    const int NPE  = CommonParameters::NPE();

    const int Nc = m_Nc;
    const int Nd = m_Nd;

    int flop_Meo = Nc * Nd * 2 * 8 * (4 * Nc - 1);            // #(mult_Field_Gn/d)

    flop_Meo += Nc * Nd * 2 * (4 * 3 + 4 * 2);                // #(mult_GMproj2)
    flop_Meo += Nc * Nd * 2 * 8;                              // #(addpart_ex)
    flop_Meo += Nc * Nd * 2;                                  // #(scal(kappa))

    const int flop_per_site = 2 * flop_Meo + Nc * Nd * 2 * 2; // #(2*Meo + axpy)

    double gflop = flop_per_site * ((Nvol / 2) * (NPE / 1.0e+9));

    if ((m_mode == "DdagD") || (m_mode == "DDdag")) {
      gflop *= 2;
    }

    return gflop;
  }


//====================================================================
}
//============================================================END=====
#endif
