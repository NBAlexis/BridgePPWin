/*!
        @file    fopr_Wilson_eo_impl.cpp

        @brief

        @author  UEDA, Satoru  (sueda)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/
#include "BridgeLib_Private.h"
#if USE_IMP
#include "fopr_Wilson_eo_impl.h"

#include "ResourceManager/threadManager_OpenMP.h"

namespace Imp {
//====================================================================
  const std::string Fopr_Wilson_eo::class_name = "Imp::Fopr_Wilson_eo";

#ifdef USE_FACTORY_AUTOREGISTER
  namespace {
    bool init = Fopr_Wilson_eo::register_factory();
  }
#endif

//====================================================================
  void Fopr_Wilson_eo::init(const std::string repr)
  {
    m_repr = repr;

    m_vl = CommonParameters::Vlevel();

    m_Nc  = CommonParameters::Nc();
    m_Nd  = CommonParameters::Nd();
    m_Nvc = 2 * m_Nc;
    m_Ndf = 2 * m_Nc * m_Nc;
    const int Nvcd = m_Nvc * m_Nd;

    m_Nvol  = CommonParameters::Nvol();
    m_Nvol2 = m_Nvol / 2;
    m_Ndim  = CommonParameters::Ndim();

    m_Nx = CommonParameters::Nx();
    m_Ny = CommonParameters::Ny();
    m_Nz = CommonParameters::Nz();
    m_Nt = CommonParameters::Nt();

    if ((m_Nx % 2) != 0) {
      vout.crucial(m_vl, "Error at %s: Nx=%d must be even.\n", class_name.c_str(), m_Nx);
      exit(EXIT_FAILURE);
    }
    m_Nx2 = m_Nx / 2;

    if ((m_Ny % 2) != 0) {
      vout.crucial(m_vl, "Error at %s: Ny=%d must be even.\n", class_name.c_str(), m_Ny);
      exit(EXIT_FAILURE);
    }

    m_Nvol = CommonParameters::Nvol();
    m_Ndim = CommonParameters::Ndim();
    m_boundary.resize(m_Ndim);
    m_boundary_each_node.resize(m_Ndim);

    m_yzt_eo.resize(m_Ny * m_Nz * m_Nt);

    for (int t = 0; t < m_Nt; ++t) {
      for (int z = 0; z < m_Nz; ++z) {
        for (int y = 0; y < m_Ny; ++y) {
          int t_global = t + Communicator::ipe(3) * m_Nt;
          int z_global = z + Communicator::ipe(2) * m_Nz;
          int y_global = y + Communicator::ipe(1) * m_Ny;

          m_yzt_eo[y + m_Ny * (z + m_Nz * t)] = (y_global + z_global + t_global) % 2;
        }
      }
    }

    m_Ueo = new Field_G(m_Nvol, m_Ndim);

    if (m_repr == "Dirac") {
      m_gm5      = &Fopr_Wilson_eo::gm5_dirac;
      m_gm5_self = &Fopr_Wilson_eo::gm5_self_dirac;
      m_mult_tp  = &Fopr_Wilson_eo::mult_tp_dirac;
      m_mult_tm  = &Fopr_Wilson_eo::mult_tm_dirac;
    } else if (m_repr == "Chiral") {
      m_gm5      = &Fopr_Wilson_eo::gm5_chiral;
      m_gm5_self = &Fopr_Wilson_eo::gm5_self_chiral;
      m_mult_tp  = &Fopr_Wilson_eo::mult_tp_chiral;
      m_mult_tm  = &Fopr_Wilson_eo::mult_tm_chiral;
    } else {
      vout.crucial(m_vl, "Error at %s: input repr is undefined.\n",
                   class_name.c_str());
      exit(EXIT_FAILURE);
    }

    const int Nvx = m_Nvc * 2 * (m_Ny / 2) * m_Nz * m_Nt;
    vcp1_xp = new double[Nvx];
    vcp2_xp = new double[Nvx];
    vcp1_xm = new double[Nvx];
    vcp2_xm = new double[Nvx];

    const int Nvy = m_Nvc * 2 * m_Nx2 * m_Nz * m_Nt;
    vcp1_yp = new double[Nvy];
    vcp2_yp = new double[Nvy];
    vcp1_ym = new double[Nvy];
    vcp2_ym = new double[Nvy];

    const int Nvz = m_Nvc * 2 * m_Nx2 * m_Ny * m_Nt;
    vcp1_zp = new double[Nvz];
    vcp2_zp = new double[Nvz];
    vcp1_zm = new double[Nvz];
    vcp2_zm = new double[Nvz];

    const int Nvt = m_Nvc * 2 * m_Nx2 * m_Ny * m_Nz;
    vcp1_tp = new double[Nvt];
    vcp2_tp = new double[Nvt];
    vcp1_tm = new double[Nvt];
    vcp2_tm = new double[Nvt];

    m_v1.reset(Nvcd, m_Nvol2, 1);
    m_v2.reset(Nvcd, m_Nvol2, 1);
    m_w1.reset(Nvcd, m_Nvol2, 1);
    m_w2.reset(Nvcd, m_Nvol2, 1);

    setup_thread();
  }


//====================================================================
  void Fopr_Wilson_eo::tidyup()
  {
    delete m_Ueo;

    delete[]  vcp1_xp;
    delete[]  vcp2_xp;
    delete[]  vcp1_xm;
    delete[]  vcp2_xm;

    delete[]  vcp1_yp;
    delete[]  vcp2_yp;
    delete[]  vcp1_ym;
    delete[]  vcp2_ym;

    delete[]  vcp1_zp;
    delete[]  vcp2_zp;
    delete[]  vcp1_zm;
    delete[]  vcp2_zm;

    delete[]  vcp1_tp;
    delete[]  vcp2_tp;
    delete[]  vcp1_tm;
    delete[]  vcp2_tm;
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
  void Fopr_Wilson_eo::set_parameters(
    const double kappa,
    const std::vector<int> bc)
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

    // m_boundary.resize(m_Ndim);  // already resized in init.
    m_boundary = bc;

    for (int idir = 0; idir < m_Ndim; ++idir) {
      m_boundary_each_node[idir] = 1.0;
      if (Communicator::ipe(idir) == 0) m_boundary_each_node[idir] = m_boundary[idir];
    }
  }


//====================================================================
  void Fopr_Wilson_eo::set_config(Field *U)
  {
    m_index.convertField(*m_Ueo, *U);
    m_U = m_Ueo;
  }


//====================================================================
  void Fopr_Wilson_eo::set_mode(std::string mode)
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
  void Fopr_Wilson_eo::prePropD(Field& Be, Field& bo,
                                const Field& b)
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

#pragma omp barrier
  }


//====================================================================
  void Fopr_Wilson_eo::postPropD(Field& x,
                                 const Field& xe, const Field& bo)
  {
    const int Nin = xe.nin();
    const int Nex = xe.nex();

    Field xo(Nin, m_Nvol2, Nex);

    copy(xo, bo);

    Field_F v1(m_Nvol2, Nex);
    Meo(v1, xe, 1);
    axpy(xo, -1.0, v1);

    x.reset(Nin, m_Nvol2 * 2, Nex);
    m_index.reverseField(x, xe, 0);
    m_index.reverseField(x, xo, 1);

#pragma omp barrier
  }


//====================================================================
  void Fopr_Wilson_eo::prePropDag(Field& Be,
                                  Field& bo, const Field& b)
  {
    const int Nin = b.nin();
    const int Nex = b.nex();

    Field be(Nin, m_Nvol2, Nex);

    m_index.convertField(be, b, 0);
    m_index.convertField(bo, b, 1);

    Be.reset(Nin, m_Nvol2, Nex);

    copy(Be, be);

    Field_F v1(m_Nvol2, Nex);
    Mdageo(v1, bo, 0);

    axpy(Be, -1.0, v1);

#pragma omp barrier
  }


//====================================================================
  void Fopr_Wilson_eo::postPropDag(Field& x,
                                   const Field& xe, const Field& bo)
  {
    const int Nin = xe.nin();
    const int Nex = xe.nex();

    Field xo(Nin, m_Nvol2, Nex);

    copy(xo, bo);

    Field_F v1(m_Nvol2, Nex);
    Mdageo(v1, xe, 1);
    axpy(xo, -1.0, v1);

    x.reset(Nin, m_Nvol2 * 2, Nex);
    m_index.reverseField(x, xe, 0);
    m_index.reverseField(x, xo, 1);

#pragma omp barrier
  }


//====================================================================
  void Fopr_Wilson_eo::DdagD(Field& v, const Field& f)
  {
    D(m_w1, f);
    Ddag(v, m_w1);
  }


//====================================================================
  void Fopr_Wilson_eo::DDdag(Field& v, const Field& f)
  {
    Ddag(m_w1, f);
    D(v, m_w1);
  }


//====================================================================
  void Fopr_Wilson_eo::H(Field& v, const Field& f)
  {
    D(v, f);
    mult_gm5(v);
  }


//====================================================================
  void Fopr_Wilson_eo::MeoMoe(Field& v, const Field& f)
  {
    Meo(m_w1, f, 1);
    Meo(m_w2, m_w1, 0);
    axpy(v, -1.0, m_w2);

#pragma omp barrier
  }


//====================================================================
  void Fopr_Wilson_eo::D(Field& v, const Field& f)
  {
    assert(f.nex() == 1);

    Meo(m_v1, f, 1);
    Meo(v, m_v1, 0);
    aypx(-1.0, v, f);

#pragma omp barrier
  }


//====================================================================
  void Fopr_Wilson_eo::Ddag(Field& v, const Field& f)
  {
    Mdageo(m_v1, f, 1);
    Mdageo(v, m_v1, 0);
    aypx(-1.0, v, f);

#pragma omp barrier
  }


//====================================================================
  void Fopr_Wilson_eo::Meo(Field& w,
                           const Field& f, const int ieo)
  {
#pragma omp barrier

    //  clear_impl(w);  //  w = 0.0;
    w.set(0.0); //  w = 0.0;

    mult_xp(w, f, ieo);
    mult_xm(w, f, ieo);

    mult_yp(w, f, ieo);
    mult_ym(w, f, ieo);

    mult_zp(w, f, ieo);
    mult_zm(w, f, ieo);

    (this->*m_mult_tp)(w, f, ieo);
    (this->*m_mult_tm)(w, f, ieo);

    scal(w, -m_kappa); //  w *= -m_kappa;  using Field function.
    //  scal_impl(w, -m_kappa);   //  w *= -m_kappa;

    ThreadManager_OpenMP::sync_barrier_all();
  }


//====================================================================
  void Fopr_Wilson_eo::Mdageo(Field& w,
                              const Field& f, const int ieo)
  {
    //  Field_F v(m_Nvol2, f.nex());
    //  mult_gm5(v, (Field)f);
    //  Meo_gm5(w, v, ieo);

    mult_gm5(m_v2, f);
    Meo(w, m_v2, ieo);
    mult_gm5(w);
  }


//====================================================================
  void Fopr_Wilson_eo::Meo_gm5(Field& v,
                               const Field& f, const int ieo)
  {
    Meo(v, f, ieo);
    mult_gm5(v);
  }


//====================================================================
  void Fopr_Wilson_eo::mult_gm5(Field& w, const Field& f)
  {
    (this->*m_gm5)(w, f);
  }


//====================================================================
  void Fopr_Wilson_eo::mult_gm5(Field& w)
  {
    (this->*m_gm5_self)(w);
  }


//====================================================================
  void Fopr_Wilson_eo::gm5_dirac(Field& w,
                                 const Field& f)
  {
    const double *v1 = f.ptr(0);
    double       *v2 = w.ptr(0);

    const int Nthread  = ThreadManager_OpenMP::get_num_threads();
    const int i_thread = ThreadManager_OpenMP::get_thread_id();

    const int is = m_Ntask * i_thread / Nthread;
    const int ns = m_Ntask * (i_thread + 1) / Nthread - is;

    for (int i = is; i < is + ns; ++i) {
      gm5_dirac_thread(i, v2, v1);
    }
#pragma omp barrier
  }


//====================================================================
  void Fopr_Wilson_eo::gm5_chiral(Field& w,
                                  const Field& f)
  {
    const double *v1 = f.ptr(0);
    double       *v2 = w.ptr(0);

    const int Nthread  = ThreadManager_OpenMP::get_num_threads();
    const int i_thread = ThreadManager_OpenMP::get_thread_id();

    const int is = m_Ntask * i_thread / Nthread;
    const int ns = m_Ntask * (i_thread + 1) / Nthread - is;

    for (int i = is; i < is + ns; ++i) {
      gm5_chiral_thread(i, v2, v1);
    }
#pragma omp barrier
  }


//====================================================================
  void Fopr_Wilson_eo::gm5_self_dirac(Field& w)
  {
    double *wp = w.ptr(0);

    const int Nthread  = ThreadManager_OpenMP::get_num_threads();
    const int i_thread = ThreadManager_OpenMP::get_thread_id();

    const int is = m_Ntask * i_thread / Nthread;
    const int ns = m_Ntask * (i_thread + 1) / Nthread - is;

    for (int i = is; i < is + ns; ++i) {
      gm5_dirac_thread(i, wp);
    }
#pragma omp barrier
  }


//====================================================================
  void Fopr_Wilson_eo::gm5_self_chiral(Field& w)
  {
    double *wp = w.ptr(0);

    const int Nthread  = ThreadManager_OpenMP::get_num_threads();
    const int i_thread = ThreadManager_OpenMP::get_thread_id();

    const int is = m_Ntask * i_thread / Nthread;
    const int ns = m_Ntask * (i_thread + 1) / Nthread - is;

    for (int i = is; i < is + ns; ++i) {
      gm5_chiral_thread(i, wp);
    }
#pragma omp barrier
  }


//====================================================================
  void Fopr_Wilson_eo::gm5p(const int mu, Field& w,
                            const Field& f)
  {
    // this function is probably not to be used.
    // determines  \gamma_5 * (1 - \gamma_\mu) v(x+\hat{x})
    vout.crucial(m_vl, "Error at %s: gm5p is undefined.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }


//====================================================================
  void Fopr_Wilson_eo::clear_impl(Field& w)
  {
    double *wp = w.ptr(0);

    const int Nthread  = ThreadManager_OpenMP::get_num_threads();
    const int i_thread = ThreadManager_OpenMP::get_thread_id();

    const int is = m_Ntask * i_thread / Nthread;
    const int ns = m_Ntask * (i_thread + 1) / Nthread - is;

    for (int i = is; i < is + ns; ++i) {
      clear_thread(i, wp);
    }
#pragma omp barrier
  }


//====================================================================
  void Fopr_Wilson_eo::scal_impl(Field& w, double a)
  {
    double *wp = w.ptr(0);

    const int Nthread  = ThreadManager_OpenMP::get_num_threads();
    const int i_thread = ThreadManager_OpenMP::get_thread_id();

    const int is = m_Ntask * i_thread / Nthread;
    const int ns = m_Ntask * (i_thread + 1) / Nthread - is;

    for (int i = is; i < is + ns; ++i) {
      scal_thread(i, wp, a);
    }
#pragma omp barrier
  }


//====================================================================
  void Fopr_Wilson_eo::mult_xp(Field& w,
                               const Field& f, const int ieo)
  {
    const double *v1 = f.ptr(0);
    double       *v2 = w.ptr(0);

    const int Nthread  = ThreadManager_OpenMP::get_num_threads();
    const int i_thread = ThreadManager_OpenMP::get_thread_id();

    const int is = m_Ntask * i_thread / Nthread;
    const int ns = m_Ntask * (i_thread + 1) / Nthread - is;

    for (int i = is; i < is + ns; ++i) {
      mult_xp1_thread(i, vcp1_xp, v1, ieo);
    }

#pragma omp barrier
#pragma omp master
    {
      const int Nv = m_Nvc * 2 * (m_Ny / 2) * m_Nz * m_Nt;
      Communicator::exchange(Nv, vcp2_xp, vcp1_xp, 0, 1, 1);
    }
#pragma omp barrier

    for (int i = is; i < is + ns; ++i) {
      mult_xpb_thread(i, v2, v1, ieo);
    }

    for (int i = is; i < is + ns; ++i) {
      mult_xp2_thread(i, v2, vcp2_xp, ieo);
    }
#pragma omp barrier
  }


//====================================================================
  void Fopr_Wilson_eo::mult_xm(Field& w,
                               const Field& f, const int ieo)
  {
    const double *v1 = f.ptr(0);
    double       *v2 = w.ptr(0);

    const int Nthread  = ThreadManager_OpenMP::get_num_threads();
    const int i_thread = ThreadManager_OpenMP::get_thread_id();

    const int is = m_Ntask * i_thread / Nthread;
    const int ns = m_Ntask * (i_thread + 1) / Nthread - is;

    for (int i = is; i < is + ns; ++i) {
      mult_xm1_thread(i, vcp1_xm, v1, ieo);
    }

#pragma omp barrier
#pragma omp master
    {
      const int Nv = m_Nvc * 2 * (m_Ny / 2) * m_Nz * m_Nt;
      Communicator::exchange(Nv, vcp2_xm, vcp1_xm, 0, -1, 2);
    }
#pragma omp barrier

    for (int i = is; i < is + ns; ++i) {
      mult_xmb_thread(i, v2, v1, ieo);
    }

    for (int i = is; i < is + ns; ++i) {
      mult_xm2_thread(i, v2, vcp2_xm, ieo);
    }
#pragma omp barrier
  }


//====================================================================
  void Fopr_Wilson_eo::mult_yp(Field& w,
                               const Field& f, const int ieo)
  {
    const double *v1 = f.ptr(0);
    double       *v2 = w.ptr(0);

    const int Nthread  = ThreadManager_OpenMP::get_num_threads();
    const int i_thread = ThreadManager_OpenMP::get_thread_id();

    const int is = m_Ntask * i_thread / Nthread;
    const int ns = m_Ntask * (i_thread + 1) / Nthread - is;

    for (int i = is; i < is + ns; ++i) {
      mult_yp1_thread(i, vcp1_yp, v1, ieo);
    }

#pragma omp barrier
#pragma omp master
    {
      const int Nv = m_Nvc * 2 * m_Nx2 * m_Nz * m_Nt;
      Communicator::exchange(Nv, vcp2_yp, vcp1_yp, 1, 1, 3);
    }
#pragma omp barrier

    for (int i = is; i < is + ns; ++i) {
      mult_ypb_thread(i, v2, v1, ieo);
    }

    for (int i = is; i < is + ns; ++i) {
      mult_yp2_thread(i, v2, vcp2_yp, ieo);
    }
#pragma omp barrier
  }


//====================================================================
  void Fopr_Wilson_eo::mult_ym(Field& w,
                               const Field& f, const int ieo)
  {
    const double *v1 = f.ptr(0);
    double       *v2 = w.ptr(0);

    const int Nthread  = ThreadManager_OpenMP::get_num_threads();
    const int i_thread = ThreadManager_OpenMP::get_thread_id();

    const int is = m_Ntask * i_thread / Nthread;
    const int ns = m_Ntask * (i_thread + 1) / Nthread - is;

    for (int i = is; i < is + ns; ++i) {
      mult_ym1_thread(i, vcp1_ym, v1, ieo);
    }

#pragma omp barrier
#pragma omp master
    {
      const int Nv = m_Nvc * 2 * m_Nx2 * m_Nz * m_Nt;
      Communicator::exchange(Nv, vcp2_ym, vcp1_ym, 1, -1, 4);
    }
#pragma omp barrier

    for (int i = is; i < is + ns; ++i) {
      mult_ymb_thread(i, v2, v1, ieo);
    }

    for (int i = is; i < is + ns; ++i) {
      mult_ym2_thread(i, v2, vcp2_ym, ieo);
    }
#pragma omp barrier
  }


//====================================================================
  void Fopr_Wilson_eo::mult_zp(Field& w,
                               const Field& f, const int ieo)
  {
    const double *v1 = f.ptr(0);
    double       *v2 = w.ptr(0);

    const int Nthread  = ThreadManager_OpenMP::get_num_threads();
    const int i_thread = ThreadManager_OpenMP::get_thread_id();

    const int is = m_Ntask * i_thread / Nthread;
    const int ns = m_Ntask * (i_thread + 1) / Nthread - is;

    for (int i = is; i < is + ns; ++i) {
      mult_zp1_thread(i, vcp1_zp, v1, ieo);
    }

#pragma omp barrier
#pragma omp master
    {
      const int Nv = m_Nvc * 2 * m_Nx2 * m_Ny * m_Nt;
      Communicator::exchange(Nv, vcp2_zp, vcp1_zp, 2, 1, 5);
    }
#pragma omp barrier

    for (int i = is; i < is + ns; ++i) {
      mult_zpb_thread(i, v2, v1, ieo);
    }

    for (int i = is; i < is + ns; ++i) {
      mult_zp2_thread(i, v2, vcp2_zp, ieo);
    }
#pragma omp barrier
  }


//====================================================================
  void Fopr_Wilson_eo::mult_zm(Field& w,
                               const Field& f, const int ieo)
  {
    const double *v1 = f.ptr(0);
    double       *v2 = w.ptr(0);

    const int Nthread  = ThreadManager_OpenMP::get_num_threads();
    const int i_thread = ThreadManager_OpenMP::get_thread_id();

    const int is = m_Ntask * i_thread / Nthread;
    const int ns = m_Ntask * (i_thread + 1) / Nthread - is;

    for (int i = is; i < is + ns; ++i) {
      mult_zm1_thread(i, vcp1_zm, v1, ieo);
    }

#pragma omp barrier
#pragma omp master
    {
      const int Nv = m_Nvc * 2 * m_Nx2 * m_Ny * m_Nt;
      Communicator::exchange(Nv, vcp2_zm, vcp1_zm, 2, -1, 6);
    }
#pragma omp barrier

    for (int i = is; i < is + ns; ++i) {
      mult_zmb_thread(i, v2, v1, ieo);
    }

    for (int i = is; i < is + ns; ++i) {
      mult_zm2_thread(i, v2, vcp2_zm, ieo);
    }
#pragma omp barrier
  }


//====================================================================
  void Fopr_Wilson_eo::mult_tp_dirac(Field& w,
                                     const Field& f, const int ieo)
  {
    const double *v1 = f.ptr(0);
    double       *v2 = w.ptr(0);

    const int Nthread  = ThreadManager_OpenMP::get_num_threads();
    const int i_thread = ThreadManager_OpenMP::get_thread_id();

    const int is = m_Ntask * i_thread / Nthread;
    const int ns = m_Ntask * (i_thread + 1) / Nthread - is;

    for (int i = is; i < is + ns; ++i) {
      mult_tp1_dirac_thread(i, vcp1_tp, v1, ieo);
    }

#pragma omp barrier
#pragma omp master
    {
      const int Nv = m_Nvc * 2 * m_Nx2 * m_Ny * m_Nz;
      Communicator::exchange(Nv, vcp2_tp, vcp1_tp, 3, 1, 7);
    }
#pragma omp barrier

    for (int i = is; i < is + ns; ++i) {
      mult_tpb_dirac_thread(i, v2, v1, ieo);
    }

    for (int i = is; i < is + ns; ++i) {
      mult_tp2_dirac_thread(i, v2, vcp2_tp, ieo);
    }
#pragma omp barrier
  }


//====================================================================
  void Fopr_Wilson_eo::mult_tm_dirac(Field& w,
                                     const Field& f, const int ieo)
  {
    const double *v1 = f.ptr(0);
    double       *v2 = w.ptr(0);

    const int Nthread  = ThreadManager_OpenMP::get_num_threads();
    const int i_thread = ThreadManager_OpenMP::get_thread_id();

    const int is = m_Ntask * i_thread / Nthread;
    const int ns = m_Ntask * (i_thread + 1) / Nthread - is;

    for (int i = is; i < is + ns; ++i) {
      mult_tm1_dirac_thread(i, vcp1_tm, v1, ieo);
    }

#pragma omp barrier
#pragma omp master
    {
      const int Nv = m_Nvc * 2 * m_Nx2 * m_Ny * m_Nz;
      Communicator::exchange(Nv, vcp2_tm, vcp1_tm, 3, -1, 8);
    }
#pragma omp barrier

    for (int i = is; i < is + ns; ++i) {
      mult_tmb_dirac_thread(i, v2, v1, ieo);
    }

    for (int i = is; i < is + ns; ++i) {
      mult_tm2_dirac_thread(i, v2, vcp2_tm, ieo);
    }
#pragma omp barrier
  }


//====================================================================
  void Fopr_Wilson_eo::mult_tp_chiral(Field& w,
                                      const Field& f, const int ieo)
  {
    const double *v1 = f.ptr(0);
    double       *v2 = w.ptr(0);

    const int Nthread  = ThreadManager_OpenMP::get_num_threads();
    const int i_thread = ThreadManager_OpenMP::get_thread_id();

    const int is = m_Ntask * i_thread / Nthread;
    const int ns = m_Ntask * (i_thread + 1) / Nthread - is;

    for (int i = is; i < is + ns; ++i) {
      mult_tp1_chiral_thread(i, vcp1_tp, v1, ieo);
    }

#pragma omp barrier
#pragma omp master
    {
      const int Nv = m_Nvc * 2 * m_Nx2 * m_Ny * m_Nz;
      Communicator::exchange(Nv, vcp2_tp, vcp1_tp, 3, 1, 7);
    }
#pragma omp barrier

    for (int i = is; i < is + ns; ++i) {
      mult_tpb_chiral_thread(i, v2, v1, ieo);
    }

    for (int i = is; i < is + ns; ++i) {
      mult_tp2_chiral_thread(i, v2, vcp2_tp, ieo);
    }
#pragma omp barrier
  }


//====================================================================
  void Fopr_Wilson_eo::mult_tm_chiral(Field& w,
                                      const Field& f, const int ieo)
  {
    const double *v1 = f.ptr(0);
    double       *v2 = w.ptr(0);

    const int Nthread  = ThreadManager_OpenMP::get_num_threads();
    const int i_thread = ThreadManager_OpenMP::get_thread_id();

    const int is = m_Ntask * i_thread / Nthread;
    const int ns = m_Ntask * (i_thread + 1) / Nthread - is;

    for (int i = is; i < is + ns; ++i) {
      mult_tm1_chiral_thread(i, vcp1_tm, v1, ieo);
    }

#pragma omp barrier
#pragma omp master
    {
      const int Nv = m_Nvc * 2 * m_Nx2 * m_Ny * m_Nz;
      Communicator::exchange(Nv, vcp2_tm, vcp1_tm, 3, -1, 8);
    }
#pragma omp barrier

    for (int i = is; i < is + ns; ++i) {
      mult_tmb_chiral_thread(i, v2, v1, ieo);
    }

    for (int i = is; i < is + ns; ++i) {
      mult_tm2_chiral_thread(i, v2, vcp2_tm, ieo);
    }
#pragma omp barrier
  }


//====================================================================
  double Fopr_Wilson_eo::flop_count()
  {
    // Counting of floating point operations in giga unit.
    // The following counting explicitly depends on the implementation.
    // It will be recalculated when the code is modified.
    // The present counting is based on rev.1107. [24 Aug 2014 H.Matsufuru]

    const int Nvol = CommonParameters::Nvol();
    const int NPE  = CommonParameters::NPE();

    int flop_site;

    if (m_repr == "Dirac") {
      flop_site = m_Nc * m_Nd * (6 * (4 * m_Nc + 2) + 2 * (4 * m_Nc + 1));
    } else if (m_repr == "Chiral") {
      flop_site = m_Nc * m_Nd * 8 * (4 * m_Nc + 2);
    } else {
      vout.crucial(m_vl, "Error at %s: input repr is undefined.\n",
                   class_name.c_str());
      exit(EXIT_FAILURE);
    }

    double gflop = flop_site * ((Nvol / 2) * (NPE / 1.0e+9));

    return gflop;
  }


//====================================================================
}
//============================================================END=====
#endif
