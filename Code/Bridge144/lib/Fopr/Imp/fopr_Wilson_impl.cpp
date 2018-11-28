#include "BridgeLib_Private.h"
#if USE_IMP

/*!
        @file    $Id:: fopr_Wilson_impl.cpp #$

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 00:49:46 #$

        @version $LastChangedRevision: 1561 $
*/

#include "fopr_Wilson_impl.h"

#include "ResourceManager/threadManager_OpenMP.h"

//#define USE_SU2

namespace Imp {
#if defined USE_GROUP_SU3
#include "fopr_Wilson_impl_SU3.inc"
#elif defined USE_GROUP_SU2
#include "fopr_Wilson_impl_SU2.inc"
#elif defined USE_GROUP_SU_N
#include "fopr_Wilson_impl_SU_N.inc"
#endif

  const std::string Fopr_Wilson::class_name = "Imp::Fopr_Wilson";

//====================================================================
  void Fopr_Wilson::init(std::string repr)
  {
    m_vl = CommonParameters::Vlevel();

    vout.general(m_vl, "%s: Construction of Wilson fermion operator.\n", class_name.c_str());

    check_Nc();

    m_Nc  = CommonParameters::Nc();
    m_Nd  = CommonParameters::Nd();
    m_Nvc = 2 * m_Nc;
    m_Ndf = 2 * m_Nc * m_Nc;

    m_Nx = CommonParameters::Nx();
    m_Ny = CommonParameters::Ny();
    m_Nz = CommonParameters::Nz();
    m_Nt = CommonParameters::Nt();

    m_Nvol = CommonParameters::Nvol();
    m_Ndim = CommonParameters::Ndim();
    m_boundary.resize(m_Ndim);
    m_boundary2.resize(m_Ndim);

    m_repr = repr;

    m_GM.resize(m_Ndim + 1);

    unique_ptr<GammaMatrixSet> gmset(GammaMatrixSet::New(m_repr));

    m_GM[0] = gmset->get_GM(GammaMatrixSet::GAMMA1);
    m_GM[1] = gmset->get_GM(GammaMatrixSet::GAMMA2);
    m_GM[2] = gmset->get_GM(GammaMatrixSet::GAMMA3);
    m_GM[3] = gmset->get_GM(GammaMatrixSet::GAMMA4);
    m_GM[4] = gmset->get_GM(GammaMatrixSet::GAMMA5);

    m_U = 0;

    m_mult     = &Fopr_Wilson::mult_undef;
    m_mult_dag = &Fopr_Wilson::mult_undef;

    if (m_repr == "Dirac") {
      m_D       = &Fopr_Wilson::D_dirac;
      m_gm5     = &Fopr_Wilson::gm5_dirac;
      m_mult_tp = &Fopr_Wilson::mult_tp_dirac;
      m_mult_tm = &Fopr_Wilson::mult_tm_dirac;
      m_D_ex    = &Fopr_Wilson::D_ex_dirac;
    } else if (m_repr == "Chiral") {
      m_D       = &Fopr_Wilson::D_chiral;
      m_gm5     = &Fopr_Wilson::gm5_chiral;
      m_mult_tp = &Fopr_Wilson::mult_tp_chiral;
      m_mult_tm = &Fopr_Wilson::mult_tm_chiral;
      m_D_ex    = &Fopr_Wilson::D_ex_chiral;
    } else {
      vout.crucial(m_vl, "Error at %s: input repr is undefined.\n", class_name.c_str());
      exit(EXIT_FAILURE);
    }

    m_w1.reset(m_Nvc * m_Nd, m_Nvol, 1);
    m_w2.reset(m_Nvc * m_Nd, m_Nvol, 1);

    int Nvx = m_Nvc * 2 * m_Ny * m_Nz * m_Nt;
    vcp1_xp = new double[Nvx];
    vcp2_xp = new double[Nvx];
    vcp1_xm = new double[Nvx];
    vcp2_xm = new double[Nvx];

    int Nvy = m_Nvc * 2 * m_Nx * m_Nz * m_Nt;
    vcp1_yp = new double[Nvy];
    vcp2_yp = new double[Nvy];
    vcp1_ym = new double[Nvy];
    vcp2_ym = new double[Nvy];

    int Nvz = m_Nvc * 2 * m_Nx * m_Ny * m_Nt;
    vcp1_zp = new double[Nvz];
    vcp2_zp = new double[Nvz];
    vcp1_zm = new double[Nvz];
    vcp2_zm = new double[Nvz];

    int Nvt = m_Nvc * 2 * m_Nx * m_Ny * m_Nz;
    vcp1_tp = new double[Nvt];
    vcp2_tp = new double[Nvt];
    vcp1_tm = new double[Nvt];
    vcp2_tm = new double[Nvt];

    setup_thread();
  }


//====================================================================
  void Fopr_Wilson::set_mode(std::string mode)
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
  void Fopr_Wilson::tidyup()
  {
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
  std::string Fopr_Wilson::get_mode() const
  {
    return m_mode;
  }


//====================================================================
  void Fopr_Wilson::set_parameters(const Parameters& params)
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
      vout.crucial(m_vl, "Error at %s: input parameter not found.\n", class_name.c_str());
      exit(EXIT_FAILURE);
    }

    set_parameters(kappa, bc);
  }


//====================================================================
  void Fopr_Wilson::set_parameters(const double kappa,
                                   const std::vector<int> bc)
  {
    assert(bc.size() == m_Ndim);

    m_kappa = kappa;

    for (int mu = 0; mu < m_Ndim; ++mu) {
      m_boundary[mu] = bc[mu];
    }

    vout.general(m_vl, "%s: input parameters\n", class_name.c_str());
    vout.general(m_vl, "  kappa  = %12.8f\n", m_kappa);
    for (int mu = 0; mu < m_Ndim; ++mu) {
      vout.general(m_vl, "  boundary[%d] = %2d\n", mu, m_boundary[mu]);
    }

    // boundary condition for each node:
    for (int idir = 0; idir < m_Ndim; ++idir) {
      m_boundary2[idir] = 1.0;
      if (Communicator::ipe(idir) == 0) m_boundary2[idir] = m_boundary[idir];
    }
  }


//====================================================================
  double Fopr_Wilson::flop_count()
  {
    // The following counting explicitly depends on the implementation.
    // It will be recalculated when the code is modified.
    // The present counting is based on rev.1107. [24 Aug 2014 H.Matsufuru]

    int    Lvol = CommonParameters::Lvol();
    double flop_site, flop;

    if (m_repr == "Dirac") {
      flop_site = static_cast<double>(
        m_Nc * m_Nd * (4 + 6 * (4 * m_Nc + 2) + 2 * (4 * m_Nc + 1)));
    } else if (m_repr == "Chiral") {
      flop_site = static_cast<double>(
        m_Nc * m_Nd * (4 + 8 * (4 * m_Nc + 2)));
    } else {
      vout.crucial(m_vl, "Error at %s: input repr is undefined.\n", class_name.c_str());
      exit(EXIT_FAILURE);
    }

    flop = flop_site * static_cast<double>(Lvol);

    if ((m_mode == "DdagD") || (m_mode == "DDdag")) flop *= 2.0;

    return flop;
  }


//====================================================================
  void Fopr_Wilson::D_dirac(Field& w, const Field& f)
  {
    D_ex_dirac(w, 0, f, 0);
  }


//====================================================================
  void Fopr_Wilson::D_chiral(Field& w, const Field& f)
  {
    D_ex_chiral(w, 0, f, 0);
  }


//====================================================================
  void Fopr_Wilson::mult_up(int mu,
                            Field& w, const Field& f)
  {
    if (mu == 0) {
      mult_xp(w, f);
    } else if (mu == 1) {
      mult_yp(w, f);
    } else if (mu == 2) {
      mult_zp(w, f);
    } else if (mu == 3) {
      (this->*m_mult_tp)(w, f);
    } else {
      vout.crucial(m_vl, "Error at %s::mult_up: illegal mu=%d\n", class_name.c_str(), mu);
      exit(EXIT_FAILURE);
    }
  }


//====================================================================
  void Fopr_Wilson::mult_dn(int mu, Field& w, const Field& f)
  {
    if (mu == 0) {
      mult_xm(w, f);
    } else if (mu == 1) {
      mult_ym(w, f);
    } else if (mu == 2) {
      mult_zm(w, f);
    } else if (mu == 3) {
      (this->*m_mult_tm)(w, f);
    } else {
      vout.crucial(m_vl, "Error at %s::mult_dn: illegal mu=%d\n", class_name.c_str(), mu);
      exit(EXIT_FAILURE);
    }
  }


//====================================================================
  void Fopr_Wilson::D_ex_dirac(Field& w, const int ex1,
                               const Field& f, const int ex2)
  {
    int          Ninvol = m_Nvc * m_Nd * m_Nvol;
    const double *v1    = f.ptr(Ninvol * ex2);
    double       *v2    = w.ptr(Ninvol * ex1);

    int Nthread  = ThreadManager_OpenMP::get_num_threads();
    int i_thread = ThreadManager_OpenMP::get_thread_id();

    int is = m_Ntask * i_thread / Nthread;
    int ns = m_Ntask * (i_thread + 1) / Nthread - is;

    for (int i = is; i < is + ns; ++i) {
      mult_xp1_thread(i, vcp1_xp, v1);
      mult_xm1_thread(i, vcp1_xm, v1);
      mult_yp1_thread(i, vcp1_yp, v1);
      mult_ym1_thread(i, vcp1_ym, v1);
      mult_zp1_thread(i, vcp1_zp, v1);
      mult_zm1_thread(i, vcp1_zm, v1);
      mult_tp1_dirac_thread(i, vcp1_tp, v1);
      mult_tm1_dirac_thread(i, vcp1_tm, v1);
    }

#pragma omp barrier
#pragma omp master
    {
      int Nvx = m_Nvc * 2 * m_Ny * m_Nz * m_Nt;
      Communicator::exchange(Nvx, vcp2_xp, vcp1_xp, 0, 1, 1);
      Communicator::exchange(Nvx, vcp2_xm, vcp1_xm, 0, -1, 2);

      int Nvy = m_Nvc * 2 * m_Nx * m_Nz * m_Nt;
      Communicator::exchange(Nvy, vcp2_yp, vcp1_yp, 1, 1, 3);
      Communicator::exchange(Nvy, vcp2_ym, vcp1_ym, 1, -1, 4);

      int Nvz = m_Nvc * 2 * m_Nx * m_Ny * m_Nt;
      Communicator::exchange(Nvz, vcp2_zp, vcp1_zp, 2, 1, 5);
      Communicator::exchange(Nvz, vcp2_zm, vcp1_zm, 2, -1, 6);

      int Nvt = m_Nvc * 2 * m_Nx * m_Ny * m_Nz;
      Communicator::exchange(Nvt, vcp2_tp, vcp1_tp, 3, 1, 7);
      Communicator::exchange(Nvt, vcp2_tm, vcp1_tm, 3, -1, 8);
    }
#pragma omp barrier

    for (int i = is; i < is + ns; ++i) {
      clear_thread(i, v2);
      mult_xpb_thread(i, v2, v1);
      mult_xmb_thread(i, v2, v1);
      mult_ypb_thread(i, v2, v1);
      mult_ymb_thread(i, v2, v1);
      mult_zpb_thread(i, v2, v1);
      mult_zmb_thread(i, v2, v1);
      mult_tpb_dirac_thread(i, v2, v1);
      mult_tmb_dirac_thread(i, v2, v1);
    }

    for (int i = is; i < is + ns; ++i) {
      mult_xp2_thread(i, v2, vcp2_xp);
      mult_xm2_thread(i, v2, vcp2_xm);
      mult_yp2_thread(i, v2, vcp2_yp);
      mult_ym2_thread(i, v2, vcp2_ym);
      mult_zp2_thread(i, v2, vcp2_zp);
      mult_zm2_thread(i, v2, vcp2_zm);
      mult_tp2_dirac_thread(i, v2, vcp2_tp);
      mult_tm2_dirac_thread(i, v2, vcp2_tm);
      daypx_thread(i, v2, -m_kappa, v1); // w = -m_kappa * w + f.
    }

    ThreadManager_OpenMP::sync_barrier_all();
  }


//====================================================================
  void Fopr_Wilson::D_ex_chiral(Field& w, const int ex1,
                                const Field& f, const int ex2)
  {
    int          Ninvol = m_Nvc * m_Nd * m_Nvol;
    const double *v1    = f.ptr(Ninvol * ex2);
    double       *v2    = w.ptr(Ninvol * ex1);

    int Nthread  = ThreadManager_OpenMP::get_num_threads();
    int i_thread = ThreadManager_OpenMP::get_thread_id();

    int is = m_Ntask * i_thread / Nthread;
    int ns = m_Ntask * (i_thread + 1) / Nthread - is;

    for (int i = is; i < is + ns; ++i) {
      mult_xp1_thread(i, vcp1_xp, v1);
      mult_xm1_thread(i, vcp1_xm, v1);
      mult_yp1_thread(i, vcp1_yp, v1);
      mult_ym1_thread(i, vcp1_ym, v1);
      mult_zp1_thread(i, vcp1_zp, v1);
      mult_zm1_thread(i, vcp1_zm, v1);
      mult_tp1_chiral_thread(i, vcp1_tp, v1);
      mult_tm1_chiral_thread(i, vcp1_tm, v1);
    }
#pragma omp barrier

#pragma omp master
    {
      int Nvx = m_Nvc * 2 * m_Ny * m_Nz * m_Nt;
      Communicator::exchange(Nvx, vcp2_xp, vcp1_xp, 0, 1, 1);
      Communicator::exchange(Nvx, vcp2_xm, vcp1_xm, 0, -1, 2);

      int Nvy = m_Nvc * 2 * m_Nx * m_Nz * m_Nt;
      Communicator::exchange(Nvy, vcp2_yp, vcp1_yp, 1, 1, 3);
      Communicator::exchange(Nvy, vcp2_ym, vcp1_ym, 1, -1, 4);

      int Nvz = m_Nvc * 2 * m_Nx * m_Ny * m_Nt;
      Communicator::exchange(Nvz, vcp2_zp, vcp1_zp, 2, 1, 5);
      Communicator::exchange(Nvz, vcp2_zm, vcp1_zm, 2, -1, 6);

      int Nvt = m_Nvc * 2 * m_Nx * m_Ny * m_Nz;
      Communicator::exchange(Nvt, vcp2_tp, vcp1_tp, 3, 1, 7);
      Communicator::exchange(Nvt, vcp2_tm, vcp1_tm, 3, -1, 8);
    }
#pragma omp barrier

    for (int i = is; i < is + ns; ++i) {
      clear_thread(i, v2);
      mult_xpb_thread(i, v2, v1);
      mult_xmb_thread(i, v2, v1);
      mult_ypb_thread(i, v2, v1);
      mult_ymb_thread(i, v2, v1);
      mult_zpb_thread(i, v2, v1);
      mult_zmb_thread(i, v2, v1);
      mult_tpb_chiral_thread(i, v2, v1);
      mult_tmb_chiral_thread(i, v2, v1);
    }

    for (int i = is; i < is + ns; ++i) {
      mult_xp2_thread(i, v2, vcp2_xp);
      mult_xm2_thread(i, v2, vcp2_xm);
      mult_yp2_thread(i, v2, vcp2_yp);
      mult_ym2_thread(i, v2, vcp2_ym);
      mult_zp2_thread(i, v2, vcp2_zp);
      mult_zm2_thread(i, v2, vcp2_zm);
      mult_tp2_chiral_thread(i, v2, vcp2_tp);
      mult_tm2_chiral_thread(i, v2, vcp2_tm);
      daypx_thread(i, v2, -m_kappa, v1); // w = -m_kappa * w + f.
    }

    ThreadManager_OpenMP::sync_barrier_all();

    /*
     clear(w);
     mult_xp(w,f);
     mult_xm(w,f);
     mult_yp(w,f);
     mult_ym(w,f);
     mult_zp(w,f);
     mult_zm(w,f);
     mult_tp_chiral(w,f);
     mult_tm_chiral(w,f);
     daypx(w, -m_kappa, f);   // w = -m_kappa * w + f.
    */
  }


//====================================================================
  void Fopr_Wilson::mult_gm5p(int mu,
                              Field& v, const Field& w)
  {
    clear(m_w2);

    mult_up(mu, m_w2, w);
    mult_gm5(v, m_w2);
  }


//====================================================================
  void Fopr_Wilson::proj_chiral(Field& w, const int ex1,
                                const Field& v, const int ex2, const int ipm)
  {
    double fpm = 0.0;

    if (ipm == 1) {
      fpm = 1.0;
    } else if (ipm == -1) {
      fpm = -1.0;
    } else {
      vout.crucial(m_vl, "Error at %s: illegal chirality = %d\n", class_name.c_str(), ipm);
      exit(EXIT_FAILURE);
    }

    m_w1.setpart_ex(0, v, ex2);
    mult_gm5(m_w2, m_w1);
    m_w1.addpart_ex(0, m_w2, 0, fpm);
    w.setpart_ex(ex1, m_w1, 0);

#pragma omp barrier
  }


//====================================================================
  void Fopr_Wilson::daypx(Field& w,
                          double fac, const Field& f)
  {
    const double *v1 = f.ptr(0);
    double       *v2 = w.ptr(0);

    int Nthread  = ThreadManager_OpenMP::get_num_threads();
    int i_thread = ThreadManager_OpenMP::get_thread_id();

    int is = m_Ntask * i_thread / Nthread;
    int ns = m_Ntask * (i_thread + 1) / Nthread - is;

    for (int i = is; i < is + ns; ++i) {
      daypx_thread(i, v2, fac, v1);
    }
#pragma omp barrier
  }


//====================================================================
  void Fopr_Wilson::clear(Field& w)
  {
    double *v2 = w.ptr(0);

    int Nthread  = ThreadManager_OpenMP::get_num_threads();
    int i_thread = ThreadManager_OpenMP::get_thread_id();

    int is = m_Ntask * i_thread / Nthread;
    int ns = m_Ntask * (i_thread + 1) / Nthread - is;

    for (int i = is; i < is + ns; ++i) {
      clear_thread(i, v2);
    }
#pragma omp barrier
  }


//====================================================================
  void Fopr_Wilson::gm5_dirac(Field& w, const Field& f)
  {
    const double *v1 = f.ptr(0);
    double       *v2 = w.ptr(0);

    int Nthread  = ThreadManager_OpenMP::get_num_threads();
    int i_thread = ThreadManager_OpenMP::get_thread_id();

    int is = m_Ntask * i_thread / Nthread;
    int ns = m_Ntask * (i_thread + 1) / Nthread - is;

    for (int i = is; i < is + ns; ++i) {
      gm5_dirac_thread(i, v2, v1);
    }
#pragma omp barrier
  }


//====================================================================
  void Fopr_Wilson::gm5_chiral(Field& w, const Field& f)
  {
    const double *v1 = f.ptr(0);
    double       *v2 = w.ptr(0);

    int Nthread  = ThreadManager_OpenMP::get_num_threads();
    int i_thread = ThreadManager_OpenMP::get_thread_id();

    int is = m_Ntask * i_thread / Nthread;
    int ns = m_Ntask * (i_thread + 1) / Nthread - is;

    for (int i = is; i < is + ns; ++i) {
      gm5_chiral_thread(i, v2, v1);
    }
#pragma omp barrier
  }


//====================================================================
  void Fopr_Wilson::mult_xp(Field& w, const Field& f)
  {
    const double *v1 = f.ptr(0);
    double       *v2 = w.ptr(0);

    int Nthread  = ThreadManager_OpenMP::get_num_threads();
    int i_thread = ThreadManager_OpenMP::get_thread_id();

    int is = m_Ntask * i_thread / Nthread;
    int ns = m_Ntask * (i_thread + 1) / Nthread - is;

#pragma omp barrier

    for (int i = is; i < is + ns; ++i) {
      mult_xp1_thread(i, vcp1_xp, v1);
    }

#pragma omp barrier
#pragma omp master
    {
      int Nv = m_Nvc * 2 * m_Ny * m_Nz * m_Nt;
      Communicator::exchange(Nv, vcp2_xp, vcp1_xp, 0, 1, 1);
    }
#pragma omp barrier

    for (int i = is; i < is + ns; ++i) {
      mult_xpb_thread(i, v2, v1);
    }

    for (int i = is; i < is + ns; ++i) {
      mult_xp2_thread(i, v2, vcp2_xp);
    }

    ThreadManager_OpenMP::sync_barrier_all();
  }


//====================================================================
  void Fopr_Wilson::mult_xm(Field& w, const Field& f)
  {
    double       *v2 = w.ptr(0);
    const double *v1 = f.ptr(0);

    int Nthread  = ThreadManager_OpenMP::get_num_threads();
    int i_thread = ThreadManager_OpenMP::get_thread_id();

    int is = m_Ntask * i_thread / Nthread;
    int ns = m_Ntask * (i_thread + 1) / Nthread - is;

#pragma omp barrier

    for (int i = is; i < is + ns; ++i) {
      mult_xm1_thread(i, vcp1_xm, v1);
    }

#pragma omp barrier
#pragma omp master
    {
      int Nv = m_Nvc * 2 * m_Ny * m_Nz * m_Nt;
      Communicator::exchange(Nv, vcp2_xm, vcp1_xm, 0, -1, 2);
    }
#pragma omp barrier

    for (int i = is; i < is + ns; ++i) {
      mult_xmb_thread(i, v2, v1);
    }

    for (int i = is; i < is + ns; ++i) {
      mult_xm2_thread(i, v2, vcp2_xm);
    }

    ThreadManager_OpenMP::sync_barrier_all();
  }


//====================================================================
  void Fopr_Wilson::mult_yp(Field& w, const Field& f)
  {
    const double *v1 = f.ptr(0);
    double       *v2 = w.ptr(0);

    int Nthread  = ThreadManager_OpenMP::get_num_threads();
    int i_thread = ThreadManager_OpenMP::get_thread_id();

    int is = m_Ntask * i_thread / Nthread;
    int ns = m_Ntask * (i_thread + 1) / Nthread - is;

#pragma omp barrier

    for (int i = is; i < is + ns; ++i) {
      mult_yp1_thread(i, vcp1_yp, v1);
    }

#pragma omp barrier
#pragma omp master
    {
      int Nv = m_Nvc * 2 * m_Nx * m_Nz * m_Nt;
      Communicator::exchange(Nv, vcp2_yp, vcp1_yp, 1, 1, 3);
    }
#pragma omp barrier

    for (int i = is; i < is + ns; ++i) {
      mult_ypb_thread(i, v2, v1);
    }

    for (int i = is; i < is + ns; ++i) {
      mult_yp2_thread(i, v2, vcp2_yp);
    }

    ThreadManager_OpenMP::sync_barrier_all();
  }


//====================================================================
  void Fopr_Wilson::mult_ym(Field& w, const Field& f)
  {
    double       *v2 = w.ptr(0);
    const double *v1 = f.ptr(0);

    int Nthread  = ThreadManager_OpenMP::get_num_threads();
    int i_thread = ThreadManager_OpenMP::get_thread_id();

    int is = m_Ntask * i_thread / Nthread;
    int ns = m_Ntask * (i_thread + 1) / Nthread - is;

#pragma omp barrier

    for (int i = is; i < is + ns; ++i) {
      mult_ym1_thread(i, vcp1_ym, v1);
    }

#pragma omp barrier
#pragma omp master
    {
      int Nv = m_Nvc * 2 * m_Nx * m_Nz * m_Nt;
      Communicator::exchange(Nv, vcp2_ym, vcp1_ym, 1, -1, 4);
    }
#pragma omp barrier

    for (int i = is; i < is + ns; ++i) {
      mult_ymb_thread(i, v2, v1);
    }

    for (int i = is; i < is + ns; ++i) {
      mult_ym2_thread(i, v2, vcp2_ym);
    }

    ThreadManager_OpenMP::sync_barrier_all();
  }


//====================================================================
  void Fopr_Wilson::mult_zp(Field& w, const Field& f)
  {
    const double *v1 = f.ptr(0);
    double       *v2 = w.ptr(0);

    int Nthread  = ThreadManager_OpenMP::get_num_threads();
    int i_thread = ThreadManager_OpenMP::get_thread_id();

    int is = m_Ntask * i_thread / Nthread;
    int ns = m_Ntask * (i_thread + 1) / Nthread - is;

#pragma omp barrier

    for (int i = is; i < is + ns; ++i) {
      mult_zp1_thread(i, vcp1_zp, v1);
    }

#pragma omp barrier
#pragma omp master
    {
      int Nv = m_Nvc * 2 * m_Nx * m_Ny * m_Nt;
      Communicator::exchange(Nv, vcp2_zp, vcp1_zp, 2, 1, 5);
    }
#pragma omp barrier

    for (int i = is; i < is + ns; ++i) {
      mult_zpb_thread(i, v2, v1);
    }

    for (int i = is; i < is + ns; ++i) {
      mult_zp2_thread(i, v2, vcp2_zp);
    }

    ThreadManager_OpenMP::sync_barrier_all();
  }


//====================================================================
  void Fopr_Wilson::mult_zm(Field& w, const Field& f)
  {
    double       *v2 = w.ptr(0);
    const double *v1 = f.ptr(0);

    int Nthread  = ThreadManager_OpenMP::get_num_threads();
    int i_thread = ThreadManager_OpenMP::get_thread_id();

    int is = m_Ntask * i_thread / Nthread;
    int ns = m_Ntask * (i_thread + 1) / Nthread - is;

#pragma omp barrier

    for (int i = is; i < is + ns; ++i) {
      mult_zm1_thread(i, vcp1_zm, v1);
    }

#pragma omp barrier
#pragma omp master
    {
      int Nv = m_Nvc * 2 * m_Nx * m_Ny * m_Nt;
      Communicator::exchange(Nv, vcp2_zm, vcp1_zm, 2, -1, 6);
    }
#pragma omp barrier

    for (int i = is; i < is + ns; ++i) {
      mult_zmb_thread(i, v2, v1);
    }

    for (int i = is; i < is + ns; ++i) {
      mult_zm2_thread(i, v2, vcp2_zm);
    }

    ThreadManager_OpenMP::sync_barrier_all();
  }


//====================================================================
  void Fopr_Wilson::mult_tp_dirac(Field& w, const Field& f)
  {
    const double *v1 = f.ptr(0);
    double       *v2 = w.ptr(0);

    int Nthread  = ThreadManager_OpenMP::get_num_threads();
    int i_thread = ThreadManager_OpenMP::get_thread_id();

    int is = m_Ntask * i_thread / Nthread;
    int ns = m_Ntask * (i_thread + 1) / Nthread - is;

#pragma omp barrier

    for (int i = is; i < is + ns; ++i) {
      mult_tp1_dirac_thread(i, vcp1_tp, v1);
    }

#pragma omp barrier
#pragma omp master
    {
      int Nv = m_Nvc * 2 * m_Nx * m_Ny * m_Nz;
      Communicator::exchange(Nv, vcp2_tp, vcp1_tp, 3, 1, 7);
    }
#pragma omp barrier

    for (int i = is; i < is + ns; ++i) {
      mult_tpb_dirac_thread(i, v2, v1);
    }

    for (int i = is; i < is + ns; ++i) {
      mult_tp2_dirac_thread(i, v2, vcp2_tp);
    }

    ThreadManager_OpenMP::sync_barrier_all();
  }


//====================================================================
  void Fopr_Wilson::mult_tm_dirac(Field& w, const Field& f)
  {
    double       *v2 = w.ptr(0);
    const double *v1 = f.ptr(0);

    int Nthread  = ThreadManager_OpenMP::get_num_threads();
    int i_thread = ThreadManager_OpenMP::get_thread_id();

    int is = m_Ntask * i_thread / Nthread;
    int ns = m_Ntask * (i_thread + 1) / Nthread - is;

#pragma omp barrier

    for (int i = is; i < is + ns; ++i) {
      mult_tm1_dirac_thread(i, vcp1_tm, v1);
    }

#pragma omp barrier
#pragma omp master
    {
      int Nv = m_Nvc * 2 * m_Nx * m_Ny * m_Nz;
      Communicator::exchange(Nv, vcp2_tm, vcp1_tm, 3, -1, 8);
    }
#pragma omp barrier

    for (int i = is; i < is + ns; ++i) {
      mult_tmb_dirac_thread(i, v2, v1);
    }

    for (int i = is; i < is + ns; ++i) {
      mult_tm2_dirac_thread(i, v2, vcp2_tm);
    }

    ThreadManager_OpenMP::sync_barrier_all();
  }


//====================================================================
  void Fopr_Wilson::mult_tp_chiral(Field& w, const Field& f)
  {
    const double *v1 = f.ptr(0);
    double       *v2 = w.ptr(0);

    int Nthread  = ThreadManager_OpenMP::get_num_threads();
    int i_thread = ThreadManager_OpenMP::get_thread_id();

    int is = m_Ntask * i_thread / Nthread;
    int ns = m_Ntask * (i_thread + 1) / Nthread - is;

#pragma omp barrier

    for (int i = is; i < is + ns; ++i) {
      mult_tp1_chiral_thread(i, vcp1_tp, v1);
    }

#pragma omp barrier
#pragma omp master
    {
      int Nv = m_Nvc * 2 * m_Nx * m_Ny * m_Nz;
      Communicator::exchange(Nv, vcp2_tp, vcp1_tp, 3, 1, 7);
    }
#pragma omp barrier

    for (int i = is; i < is + ns; ++i) {
      mult_tpb_chiral_thread(i, v2, v1);
    }

    for (int i = is; i < is + ns; ++i) {
      mult_tp2_chiral_thread(i, v2, vcp2_tp);
    }

    ThreadManager_OpenMP::sync_barrier_all();
  }


//====================================================================
  void Fopr_Wilson::mult_tm_chiral(Field& w, const Field& f)
  {
    double       *v2 = w.ptr(0);
    const double *v1 = f.ptr(0);

    int Nthread  = ThreadManager_OpenMP::get_num_threads();
    int i_thread = ThreadManager_OpenMP::get_thread_id();

    int is = m_Ntask * i_thread / Nthread;
    int ns = m_Ntask * (i_thread + 1) / Nthread - is;

#pragma omp barrier

    for (int i = is; i < is + ns; ++i) {
      mult_tm1_chiral_thread(i, vcp1_tm, v1);
    }

#pragma omp barrier
#pragma omp master
    {
      int Nv = m_Nvc * 2 * m_Nx * m_Ny * m_Nz;
      Communicator::exchange(Nv, vcp2_tm, vcp1_tm, 3, -1, 8);
    }
#pragma omp barrier

    for (int i = is; i < is + ns; ++i) {
      mult_tmb_chiral_thread(i, v2, v1);
    }

    for (int i = is; i < is + ns; ++i) {
      mult_tm2_chiral_thread(i, v2, vcp2_tm);
    }

    ThreadManager_OpenMP::sync_barrier_all();
  }


//====================================================================
}
//============================================================END=====

#endif
