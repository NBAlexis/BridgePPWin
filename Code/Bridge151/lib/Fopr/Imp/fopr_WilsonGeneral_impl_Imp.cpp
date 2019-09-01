/*!
        @file    fopr_WilsonGeneral_impl.cpp

        @brief

        @author  Yusuke Namekawa (namekawa)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/
#include "BridgeLib_Private.h"
#if USE_IMP
#include "fopr_WilsonGeneral_impl.h"

namespace Imp {
#if defined USE_GROUP_SU3
#include "fopr_Wilson_impl_SU3-inc.h"
#elif defined USE_GROUP_SU2
#include "fopr_Wilson_impl_SU2-inc.h"
#elif defined USE_GROUP_SU_N
#include "fopr_Wilson_impl_SU_N-inc.h"
#endif

  const std::string Fopr_WilsonGeneral::class_name = "Imp::Fopr_WilsonGeneral";

#ifdef USE_FACTORY_AUTOREGISTER
  namespace {
    bool init = Fopr_WilsonGeneral::register_factory();
  }
#endif

//====================================================================
  void Fopr_WilsonGeneral::init(const std::string repr)
  {
    m_vl = CommonParameters::Vlevel();

    vout.general(m_vl, "%s: Construction of fermion operator(imp).\n", class_name.c_str());

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
    m_boundary_each_node.resize(m_Ndim);

    m_repr = repr;

    m_GM.resize(m_Ndim + 1);

    unique_ptr<GammaMatrixSet> gmset(GammaMatrixSet::New(m_repr));

    m_GM[0] = gmset->get_GM(GammaMatrixSet::GAMMA1);
    m_GM[1] = gmset->get_GM(GammaMatrixSet::GAMMA2);
    m_GM[2] = gmset->get_GM(GammaMatrixSet::GAMMA3);
    m_GM[3] = gmset->get_GM(GammaMatrixSet::GAMMA4);
    m_GM[4] = gmset->get_GM(GammaMatrixSet::GAMMA5);

    m_U = 0;

    m_mult     = &Fopr_WilsonGeneral::mult_undef;
    m_mult_dag = &Fopr_WilsonGeneral::mult_undef;

    if (m_repr == "Dirac") {
      m_D            = &Fopr_WilsonGeneral::D_dirac;
      m_gm5          = &Fopr_WilsonGeneral::gm5_dirac;
      m_mult_t_plus  = &Fopr_WilsonGeneral::mult_t_plus_dirac;
      m_mult_t_minus = &Fopr_WilsonGeneral::mult_t_minus_dirac;
      m_D_ex         = &Fopr_WilsonGeneral::D_ex_dirac;
    } else if (m_repr == "Chiral") {
      m_D            = &Fopr_WilsonGeneral::D_chiral;
      m_gm5          = &Fopr_WilsonGeneral::gm5_chiral;
      m_mult_t_plus  = &Fopr_WilsonGeneral::mult_t_plus_chiral;
      m_mult_t_minus = &Fopr_WilsonGeneral::mult_t_minus_chiral;
      m_D_ex         = &Fopr_WilsonGeneral::D_ex_chiral;
    } else {
      vout.crucial(m_vl, "Error at %s: input repr is undefined.\n", class_name.c_str());
      exit(EXIT_FAILURE);
    }

    m_w1.reset(m_Nvc * m_Nd, m_Nvol, 1);
    m_w2.reset(m_Nvc * m_Nd, m_Nvol, 1);

    const int Nvx = 2 * m_Nvc * m_Nd * m_Ny * m_Nz * m_Nt;
    vcp1_x_plus  = new double[Nvx];
    vcp2_x_plus  = new double[Nvx];
    vcp1_x_minus = new double[Nvx];
    vcp2_x_minus = new double[Nvx];

    const int Nvy = 2 * m_Nvc * m_Nd * m_Nx * m_Nz * m_Nt;
    vcp1_y_plus  = new double[Nvy];
    vcp2_y_plus  = new double[Nvy];
    vcp1_y_minus = new double[Nvy];
    vcp2_y_minus = new double[Nvy];

    const int Nvz = 2 * m_Nvc * m_Nd * m_Nx * m_Ny * m_Nt;
    vcp1_z_plus  = new double[Nvz];
    vcp2_z_plus  = new double[Nvz];
    vcp1_z_minus = new double[Nvz];
    vcp2_z_minus = new double[Nvz];

    // const int Nvt = 2 * m_Nvc * m_Nd * m_Nx * m_Ny * m_Nz;
    const int Nvt = m_Nvc * m_Nd * m_Nx * m_Ny * m_Nz;
    vcp1_t_plus  = new double[Nvt];
    vcp2_t_plus  = new double[Nvt];
    vcp1_t_minus = new double[Nvt];
    vcp2_t_minus = new double[Nvt];

    setup_thread();
  }


//====================================================================
  void Fopr_WilsonGeneral::set_mode(const std::string mode)
  {
    m_mode = mode;

    if (m_mode == "D") {
      m_mult     = &Fopr_WilsonGeneral::D;
      m_mult_dag = &Fopr_WilsonGeneral::Ddag;
    } else if (m_mode == "Ddag") {
      m_mult     = &Fopr_WilsonGeneral::Ddag;
      m_mult_dag = &Fopr_WilsonGeneral::D;
    } else if (m_mode == "DdagD") {
      m_mult     = &Fopr_WilsonGeneral::DdagD;
      m_mult_dag = &Fopr_WilsonGeneral::DdagD;
    } else if (m_mode == "DDdag") {
      m_mult     = &Fopr_WilsonGeneral::DDdag;
      m_mult_dag = &Fopr_WilsonGeneral::DDdag;
    } else if (m_mode == "H") {
      m_mult     = &Fopr_WilsonGeneral::H;
      m_mult_dag = &Fopr_WilsonGeneral::H;
    } else {
      vout.crucial(m_vl, "Error at %s: input mode is undefined.\n", class_name.c_str());
      exit(EXIT_FAILURE);
    }
  }


//====================================================================
  void Fopr_WilsonGeneral::tidyup()
  {
    delete[]  vcp1_x_plus;
    delete[]  vcp2_x_plus;
    delete[]  vcp1_x_minus;
    delete[]  vcp2_x_minus;

    delete[]  vcp1_y_plus;
    delete[]  vcp2_y_plus;
    delete[]  vcp1_y_minus;
    delete[]  vcp2_y_minus;

    delete[]  vcp1_z_plus;
    delete[]  vcp2_z_plus;
    delete[]  vcp1_z_minus;
    delete[]  vcp2_z_minus;

    delete[]  vcp1_t_plus;
    delete[]  vcp2_t_plus;
    delete[]  vcp1_t_minus;
    delete[]  vcp2_t_minus;
  }


//====================================================================
  std::string Fopr_WilsonGeneral::get_mode() const
  {
    return m_mode;
  }


//====================================================================
  void Fopr_WilsonGeneral::set_parameters(const Parameters& params)
  {
    const string str_vlevel = params.get_string("verbose_level");

    m_vl = vout.set_verbose_level(str_vlevel);

    //- fetch and check input parameters
    double           kappa_s, kappa_t;
    double           nu_s, r_s;
    std::vector<int> bc;

    int err = 0;
    err += params.fetch_double("hopping_parameter_spatial", kappa_s);
    err += params.fetch_double("hopping_parameter_temporal", kappa_t);
    err += params.fetch_double("dispersion_parameter_spatial", nu_s);
    err += params.fetch_double("Wilson_parameter_spatial", r_s);
    err += params.fetch_int_vector("boundary_condition", bc);

    if (err) {
      vout.crucial(m_vl, "Error at %s: input parameter not found.\n", class_name.c_str());
      exit(EXIT_FAILURE);
    }

    set_parameters(kappa_s, kappa_t, nu_s, r_s, bc);
  }


//====================================================================
  void Fopr_WilsonGeneral::set_parameters(const double kappa_s,
                                          const double kappa_t,
                                          const double nu_s,
                                          const double r_s,
                                          const std::vector<int> bc)
  {
    //- print input parameters
    vout.general(m_vl, "%s: input parameters\n", class_name.c_str());
    vout.general(m_vl, "  kappa_s  = %12.8f\n", kappa_s);
    vout.general(m_vl, "  kappa_t  = %12.8f\n", kappa_t);
    vout.general(m_vl, "  nu_s     = %12.8f\n", nu_s);
    vout.general(m_vl, "  r_s      = %12.8f\n", r_s);
    for (int mu = 0; mu < m_Ndim; ++mu) {
      vout.general(m_vl, "  boundary[%d] = %2d\n", mu, bc[mu]);
    }

    //- range check
    // NB. kappa,cSW == 0 is allowed.
    assert(bc.size() == m_Ndim);

    //- store values
    m_kappa_s = kappa_s;
    m_kappa_t = kappa_t;
    m_nu_s    = nu_s;
    m_r_s     = r_s;

    // m_boundary.resize(m_Ndim);  // already resized in init.
    m_boundary = bc;

    for (int idir = 0; idir < m_Ndim; ++idir) {
      m_boundary_each_node[idir] = 1.0;
      if (Communicator::ipe(idir) == 0) m_boundary_each_node[idir] = m_boundary[idir];
    }
  }


//====================================================================
  void Fopr_WilsonGeneral::D_dirac(Field& w, const Field& f)
  {
    // D_ex_dirac(w, 0, f, 0);
    // daypx(w, -m_kappa_s, f);

    clear(w);

    clear(m_w1);
    for (int mu = 0; mu < (m_Ndim - 1); ++mu) {
      mult_up(mu, m_w1, f);
      mult_dn(mu, m_w1, f);
    }
    daxpy(w, -m_kappa_s, m_w1);

    clear(m_w1);
    const int mu = (m_Ndim - 1);
    {
      mult_up(mu, m_w1, f);
      mult_dn(mu, m_w1, f);
    }
    daxpy(w, -m_kappa_t, m_w1);


    daxpy(w, 1.0, f); // w += f;
  }


//====================================================================
  void Fopr_WilsonGeneral::D_chiral(Field& w, const Field& f)
  {
    // D_ex_chiral(w, 0, f, 0);
    // daypx(w, -m_kappa_s, f);

    clear(w);

    clear(m_w1);
    for (int mu = 0; mu < (m_Ndim - 1); ++mu) {
      mult_up(mu, m_w1, f);
      mult_dn(mu, m_w1, f);
    }
    daxpy(w, -m_kappa_s, m_w1);

    clear(m_w1);
    const int mu = (m_Ndim - 1);
    {
      mult_up(mu, m_w1, f);
      mult_dn(mu, m_w1, f);
    }
    daxpy(w, -m_kappa_t, m_w1);


    daxpy(w, 1.0, f); // w += f;
  }


//====================================================================
  void Fopr_WilsonGeneral::mult_up(const int mu,
                                   Field& w, const Field& f)
  {
    if (mu == 0) {
      mult_x_plus(w, f);
    } else if (mu == 1) {
      mult_y_plus(w, f);
    } else if (mu == 2) {
      mult_z_plus(w, f);
    } else if (mu == 3) {
      (this->*m_mult_t_plus)(w, f);
    } else {
      vout.crucial(m_vl, "Error at %s: illegal mu=%d in mult_up.\n", class_name.c_str(), mu);
      exit(EXIT_FAILURE);
    }
  }


//====================================================================
  void Fopr_WilsonGeneral::mult_dn(const int mu,
                                   Field& w, const Field& f)
  {
    if (mu == 0) {
      mult_x_minus(w, f);
    } else if (mu == 1) {
      mult_y_minus(w, f);
    } else if (mu == 2) {
      mult_z_minus(w, f);
    } else if (mu == 3) {
      (this->*m_mult_t_minus)(w, f);
    } else {
      vout.crucial(m_vl, "Error at %s: illegal mu=%d in mult_dn.\n", class_name.c_str(), mu);
      exit(EXIT_FAILURE);
    }
  }


//====================================================================
  void Fopr_WilsonGeneral::D_ex_dirac(Field& w, const int ex1,
                                      const Field& f, const int ex2)
  {
    const int    Ninvol = m_Nvc * m_Nd * m_Nvol;
    const double *v1    = f.ptr(Ninvol * ex2);
    double       *v2    = w.ptr(Ninvol * ex1);

    const int Nthread  = ThreadManager_OpenMP::get_num_threads();
    const int i_thread = ThreadManager_OpenMP::get_thread_id();

    const int is = m_Ntask * i_thread / Nthread;
    const int ns = m_Ntask * (i_thread + 1) / Nthread - is;

    for (int i = is; i < is + ns; ++i) {
      mult_x_plus1_thread(i, vcp1_x_plus, v1);
      mult_x_minus1_thread(i, vcp1_x_minus, v1);
      mult_y_plus1_thread(i, vcp1_y_plus, v1);
      mult_y_minus1_thread(i, vcp1_y_minus, v1);
      mult_z_plus1_thread(i, vcp1_z_plus, v1);
      mult_z_minus1_thread(i, vcp1_z_minus, v1);
      mult_t_plus1_dirac_thread(i, vcp1_t_plus, v1);
      mult_t_minus1_dirac_thread(i, vcp1_t_minus, v1);
    }

#pragma omp barrier
#pragma omp master
    {
      const int idir_x = 0;
      const int Nvx    = 2 * m_Nvc * m_Nd * m_Ny * m_Nz * m_Nt;
      Communicator::exchange(Nvx, vcp2_x_plus, vcp1_x_plus, idir_x, 1, 1);
      Communicator::exchange(Nvx, vcp2_x_minus, vcp1_x_minus, idir_x, -1, 2);

      const int idir_y = 1;
      const int Nvy    = 2 * m_Nvc * m_Nd * m_Nx * m_Nz * m_Nt;
      Communicator::exchange(Nvy, vcp2_y_plus, vcp1_y_plus, idir_y, 1, 3);
      Communicator::exchange(Nvy, vcp2_y_minus, vcp1_y_minus, idir_y, -1, 4);

      const int idir_z = 2;
      const int Nvz    = 2 * m_Nvc * m_Nd * m_Nx * m_Ny * m_Nt;
      Communicator::exchange(Nvz, vcp2_z_plus, vcp1_z_plus, idir_z, 1, 5);
      Communicator::exchange(Nvz, vcp2_z_minus, vcp1_z_minus, idir_z, -1, 6);

      const int idir_t = 3;
      // const int Nvt = 2 * m_Nvc * m_Nd * m_Nx * m_Ny * m_Nz;
      const int Nvt = m_Nvc * m_Nd * m_Nx * m_Ny * m_Nz;
      Communicator::exchange(Nvt, vcp2_t_plus, vcp1_t_plus, idir_t, 1, 7);
      Communicator::exchange(Nvt, vcp2_t_minus, vcp1_t_minus, idir_t, -1, 8);
    }
#pragma omp barrier

    for (int i = is; i < is + ns; ++i) {
      clear_thread(i, v2);

      mult_x_plus_bulk_thread(i, v2, v1);
      mult_x_minus_bulk_thread(i, v2, v1);
      mult_y_plus_bulk_thread(i, v2, v1);
      mult_y_minus_bulk_thread(i, v2, v1);
      mult_z_plus_bulk_thread(i, v2, v1);
      mult_z_minus_bulk_thread(i, v2, v1);
      mult_t_plus_bulk_dirac_thread(i, v2, v1);
      mult_t_minus_bulk_dirac_thread(i, v2, v1);
    }

    for (int i = is; i < is + ns; ++i) {
      mult_x_plus2_thread(i, v2, vcp2_x_plus);
      mult_x_minus2_thread(i, v2, vcp2_x_minus);
      mult_y_plus2_thread(i, v2, vcp2_y_plus);
      mult_y_minus2_thread(i, v2, vcp2_y_minus);
      mult_z_plus2_thread(i, v2, vcp2_z_plus);
      mult_z_minus2_thread(i, v2, vcp2_z_minus);
      mult_t_plus2_dirac_thread(i, v2, vcp2_t_plus);
      mult_t_minus2_dirac_thread(i, v2, vcp2_t_minus);
    }

    // for (int i = is; i < is + ns; ++i) {
    // daypx_thread(i, v2, -m_kappa, v1);   // w = -m_kappa * w + f.
    //   daypx_thread(i, v2, -m_kappa_s, v1);
    // }


    ThreadManager_OpenMP::sync_barrier_all();
  }


//====================================================================
  void Fopr_WilsonGeneral::D_ex_chiral(Field& w, const int ex1,
                                       const Field& f, const int ex2)
  {
    const int    Ninvol = m_Nvc * m_Nd * m_Nvol;
    const double *v1    = f.ptr(Ninvol * ex2);
    double       *v2    = w.ptr(Ninvol * ex1);

    const int Nthread  = ThreadManager_OpenMP::get_num_threads();
    const int i_thread = ThreadManager_OpenMP::get_thread_id();

    const int is = m_Ntask * i_thread / Nthread;
    const int ns = m_Ntask * (i_thread + 1) / Nthread - is;

    for (int i = is; i < is + ns; ++i) {
      mult_x_plus1_thread(i, vcp1_x_plus, v1);
      mult_x_minus1_thread(i, vcp1_x_minus, v1);
      mult_y_plus1_thread(i, vcp1_y_plus, v1);
      mult_y_minus1_thread(i, vcp1_y_minus, v1);
      mult_z_plus1_thread(i, vcp1_z_plus, v1);
      mult_z_minus1_thread(i, vcp1_z_minus, v1);
      mult_t_plus1_chiral_thread(i, vcp1_t_plus, v1);
      mult_t_minus1_chiral_thread(i, vcp1_t_minus, v1);
    }
#pragma omp barrier

#pragma omp master
    {
      const int idir_x = 0;
      const int Nvx    = 2 * m_Nvc * m_Nd * m_Ny * m_Nz * m_Nt;
      Communicator::exchange(Nvx, vcp2_x_plus, vcp1_x_plus, idir_x, 1, 1);
      Communicator::exchange(Nvx, vcp2_x_minus, vcp1_x_minus, idir_x, -1, 2);

      const int idir_y = 1;
      const int Nvy    = 2 * m_Nvc * m_Nd * m_Nx * m_Nz * m_Nt;
      Communicator::exchange(Nvy, vcp2_y_plus, vcp1_y_plus, idir_y, 1, 3);
      Communicator::exchange(Nvy, vcp2_y_minus, vcp1_y_minus, idir_y, -1, 4);

      const int idir_z = 2;
      const int Nvz    = 2 * m_Nvc * m_Nd * m_Nx * m_Ny * m_Nt;
      Communicator::exchange(Nvz, vcp2_z_plus, vcp1_z_plus, idir_z, 1, 5);
      Communicator::exchange(Nvz, vcp2_z_minus, vcp1_z_minus, idir_z, -1, 6);

      const int idir_t = 3;
      // const int Nvt = 2 * m_Nvc * m_Nd * m_Nx * m_Ny * m_Nz;
      const int Nvt = m_Nvc * m_Nd * m_Nx * m_Ny * m_Nz;
      Communicator::exchange(Nvt, vcp2_t_plus, vcp1_t_plus, idir_t, 1, 7);
      Communicator::exchange(Nvt, vcp2_t_minus, vcp1_t_minus, idir_t, -1, 8);
    }
#pragma omp barrier

    for (int i = is; i < is + ns; ++i) {
      clear_thread(i, v2);

      mult_x_plus_bulk_thread(i, v2, v1);
      mult_x_minus_bulk_thread(i, v2, v1);
      mult_y_plus_bulk_thread(i, v2, v1);
      mult_y_minus_bulk_thread(i, v2, v1);
      mult_z_plus_bulk_thread(i, v2, v1);
      mult_z_minus_bulk_thread(i, v2, v1);
      mult_t_plus_bulk_chiral_thread(i, v2, v1);
      mult_t_minus_bulk_chiral_thread(i, v2, v1);
    }

    for (int i = is; i < is + ns; ++i) {
      mult_x_plus2_thread(i, v2, vcp2_x_plus);
      mult_x_minus2_thread(i, v2, vcp2_x_minus);
      mult_y_plus2_thread(i, v2, vcp2_y_plus);
      mult_y_minus2_thread(i, v2, vcp2_y_minus);
      mult_z_plus2_thread(i, v2, vcp2_z_plus);
      mult_z_minus2_thread(i, v2, vcp2_z_minus);
      mult_t_plus2_chiral_thread(i, v2, vcp2_t_plus);
      mult_t_minus2_chiral_thread(i, v2, vcp2_t_minus);
    }

    // for (int i = is; i < is + ns; ++i) {
    // daypx_thread(i, v2, -m_kappa, v1);   // w = -m_kappa * w + f.
    //  daypx_thread(i, v2, -m_kappa_s, v1);   // w = -m_kappa * w + f.
    // }


    ThreadManager_OpenMP::sync_barrier_all();

    // clear(w);
    // mult_x_plus(w,f);
    // mult_x_minus(w,f);
    // mult_y_plus(w,f);
    // mult_y_minus(w,f);
    // mult_z_plus(w,f);
    // mult_z_minus(w,f);
    // mult_t_plus_chiral(w,f);
    // mult_t_minus_chiral(w,f);
    // daypx(w, -m_kappa, f);   // w = -m_kappa * w + f.
  }


//====================================================================
  void Fopr_WilsonGeneral::mult_gm5p(const int mu,
                                     Field& v, const Field& w)
  {
    clear(m_w2);

    mult_up(mu, m_w2, w);
    mult_gm5(v, m_w2);
  }


//====================================================================
  void Fopr_WilsonGeneral::proj_chiral(Field& w, const int ex1,
                                       const Field& v, const int ex2,
                                       const int ipm)
  {
    double fpm = 0.0;

    if (ipm == 1) {
      fpm = 1.0;
    } else if (ipm == -1) {
      fpm = -1.0;
    } else {
      vout.crucial(m_vl, "Error at %s: illegal chirality = %d.\n", class_name.c_str(), ipm);
      exit(EXIT_FAILURE);
    }

    m_w1.setpart_ex(0, v, ex2);
    mult_gm5(m_w2, m_w1);
    m_w1.addpart_ex(0, m_w2, 0, fpm);
    w.setpart_ex(ex1, m_w1, 0);

#pragma omp barrier
  }


//====================================================================
  void Fopr_WilsonGeneral::daxpy(Field& w,
                                 const double fac, const Field& f)
  {
    const double *v1 = f.ptr(0);
    double       *v2 = w.ptr(0);

    const int Nthread  = ThreadManager_OpenMP::get_num_threads();
    const int i_thread = ThreadManager_OpenMP::get_thread_id();

    const int is = m_Ntask * i_thread / Nthread;
    const int ns = m_Ntask * (i_thread + 1) / Nthread - is;

    for (int i = is; i < is + ns; ++i) {
      daxpy_thread(i, v2, fac, v1);
    }
#pragma omp barrier
  }


//====================================================================
  void Fopr_WilsonGeneral::daypx(Field& w,
                                 const double fac, const Field& f)
  {
    const double *v1 = f.ptr(0);
    double       *v2 = w.ptr(0);

    const int Nthread  = ThreadManager_OpenMP::get_num_threads();
    const int i_thread = ThreadManager_OpenMP::get_thread_id();

    const int is = m_Ntask * i_thread / Nthread;
    const int ns = m_Ntask * (i_thread + 1) / Nthread - is;

    for (int i = is; i < is + ns; ++i) {
      daypx_thread(i, v2, fac, v1);
    }
#pragma omp barrier
  }


//====================================================================
  void Fopr_WilsonGeneral::scal(Field& w,
                                const double fac)
  {
    double *v2 = w.ptr(0);

    const int Nthread  = ThreadManager_OpenMP::get_num_threads();
    const int i_thread = ThreadManager_OpenMP::get_thread_id();

    const int is = m_Ntask * i_thread / Nthread;
    const int ns = m_Ntask * (i_thread + 1) / Nthread - is;

    for (int i = is; i < is + ns; ++i) {
      scal_thread(i, v2, fac);
    }
#pragma omp barrier
  }


//====================================================================
  void Fopr_WilsonGeneral::clear(Field& w)
  {
    double *v2 = w.ptr(0);

    const int Nthread  = ThreadManager_OpenMP::get_num_threads();
    const int i_thread = ThreadManager_OpenMP::get_thread_id();

    const int is = m_Ntask * i_thread / Nthread;
    const int ns = m_Ntask * (i_thread + 1) / Nthread - is;

    for (int i = is; i < is + ns; ++i) {
      clear_thread(i, v2);
    }
#pragma omp barrier
  }


//====================================================================
  void Fopr_WilsonGeneral::gm5_dirac(Field& w, const Field& f)
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
  void Fopr_WilsonGeneral::gm5_chiral(Field& w, const Field& f)
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
  void Fopr_WilsonGeneral::mult_x_plus(Field& w, const Field& f)
  {
    const double *v1 = f.ptr(0);
    double       *v2 = w.ptr(0);

    const int Nthread  = ThreadManager_OpenMP::get_num_threads();
    const int i_thread = ThreadManager_OpenMP::get_thread_id();

    const int is = m_Ntask * i_thread / Nthread;
    const int ns = m_Ntask * (i_thread + 1) / Nthread - is;

#pragma omp barrier

    for (int i = is; i < is + ns; ++i) {
      mult_x_plus1_thread(i, vcp1_x_plus, v1);
    }

#pragma omp barrier
#pragma omp master
    {
      const int idir_x = 0;
      const int Nv     = 2 * m_Nvc * m_Nd * m_Ny * m_Nz * m_Nt;
      Communicator::exchange(Nv, vcp2_x_plus, vcp1_x_plus, idir_x, 1, 1);
    }
#pragma omp barrier

    for (int i = is; i < is + ns; ++i) {
      mult_x_plus_bulk_thread(i, v2, v1);
    }

    for (int i = is; i < is + ns; ++i) {
      mult_x_plus2_thread(i, v2, vcp2_x_plus);
    }

    ThreadManager_OpenMP::sync_barrier_all();
  }


//====================================================================
  void Fopr_WilsonGeneral::mult_x_minus(Field& w, const Field& f)
  {
    const double *v1 = f.ptr(0);
    double       *v2 = w.ptr(0);

    const int Nthread  = ThreadManager_OpenMP::get_num_threads();
    const int i_thread = ThreadManager_OpenMP::get_thread_id();

    const int is = m_Ntask * i_thread / Nthread;
    const int ns = m_Ntask * (i_thread + 1) / Nthread - is;

#pragma omp barrier

    for (int i = is; i < is + ns; ++i) {
      mult_x_minus1_thread(i, vcp1_x_minus, v1);
    }

#pragma omp barrier
#pragma omp master
    {
      const int idir_x = 0;
      const int Nv     = 2 * m_Nvc * m_Nd * m_Ny * m_Nz * m_Nt;
      Communicator::exchange(Nv, vcp2_x_minus, vcp1_x_minus, idir_x, -1, 2);
    }
#pragma omp barrier

    for (int i = is; i < is + ns; ++i) {
      mult_x_minus_bulk_thread(i, v2, v1);
    }

    for (int i = is; i < is + ns; ++i) {
      mult_x_minus2_thread(i, v2, vcp2_x_minus);
    }

    ThreadManager_OpenMP::sync_barrier_all();
  }


//====================================================================
  void Fopr_WilsonGeneral::mult_y_plus(Field& w, const Field& f)
  {
    const double *v1 = f.ptr(0);
    double       *v2 = w.ptr(0);

    const int Nthread  = ThreadManager_OpenMP::get_num_threads();
    const int i_thread = ThreadManager_OpenMP::get_thread_id();

    const int is = m_Ntask * i_thread / Nthread;
    const int ns = m_Ntask * (i_thread + 1) / Nthread - is;

#pragma omp barrier

    for (int i = is; i < is + ns; ++i) {
      mult_y_plus1_thread(i, vcp1_y_plus, v1);
    }

#pragma omp barrier
#pragma omp master
    {
      const int idir_y = 1;
      const int Nv     = 2 * m_Nvc * m_Nd * m_Nx * m_Nz * m_Nt;
      Communicator::exchange(Nv, vcp2_y_plus, vcp1_y_plus, idir_y, 1, 3);
    }
#pragma omp barrier

    for (int i = is; i < is + ns; ++i) {
      mult_y_plus_bulk_thread(i, v2, v1);
    }

    for (int i = is; i < is + ns; ++i) {
      mult_y_plus2_thread(i, v2, vcp2_y_plus);
    }

    ThreadManager_OpenMP::sync_barrier_all();
  }


//====================================================================
  void Fopr_WilsonGeneral::mult_y_minus(Field& w, const Field& f)
  {
    const double *v1 = f.ptr(0);
    double       *v2 = w.ptr(0);

    const int Nthread  = ThreadManager_OpenMP::get_num_threads();
    const int i_thread = ThreadManager_OpenMP::get_thread_id();

    const int is = m_Ntask * i_thread / Nthread;
    const int ns = m_Ntask * (i_thread + 1) / Nthread - is;

#pragma omp barrier

    for (int i = is; i < is + ns; ++i) {
      mult_y_minus1_thread(i, vcp1_y_minus, v1);
    }

#pragma omp barrier
#pragma omp master
    {
      const int idir_y = 1;
      const int Nv     = 2 * m_Nvc * m_Nd * m_Nx * m_Nz * m_Nt;
      Communicator::exchange(Nv, vcp2_y_minus, vcp1_y_minus, idir_y, -1, 4);
    }
#pragma omp barrier

    for (int i = is; i < is + ns; ++i) {
      mult_y_minus_bulk_thread(i, v2, v1);
    }

    for (int i = is; i < is + ns; ++i) {
      mult_y_minus2_thread(i, v2, vcp2_y_minus);
    }

    ThreadManager_OpenMP::sync_barrier_all();
  }


//====================================================================
  void Fopr_WilsonGeneral::mult_z_plus(Field& w, const Field& f)
  {
    const double *v1 = f.ptr(0);
    double       *v2 = w.ptr(0);

    const int Nthread  = ThreadManager_OpenMP::get_num_threads();
    const int i_thread = ThreadManager_OpenMP::get_thread_id();

    const int is = m_Ntask * i_thread / Nthread;
    const int ns = m_Ntask * (i_thread + 1) / Nthread - is;

#pragma omp barrier

    for (int i = is; i < is + ns; ++i) {
      mult_z_plus1_thread(i, vcp1_z_plus, v1);
    }

#pragma omp barrier
#pragma omp master
    {
      const int idir_z = 2;
      const int Nv     = 2 * m_Nvc * m_Nd * m_Nx * m_Ny * m_Nt;
      Communicator::exchange(Nv, vcp2_z_plus, vcp1_z_plus, idir_z, 1, 5);
    }
#pragma omp barrier

    for (int i = is; i < is + ns; ++i) {
      mult_z_plus_bulk_thread(i, v2, v1);
    }

    for (int i = is; i < is + ns; ++i) {
      mult_z_plus2_thread(i, v2, vcp2_z_plus);
    }

    ThreadManager_OpenMP::sync_barrier_all();
  }


//====================================================================
  void Fopr_WilsonGeneral::mult_z_minus(Field& w, const Field& f)
  {
    const double *v1 = f.ptr(0);
    double       *v2 = w.ptr(0);

    const int Nthread  = ThreadManager_OpenMP::get_num_threads();
    const int i_thread = ThreadManager_OpenMP::get_thread_id();

    const int is = m_Ntask * i_thread / Nthread;
    const int ns = m_Ntask * (i_thread + 1) / Nthread - is;

#pragma omp barrier

    for (int i = is; i < is + ns; ++i) {
      mult_z_minus1_thread(i, vcp1_z_minus, v1);
    }

#pragma omp barrier
#pragma omp master
    {
      const int idir_z = 2;
      const int Nv     = 2 * m_Nvc * m_Nd * m_Nx * m_Ny * m_Nt;
      Communicator::exchange(Nv, vcp2_z_minus, vcp1_z_minus, idir_z, -1, 6);
    }
#pragma omp barrier

    for (int i = is; i < is + ns; ++i) {
      mult_z_minus_bulk_thread(i, v2, v1);
    }

    for (int i = is; i < is + ns; ++i) {
      mult_z_minus2_thread(i, v2, vcp2_z_minus);
    }

    ThreadManager_OpenMP::sync_barrier_all();
  }


//====================================================================
  void Fopr_WilsonGeneral::mult_t_plus_dirac(Field& w, const Field& f)
  {
    const double *v1 = f.ptr(0);
    double       *v2 = w.ptr(0);

    const int Nthread  = ThreadManager_OpenMP::get_num_threads();
    const int i_thread = ThreadManager_OpenMP::get_thread_id();

    const int is = m_Ntask * i_thread / Nthread;
    const int ns = m_Ntask * (i_thread + 1) / Nthread - is;

#pragma omp barrier

    for (int i = is; i < is + ns; ++i) {
      mult_t_plus1_dirac_thread(i, vcp1_t_plus, v1);
    }

#pragma omp barrier
#pragma omp master
    {
      const int idir_t = 3;
      // const int Nv = 2 * m_Nvc * m_Nd * m_Nx * m_Ny * m_Nz;
      const int Nv = m_Nvc * m_Nd * m_Nx * m_Ny * m_Nz;
      Communicator::exchange(Nv, vcp2_t_plus, vcp1_t_plus, idir_t, 1, 7);
    }
#pragma omp barrier

    for (int i = is; i < is + ns; ++i) {
      mult_t_plus_bulk_dirac_thread(i, v2, v1);
    }

    for (int i = is; i < is + ns; ++i) {
      mult_t_plus2_dirac_thread(i, v2, vcp2_t_plus);
    }

    ThreadManager_OpenMP::sync_barrier_all();
  }


//====================================================================
  void Fopr_WilsonGeneral::mult_t_minus_dirac(Field& w, const Field& f)
  {
    const double *v1 = f.ptr(0);
    double       *v2 = w.ptr(0);

    const int Nthread  = ThreadManager_OpenMP::get_num_threads();
    const int i_thread = ThreadManager_OpenMP::get_thread_id();

    const int is = m_Ntask * i_thread / Nthread;
    const int ns = m_Ntask * (i_thread + 1) / Nthread - is;

#pragma omp barrier

    for (int i = is; i < is + ns; ++i) {
      mult_t_minus1_dirac_thread(i, vcp1_t_minus, v1);
    }

#pragma omp barrier
#pragma omp master
    {
      const int idir_t = 3;
      // const int Nv = 2 * m_Nvc * m_Nd * m_Nx * m_Ny * m_Nz;
      const int Nv = m_Nvc * m_Nd * m_Nx * m_Ny * m_Nz;
      Communicator::exchange(Nv, vcp2_t_minus, vcp1_t_minus, idir_t, -1, 8);
    }
#pragma omp barrier

    for (int i = is; i < is + ns; ++i) {
      mult_t_minus_bulk_dirac_thread(i, v2, v1);
    }

    for (int i = is; i < is + ns; ++i) {
      mult_t_minus2_dirac_thread(i, v2, vcp2_t_minus);
    }

    ThreadManager_OpenMP::sync_barrier_all();
  }


//====================================================================
  void Fopr_WilsonGeneral::mult_t_plus_chiral(Field& w, const Field& f)
  {
    const double *v1 = f.ptr(0);
    double       *v2 = w.ptr(0);

    const int Nthread  = ThreadManager_OpenMP::get_num_threads();
    const int i_thread = ThreadManager_OpenMP::get_thread_id();

    const int is = m_Ntask * i_thread / Nthread;
    const int ns = m_Ntask * (i_thread + 1) / Nthread - is;

#pragma omp barrier

    for (int i = is; i < is + ns; ++i) {
      mult_t_plus1_chiral_thread(i, vcp1_t_plus, v1);
    }

#pragma omp barrier
#pragma omp master
    {
      const int idir_t = 3;
      // const int Nv = 2 * m_Nvc * m_Nd * m_Nx * m_Ny * m_Nz;
      const int Nv = m_Nvc * m_Nd * m_Nx * m_Ny * m_Nz;
      Communicator::exchange(Nv, vcp2_t_plus, vcp1_t_plus, idir_t, 1, 7);
    }
#pragma omp barrier

    for (int i = is; i < is + ns; ++i) {
      mult_t_plus_bulk_chiral_thread(i, v2, v1);
    }

    for (int i = is; i < is + ns; ++i) {
      mult_t_plus2_chiral_thread(i, v2, vcp2_t_plus);
    }

    ThreadManager_OpenMP::sync_barrier_all();
  }


//====================================================================
  void Fopr_WilsonGeneral::mult_t_minus_chiral(Field& w, const Field& f)
  {
    const double *v1 = f.ptr(0);
    double       *v2 = w.ptr(0);

    const int Nthread  = ThreadManager_OpenMP::get_num_threads();
    const int i_thread = ThreadManager_OpenMP::get_thread_id();

    const int is = m_Ntask * i_thread / Nthread;
    const int ns = m_Ntask * (i_thread + 1) / Nthread - is;

#pragma omp barrier

    for (int i = is; i < is + ns; ++i) {
      mult_t_minus1_chiral_thread(i, vcp1_t_minus, v1);
    }

#pragma omp barrier
#pragma omp master
    {
      const int idir_t = 3;
      // const int Nv = 2 * m_Nvc * m_Nd * m_Nx * m_Ny * m_Nz;
      const int Nv = m_Nvc * m_Nd * m_Nx * m_Ny * m_Nz;
      Communicator::exchange(Nv, vcp2_t_minus, vcp1_t_minus, idir_t, -1, 8);
    }
#pragma omp barrier

    for (int i = is; i < is + ns; ++i) {
      mult_t_minus_bulk_chiral_thread(i, v2, v1);
    }

    for (int i = is; i < is + ns; ++i) {
      mult_t_minus2_chiral_thread(i, v2, vcp2_t_minus);
    }

    ThreadManager_OpenMP::sync_barrier_all();
  }


//====================================================================
  double Fopr_WilsonGeneral::flop_count()
  {
    // Counting of floating point operations in giga unit.
    // The following counting explicitly depends on the implementation.
    // It will be recalculated when the code is modified.
    // The present counting is based on rev.1107. [24 Aug 2014 H.Matsufuru]

    const int Nvol = CommonParameters::Nvol();
    const int NPE  = CommonParameters::NPE();

    int flop_site;

    if (m_repr == "Dirac") {
      flop_site = m_Nc * m_Nd * (4 + 6 * (4 * m_Nc + 2) + 2 * (4 * m_Nc + 1));
    } else if (m_repr == "Chiral") {
      flop_site = m_Nc * m_Nd * (4 + 8 * (4 * m_Nc + 2));
    } else {
      vout.crucial(m_vl, "Error at %s: input repr is undefined.\n",
                   class_name.c_str());
      exit(EXIT_FAILURE);
    }

    double gflop = flop_site * (Nvol * (NPE / 1.0e+9));

    if ((m_mode == "DdagD") || (m_mode == "DDdag")) gflop *= 2;

    return gflop;
  }


//====================================================================
}
//============================================================END=====
#endif
