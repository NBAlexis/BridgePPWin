#include "BridgeLib_Private.h"
#if USE_IMP_BGQ

/*!
        @file    $Id:: fopr_CloverTerm.cpp #$

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 00:49:46 #$

        @version $LastChangedRevision: 1561 $
*/

#include "Fopr/fopr_CloverTerm.h"

#include "ResourceManager/threadManager_OpenMP.h"


using std::string;

namespace Imp_BGQ {
#if defined USE_GROUP_SU3
#include "fopr_Wilson_impl_SU3.inc"
#elif defined USE_GROUP_SU2
#include "fopr_Wilson_impl_SU2.inc"
#elif defined USE_GROUP_SU_N
#include "fopr_Wilson_impl_SU_N.inc"
#endif

//====================================================================

  const std::string Fopr_CloverTerm::class_name = "Imp_BGQ::Fopr_CloverTerm";

//====================================================================
  void Fopr_CloverTerm::set_parameters(const Parameters& params)
  {
    const string str_vlevel = params.get_string("verbose_level");

    m_vl = vout.set_verbose_level(str_vlevel);

    //- fetch and check input parameters
    double           kappa, cSW;
    std::vector<int> bc;

    int err = 0;
    err += params.fetch_double("hopping_parameter", kappa);
    err += params.fetch_double("clover_coefficient", cSW);
    err += params.fetch_int_vector("boundary_condition", bc);

    if (err) {
      vout.crucial(m_vl, "Error at %s: input parameter not found.\n", class_name.c_str());
      exit(EXIT_FAILURE);
    }


    set_parameters(kappa, cSW, bc);
  }


//====================================================================
  void Fopr_CloverTerm::set_parameters(double kappa, double cSW,
                                       std::vector<int> bc)
  {
    //- print input parameters
    vout.general(m_vl, "%s:\n", class_name.c_str());
    vout.general(m_vl, "  kappa = %12.8f\n", kappa);
    vout.general(m_vl, "  cSW   = %12.8f\n", cSW);
    for (int mu = 0; mu < m_Ndim; ++mu) {
      vout.general(m_vl, "  boundary[%d] = %2d\n", mu, bc[mu]);
    }

    //- range check
    // NB. kappa,cSW == 0 is allowed.
    assert(bc.size() == m_Ndim);

    //- store values
    m_kappa = kappa;
    m_cSW   = cSW;

    // m_boundary.resize(m_Ndim);  // already resized in init.
    for (int mu = 0; mu < m_Ndim; ++mu) {
      m_boundary[mu] = bc[mu];
    }
  }


//====================================================================
  void Fopr_CloverTerm::set_config(Field *U)
  {
    m_U = (Field_G *)U;
    set_csw();
  }


//====================================================================
  void Fopr_CloverTerm::init(string repr)
  {
    m_Nvol = CommonParameters::Nvol();
    m_Ndim = CommonParameters::Ndim();
    m_Nc   = CommonParameters::Nc();
    m_Nd   = CommonParameters::Nd();
    m_NinF = 2 * m_Nc * m_Nd;

    m_U = 0;

    m_repr = repr;

    m_boundary.resize(m_Ndim);
    m_SG.resize(m_Ndim * m_Ndim);

    GammaMatrixSet *gmset = GammaMatrixSet::New(m_repr);

    m_GM5 = gmset->get_GM(gmset->GAMMA5);

    m_SG[sg_index(0, 1)] = gmset->get_GM(gmset->SIGMA12);
    m_SG[sg_index(1, 2)] = gmset->get_GM(gmset->SIGMA23);
    m_SG[sg_index(2, 0)] = gmset->get_GM(gmset->SIGMA31);
    m_SG[sg_index(3, 0)] = gmset->get_GM(gmset->SIGMA41);
    m_SG[sg_index(3, 1)] = gmset->get_GM(gmset->SIGMA42);
    m_SG[sg_index(3, 2)] = gmset->get_GM(gmset->SIGMA43);

    m_SG[sg_index(1, 0)] = m_SG[sg_index(0, 1)].mult(-1);
    m_SG[sg_index(2, 1)] = m_SG[sg_index(1, 2)].mult(-1);
    m_SG[sg_index(0, 2)] = m_SG[sg_index(2, 0)].mult(-1);
    m_SG[sg_index(0, 3)] = m_SG[sg_index(3, 0)].mult(-1);
    m_SG[sg_index(1, 3)] = m_SG[sg_index(3, 1)].mult(-1);
    m_SG[sg_index(2, 3)] = m_SG[sg_index(3, 2)].mult(-1);

    m_SG[sg_index(0, 0)] = gmset->get_GM(gmset->UNITY);
    m_SG[sg_index(1, 1)] = gmset->get_GM(gmset->UNITY);
    m_SG[sg_index(2, 2)] = gmset->get_GM(gmset->UNITY);
    m_SG[sg_index(3, 3)] = gmset->get_GM(gmset->UNITY);
    // these 4 gamma matrices are actually not used.

    //  m_fopr_w = new Fopr_Wilson(repr);
    if (m_repr == "Dirac") {
      m_csw = &Fopr_CloverTerm::mult_csw_dirac;
      m_gm5 = &Fopr_CloverTerm::gm5_dirac;
    } else if (m_repr == "Chiral") {
      m_csw = &Fopr_CloverTerm::mult_csw_chiral;
      m_gm5 = &Fopr_CloverTerm::gm5_chiral;
    }

    delete gmset;

    //  m_w2.reset(m_NinF,m_Nvol,1);
  }


//====================================================================
  void Fopr_CloverTerm::tidyup()
  {
    //  nothing to do.
  }


//====================================================================
  void Fopr_CloverTerm::mult_gm5(Field& v, const Field& f)
  {
    (this->*m_gm5)(v, f);
  }


//====================================================================
  void Fopr_CloverTerm::gm5_dirac(Field& w, const Field& f)
  {
    int Nvc = 2 * CommonParameters::Nc();
    int Nd  = CommonParameters::Nd();

    const double *v1 = f.ptr(0);
    double       *v2 = w.ptr(0);

    int id1 = 0;
    int id2 = Nvc;
    int id3 = Nvc * 2;
    int id4 = Nvc * 3;

    // threadding applied.
    int nth = ThreadManager_OpenMP::get_num_threads();
    int ith = ThreadManager_OpenMP::get_thread_id();

    int is = m_Nvol * ith / nth;
    int ns = m_Nvol * (ith + 1) / nth - is;

    for (int site = is; site < is + ns; ++site) {
      // for (int site = 0; site < m_Nvol; ++site) {
      for (int icc = 0; icc < Nvc; icc++) {
        int in = Nvc * Nd * site;
        v2[icc + id1 + in] = v1[icc + id3 + in];
        v2[icc + id2 + in] = v1[icc + id4 + in];
        v2[icc + id3 + in] = v1[icc + id1 + in];
        v2[icc + id4 + in] = v1[icc + id2 + in];
      }
    }
  }


//====================================================================
  void Fopr_CloverTerm::gm5_chiral(Field& w, const Field& f)
  {
    int Nvc = 2 * CommonParameters::Nc();
    int Nd  = CommonParameters::Nd();

    const double *v1 = f.ptr(0);
    double       *v2 = w.ptr(0);

    int id1 = 0;
    int id2 = Nvc;
    int id3 = Nvc * 2;
    int id4 = Nvc * 3;

    // threadding applied.
    int nth = ThreadManager_OpenMP::get_num_threads();
    int ith = ThreadManager_OpenMP::get_thread_id();

    int is = m_Nvol * ith / nth;
    int ns = m_Nvol * (ith + 1) / nth - is;

    for (int site = is; site < is + ns; ++site) {
      // for (int site = 0; site < m_Nvol; ++site) {
      for (int icc = 0; icc < Nvc; icc++) {
        int in = Nvc * Nd * site;
        v2[icc + id1 + in] = v1[icc + id1 + in];
        v2[icc + id2 + in] = v1[icc + id2 + in];
        v2[icc + id3 + in] = -v1[icc + id3 + in];
        v2[icc + id4 + in] = -v1[icc + id4 + in];
      }
    }
  }


//====================================================================
  void Fopr_CloverTerm::mult_isigma(Field_F& v, const Field_F& w,
                                    const int mu, const int nu)
  {
    assert(mu != nu);
    mult_iGM(v, m_SG[sg_index(mu, nu)], w);
  }


//====================================================================
  void Fopr_CloverTerm::mult_sigmaF(Field& v, const Field& f)
  {
    mult_csw(v, f);
  }


//====================================================================
  void Fopr_CloverTerm::mult_csw(Field& v, const Field& w)
  {
    (this->*m_csw)(v, w);
  }


//====================================================================
  void Fopr_CloverTerm::mult_csw_chiral(Field& v, const Field& w)
  {
    assert(w.nex() == 1);

    int Nc   = CommonParameters::Nc();
    int Nvc  = 2 * Nc;
    int Ndf  = 2 * Nc * Nc;
    int Nd   = CommonParameters::Nd();
    int Nvol = w.nvol();

    int id1 = 0;
    int id2 = Nvc;
    int id3 = Nvc * 2;
    int id4 = Nvc * 3;

    double kappa_cSW = m_kappa * m_cSW;

    const double *w2 = w.ptr(0);
    double       *v2 = v.ptr(0);

    double *Bx = m_Bx.ptr(0);
    double *By = m_By.ptr(0);
    double *Bz = m_Bz.ptr(0);
    double *Ex = m_Ex.ptr(0);
    double *Ey = m_Ey.ptr(0);
    double *Ez = m_Ez.ptr(0);

    // threadding applied.
    int nth = ThreadManager_OpenMP::get_num_threads();
    int ith = ThreadManager_OpenMP::get_thread_id();

    int is = m_Nvol * ith / nth;
    int ns = m_Nvol * (ith + 1) / nth - is;

    for (int site = is; site < is + ns; ++site) {
      int iv = Nvc * Nd * site;
      int ig = Ndf * site;

      for (int ic = 0; ic < Nc; ++ic) {
        int icr = 2 * ic;
        int ici = icr + 1;
        int icg = ic * Nvc + ig;

        v2[icr + id1 + iv] = 0.0;
        v2[ici + id1 + iv] = 0.0;
        v2[icr + id2 + iv] = 0.0;
        v2[ici + id2 + iv] = 0.0;

        v2[icr + id3 + iv] = 0.0;
        v2[ici + id3 + iv] = 0.0;
        v2[icr + id4 + iv] = 0.0;
        v2[ici + id4 + iv] = 0.0;

        // isigma_23 * Bx
        v2[icr + id1 + iv] -= mult_uv_i(&Bx[icg], &w2[id2 + iv], Nc);
        v2[ici + id1 + iv] += mult_uv_r(&Bx[icg], &w2[id2 + iv], Nc);
        v2[icr + id2 + iv] -= mult_uv_i(&Bx[icg], &w2[id1 + iv], Nc);
        v2[ici + id2 + iv] += mult_uv_r(&Bx[icg], &w2[id1 + iv], Nc);

        v2[icr + id3 + iv] -= mult_uv_i(&Bx[icg], &w2[id4 + iv], Nc);
        v2[ici + id3 + iv] += mult_uv_r(&Bx[icg], &w2[id4 + iv], Nc);
        v2[icr + id4 + iv] -= mult_uv_i(&Bx[icg], &w2[id3 + iv], Nc);
        v2[ici + id4 + iv] += mult_uv_r(&Bx[icg], &w2[id3 + iv], Nc);

        // isigma_31 * By
        v2[icr + id1 + iv] += mult_uv_r(&By[icg], &w2[id2 + iv], Nc);
        v2[ici + id1 + iv] += mult_uv_i(&By[icg], &w2[id2 + iv], Nc);
        v2[icr + id2 + iv] -= mult_uv_r(&By[icg], &w2[id1 + iv], Nc);
        v2[ici + id2 + iv] -= mult_uv_i(&By[icg], &w2[id1 + iv], Nc);

        v2[icr + id3 + iv] += mult_uv_r(&By[icg], &w2[id4 + iv], Nc);
        v2[ici + id3 + iv] += mult_uv_i(&By[icg], &w2[id4 + iv], Nc);
        v2[icr + id4 + iv] -= mult_uv_r(&By[icg], &w2[id3 + iv], Nc);
        v2[ici + id4 + iv] -= mult_uv_i(&By[icg], &w2[id3 + iv], Nc);

        // isigma_12 * Bz
        v2[icr + id1 + iv] -= mult_uv_i(&Bz[icg], &w2[id1 + iv], Nc);
        v2[ici + id1 + iv] += mult_uv_r(&Bz[icg], &w2[id1 + iv], Nc);
        v2[icr + id2 + iv] += mult_uv_i(&Bz[icg], &w2[id2 + iv], Nc);
        v2[ici + id2 + iv] -= mult_uv_r(&Bz[icg], &w2[id2 + iv], Nc);

        v2[icr + id3 + iv] -= mult_uv_i(&Bz[icg], &w2[id3 + iv], Nc);
        v2[ici + id3 + iv] += mult_uv_r(&Bz[icg], &w2[id3 + iv], Nc);
        v2[icr + id4 + iv] += mult_uv_i(&Bz[icg], &w2[id4 + iv], Nc);
        v2[ici + id4 + iv] -= mult_uv_r(&Bz[icg], &w2[id4 + iv], Nc);

        // isigma_41 * Ex
        v2[icr + id1 + iv] += mult_uv_i(&Ex[icg], &w2[id2 + iv], Nc);
        v2[ici + id1 + iv] -= mult_uv_r(&Ex[icg], &w2[id2 + iv], Nc);
        v2[icr + id2 + iv] += mult_uv_i(&Ex[icg], &w2[id1 + iv], Nc);
        v2[ici + id2 + iv] -= mult_uv_r(&Ex[icg], &w2[id1 + iv], Nc);

        v2[icr + id3 + iv] -= mult_uv_i(&Ex[icg], &w2[id4 + iv], Nc);
        v2[ici + id3 + iv] += mult_uv_r(&Ex[icg], &w2[id4 + iv], Nc);
        v2[icr + id4 + iv] -= mult_uv_i(&Ex[icg], &w2[id3 + iv], Nc);
        v2[ici + id4 + iv] += mult_uv_r(&Ex[icg], &w2[id3 + iv], Nc);

        // isigma_42 * Ey
        v2[icr + id1 + iv] -= mult_uv_r(&Ey[icg], &w2[id2 + iv], Nc);
        v2[ici + id1 + iv] -= mult_uv_i(&Ey[icg], &w2[id2 + iv], Nc);
        v2[icr + id2 + iv] += mult_uv_r(&Ey[icg], &w2[id1 + iv], Nc);
        v2[ici + id2 + iv] += mult_uv_i(&Ey[icg], &w2[id1 + iv], Nc);

        v2[icr + id3 + iv] += mult_uv_r(&Ey[icg], &w2[id4 + iv], Nc);
        v2[ici + id3 + iv] += mult_uv_i(&Ey[icg], &w2[id4 + iv], Nc);
        v2[icr + id4 + iv] -= mult_uv_r(&Ey[icg], &w2[id3 + iv], Nc);
        v2[ici + id4 + iv] -= mult_uv_i(&Ey[icg], &w2[id3 + iv], Nc);

        // isigma_43 * Ez
        v2[icr + id1 + iv] += mult_uv_i(&Ez[icg], &w2[id1 + iv], Nc);
        v2[ici + id1 + iv] -= mult_uv_r(&Ez[icg], &w2[id1 + iv], Nc);
        v2[icr + id2 + iv] -= mult_uv_i(&Ez[icg], &w2[id2 + iv], Nc);
        v2[ici + id2 + iv] += mult_uv_r(&Ez[icg], &w2[id2 + iv], Nc);

        v2[icr + id3 + iv] -= mult_uv_i(&Ez[icg], &w2[id3 + iv], Nc);
        v2[ici + id3 + iv] += mult_uv_r(&Ez[icg], &w2[id3 + iv], Nc);
        v2[icr + id4 + iv] += mult_uv_i(&Ez[icg], &w2[id4 + iv], Nc);
        v2[ici + id4 + iv] -= mult_uv_r(&Ez[icg], &w2[id4 + iv], Nc);

        // v *= m_kappa * m_cSW;
        v2[icr + id1 + iv] *= kappa_cSW;
        v2[ici + id1 + iv] *= kappa_cSW;
        v2[icr + id2 + iv] *= kappa_cSW;
        v2[ici + id2 + iv] *= kappa_cSW;

        v2[icr + id3 + iv] *= kappa_cSW;
        v2[ici + id3 + iv] *= kappa_cSW;
        v2[icr + id4 + iv] *= kappa_cSW;
        v2[ici + id4 + iv] *= kappa_cSW;
      }
    }
#pragma omp barrier
  }


//====================================================================
  void Fopr_CloverTerm::mult_csw_dirac(Field& v, const Field& w)
  {
    assert(w.nex() == 1);

    int Nc   = CommonParameters::Nc();
    int Nvc  = 2 * Nc;
    int Ndf  = 2 * Nc * Nc;
    int Nd   = CommonParameters::Nd();
    int Nvol = w.nvol();

    int id1 = 0;
    int id2 = Nvc;
    int id3 = Nvc * 2;
    int id4 = Nvc * 3;

    double kappa_cSW = m_kappa * m_cSW;

    const double *w2 = w.ptr(0);
    double       *v2 = v.ptr(0);

    double *Bx = m_Bx.ptr(0);
    double *By = m_By.ptr(0);
    double *Bz = m_Bz.ptr(0);
    double *Ex = m_Ex.ptr(0);
    double *Ey = m_Ey.ptr(0);
    double *Ez = m_Ez.ptr(0);

    // threadding applied.
    int nth = ThreadManager_OpenMP::get_num_threads();
    int ith = ThreadManager_OpenMP::get_thread_id();

    int is = m_Nvol * ith / nth;
    int ns = m_Nvol * (ith + 1) / nth - is;

    for (int site = is; site < is + ns; ++site) {
      int iv = Nvc * Nd * site;
      int ig = Ndf * site;

      for (int ic = 0; ic < Nc; ++ic) {
        int icr = 2 * ic;
        int ici = icr + 1;
        int icg = ic * Nvc + ig;

        v2[icr + id1 + iv] = 0.0;
        v2[ici + id1 + iv] = 0.0;
        v2[icr + id2 + iv] = 0.0;
        v2[ici + id2 + iv] = 0.0;

        v2[icr + id3 + iv] = 0.0;
        v2[ici + id3 + iv] = 0.0;
        v2[icr + id4 + iv] = 0.0;
        v2[ici + id4 + iv] = 0.0;

        // isigma_23 * Bx
        v2[icr + id1 + iv] -= mult_uv_i(&Bx[icg], &w2[id2 + iv], Nc);
        v2[ici + id1 + iv] += mult_uv_r(&Bx[icg], &w2[id2 + iv], Nc);
        v2[icr + id2 + iv] -= mult_uv_i(&Bx[icg], &w2[id1 + iv], Nc);
        v2[ici + id2 + iv] += mult_uv_r(&Bx[icg], &w2[id1 + iv], Nc);

        v2[icr + id3 + iv] -= mult_uv_i(&Bx[icg], &w2[id4 + iv], Nc);
        v2[ici + id3 + iv] += mult_uv_r(&Bx[icg], &w2[id4 + iv], Nc);
        v2[icr + id4 + iv] -= mult_uv_i(&Bx[icg], &w2[id3 + iv], Nc);
        v2[ici + id4 + iv] += mult_uv_r(&Bx[icg], &w2[id3 + iv], Nc);

        // isigma_31 * By
        v2[icr + id1 + iv] += mult_uv_r(&By[icg], &w2[id2 + iv], Nc);
        v2[ici + id1 + iv] += mult_uv_i(&By[icg], &w2[id2 + iv], Nc);
        v2[icr + id2 + iv] -= mult_uv_r(&By[icg], &w2[id1 + iv], Nc);
        v2[ici + id2 + iv] -= mult_uv_i(&By[icg], &w2[id1 + iv], Nc);

        v2[icr + id3 + iv] += mult_uv_r(&By[icg], &w2[id4 + iv], Nc);
        v2[ici + id3 + iv] += mult_uv_i(&By[icg], &w2[id4 + iv], Nc);
        v2[icr + id4 + iv] -= mult_uv_r(&By[icg], &w2[id3 + iv], Nc);
        v2[ici + id4 + iv] -= mult_uv_i(&By[icg], &w2[id3 + iv], Nc);

        // isigma_12 * Bz
        v2[icr + id1 + iv] -= mult_uv_i(&Bz[icg], &w2[id1 + iv], Nc);
        v2[ici + id1 + iv] += mult_uv_r(&Bz[icg], &w2[id1 + iv], Nc);
        v2[icr + id2 + iv] += mult_uv_i(&Bz[icg], &w2[id2 + iv], Nc);
        v2[ici + id2 + iv] -= mult_uv_r(&Bz[icg], &w2[id2 + iv], Nc);

        v2[icr + id3 + iv] -= mult_uv_i(&Bz[icg], &w2[id3 + iv], Nc);
        v2[ici + id3 + iv] += mult_uv_r(&Bz[icg], &w2[id3 + iv], Nc);
        v2[icr + id4 + iv] += mult_uv_i(&Bz[icg], &w2[id4 + iv], Nc);
        v2[ici + id4 + iv] -= mult_uv_r(&Bz[icg], &w2[id4 + iv], Nc);

        // isigma_41 * Ex
        v2[icr + id1 + iv] += mult_uv_i(&Ex[icg], &w2[id4 + iv], Nc);
        v2[ici + id1 + iv] -= mult_uv_r(&Ex[icg], &w2[id4 + iv], Nc);
        v2[icr + id2 + iv] += mult_uv_i(&Ex[icg], &w2[id3 + iv], Nc);
        v2[ici + id2 + iv] -= mult_uv_r(&Ex[icg], &w2[id3 + iv], Nc);

        v2[icr + id3 + iv] += mult_uv_i(&Ex[icg], &w2[id2 + iv], Nc);
        v2[ici + id3 + iv] -= mult_uv_r(&Ex[icg], &w2[id2 + iv], Nc);
        v2[icr + id4 + iv] += mult_uv_i(&Ex[icg], &w2[id1 + iv], Nc);
        v2[ici + id4 + iv] -= mult_uv_r(&Ex[icg], &w2[id1 + iv], Nc);

        // isigma_42 * Ey
        v2[icr + id1 + iv] -= mult_uv_r(&Ey[icg], &w2[id4 + iv], Nc);
        v2[ici + id1 + iv] -= mult_uv_i(&Ey[icg], &w2[id4 + iv], Nc);
        v2[icr + id2 + iv] += mult_uv_r(&Ey[icg], &w2[id3 + iv], Nc);
        v2[ici + id2 + iv] += mult_uv_i(&Ey[icg], &w2[id3 + iv], Nc);

        v2[icr + id3 + iv] -= mult_uv_r(&Ey[icg], &w2[id2 + iv], Nc);
        v2[ici + id3 + iv] -= mult_uv_i(&Ey[icg], &w2[id2 + iv], Nc);
        v2[icr + id4 + iv] += mult_uv_r(&Ey[icg], &w2[id1 + iv], Nc);
        v2[ici + id4 + iv] += mult_uv_i(&Ey[icg], &w2[id1 + iv], Nc);

        // isigma_43 * Ez
        v2[icr + id1 + iv] += mult_uv_i(&Ez[icg], &w2[id3 + iv], Nc);
        v2[ici + id1 + iv] -= mult_uv_r(&Ez[icg], &w2[id3 + iv], Nc);
        v2[icr + id2 + iv] -= mult_uv_i(&Ez[icg], &w2[id4 + iv], Nc);
        v2[ici + id2 + iv] += mult_uv_r(&Ez[icg], &w2[id4 + iv], Nc);

        v2[icr + id3 + iv] += mult_uv_i(&Ez[icg], &w2[id1 + iv], Nc);
        v2[ici + id3 + iv] -= mult_uv_r(&Ez[icg], &w2[id1 + iv], Nc);
        v2[icr + id4 + iv] -= mult_uv_i(&Ez[icg], &w2[id2 + iv], Nc);
        v2[ici + id4 + iv] += mult_uv_r(&Ez[icg], &w2[id2 + iv], Nc);

        // v *= m_kappa * m_cSW;
        v2[icr + id1 + iv] *= kappa_cSW;
        v2[ici + id1 + iv] *= kappa_cSW;
        v2[icr + id2 + iv] *= kappa_cSW;
        v2[ici + id2 + iv] *= kappa_cSW;

        v2[icr + id3 + iv] *= kappa_cSW;
        v2[ici + id3 + iv] *= kappa_cSW;
        v2[icr + id4 + iv] *= kappa_cSW;
        v2[ici + id4 + iv] *= kappa_cSW;
      }
    }
#pragma omp barrier
  }


//====================================================================
  void Fopr_CloverTerm::set_csw()
  {
    set_fieldstrength(m_Bx, 1, 2);
    set_fieldstrength(m_By, 2, 0);
    set_fieldstrength(m_Bz, 0, 1);
    set_fieldstrength(m_Ex, 3, 0);
    set_fieldstrength(m_Ey, 3, 1);
    set_fieldstrength(m_Ez, 3, 2);
  }


//====================================================================
  void Fopr_CloverTerm::set_fieldstrength(Field_G& Fst,
                                          const int mu, const int nu)
  {
    int nth = ThreadManager_OpenMP::get_num_threads();

    assert(nth == 1);
    // this function must be called in single thread region.

    m_staple.upper(m_Cup, *m_U, mu, nu); // these staple constructions
    m_staple.lower(m_Cdn, *m_U, mu, nu); // are multi-threaded.

    //#pragma omp parallel
    {
      mult_Field_Gnd(Fst, 0, *m_U, mu, m_Cup, 0);
      multadd_Field_Gnd(Fst, 0, *m_U, mu, m_Cdn, 0, -1.0);

      mult_Field_Gdn(m_v1, 0, m_Cup, 0, *m_U, mu);
      multadd_Field_Gdn(m_v1, 0, m_Cdn, 0, *m_U, mu, -1.0);
    }

    m_shift.forward(m_v2, m_v1, mu);

    //#pragma omp parallel
    {
      axpy(Fst, 1.0, m_v2);
      ah_Field_G(Fst, 0);
      scal(Fst, 0.25);
    }
  }


//====================================================================
  double Fopr_CloverTerm::flop_count()
  {
    // The following counting explicitly depends on the implementation
    // and to be recalculated when the code is modified.
    // Present counting is based on rev.1107. [24 Aug 2014 H.Matsufuru]

    int Lvol = CommonParameters::Lvol();

    double flop_site
      = static_cast<double>(m_Nc * m_Nd * (2 + 48 * m_Nc));
    double flop = flop_site * static_cast<double>(Lvol);

    return flop;
  }


//====================================================================
}
//============================================================END=====
#endif
