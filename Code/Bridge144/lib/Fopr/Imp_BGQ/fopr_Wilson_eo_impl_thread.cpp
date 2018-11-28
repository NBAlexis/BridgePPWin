#include "BridgeLib_Private.h"
#if USE_IMP_BGQ

/*!
        @file    $Id:: fopr_Wilson_eo_impl_thread.cpp #$

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 00:49:46 #$

        @version $LastChangedRevision: 1561 $
*/

#include "fopr_Wilson_eo_impl.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

#include "ResourceManager/threadManager_OpenMP.h"

namespace Imp_BGQ {
#if defined USE_GROUP_SU3
#include "fopr_Wilson_impl_SU3.inc"
#elif defined USE_GROUP_SU2
#include "fopr_Wilson_impl_SU2.inc"
#elif defined USE_GROUP_SU_N
#include "fopr_Wilson_impl_SU_N.inc"
#endif

//====================================================================
  void Fopr_Wilson_eo::setup_thread()
  {
    m_Nthread = ThreadManager_OpenMP::get_num_threads_available();

    // The following setup corresponds to uniform division of volume.
    if (m_Nthread <= m_Nt) {
      m_Ntask_t = m_Nthread;
    } else if (m_Nthread <= m_Nz * m_Nt) {
      m_Ntask_t = m_Nt;
    } else {
      vout.crucial(m_vl, "Error at %s: Too large Nthread: %d\n",
                   class_name.c_str(), m_Nthread);
      exit(EXIT_FAILURE);
    }
    m_Ntask_z = m_Nthread / m_Ntask_t;
    if (m_Ntask_z * m_Ntask_t != m_Nthread) {
      vout.crucial(m_vl, "Error at %s: Nz(%d) and Nt(%d) \neq Nthread(%d)\n",
                   class_name.c_str(), m_Nz, m_Nt, m_Nthread);
      exit(EXIT_FAILURE);
    }
    m_Ntask = m_Ntask_t * m_Ntask_z;
    m_Mz    = m_Nz / m_Ntask_z;
    m_Mt    = m_Nt / m_Ntask_t;

    vout.general(m_vl, "  Nthread = %d\n", m_Nthread);
    vout.general(m_vl, "  Ntask   = %d\n", m_Ntask);
    vout.general(m_vl, "  Ntask_z = %d  Ntask_t = %d\n", m_Ntask_z, m_Ntask_t);
    vout.general(m_vl, "  Mz      = %d  Mt      = %d\n", m_Mz, m_Mt);

    // setup of arguments
    int Nxy2 = m_Nx2 * m_Ny;
    m_arg.resize(m_Ntask);
    for (int ith_t = 0; ith_t < m_Ntask_t; ++ith_t) {
      for (int ith_z = 0; ith_z < m_Ntask_z; ++ith_z) {
        int itask = ith_z + m_Ntask_z * ith_t;

        m_arg[itask].isite = (ith_z * m_Mz + ith_t * (m_Nz * m_Mt)) * Nxy2;

        m_arg[itask].kt0 = 0;
        m_arg[itask].kt1 = 0;
        m_arg[itask].kz0 = 0;
        m_arg[itask].kz1 = 0;
        if (ith_t == 0) m_arg[itask].kt0 = 1;
        if (ith_z == 0) m_arg[itask].kz0 = 1;
        if (ith_t == m_Ntask_t - 1) m_arg[itask].kt1 = 1;
        if (ith_z == m_Ntask_z - 1) m_arg[itask].kz1 = 1;

        m_arg[itask].isite_cpx = itask * m_Mz * m_Mt * (m_Ny / 2);
        m_arg[itask].isite_cpy = itask * m_Mz * m_Mt * m_Nx2;
        m_arg[itask].isite_cpz = ith_t * m_Mt * Nxy2;
        m_arg[itask].isite_cpt = ith_z * m_Mz * Nxy2;
      }
    }

    // setup for async data transfer
    int Nc    = CommonParameters::Nc();
    int Nd    = CommonParameters::Nd();
    int Nvcd2 = 2 * Nc * Nd / 2;

    std::vector<int> destid(m_Ntask);
    std::vector<int> offset(m_Ntask);
    std::vector<int> datasize(m_Ntask);
    std::vector<int> offset_up(m_Ntask);
    std::vector<int> offset_lw(m_Ntask);
    std::vector<int> datasize_up(m_Ntask);
    std::vector<int> datasize_lw(m_Ntask);

    int imu = 0;
    for (int ith_t = 0; ith_t < m_Ntask_t; ++ith_t) {
      for (int ith_z = 0; ith_z < m_Ntask_z; ++ith_z) {
        int itask    = ith_z + ith_t * m_Ntask_z;
        int isite_cp = itask * m_Mz * m_Mt * (m_Ny / 2);
        destid[itask]   = itask;
        offset[itask]   = sizeof(double) * Nvcd2 * isite_cp;
        datasize[itask] = sizeof(double) * Nvcd2 * m_Mz * m_Mt * (m_Ny / 2);
      }
    }
    m_bw_send[imu]->set_thread(m_Ntask, destid, offset, datasize);
    m_fw_send[imu]->set_thread(m_Ntask, destid, offset, datasize);
    m_bw_recv[imu]->set_thread(m_Ntask, destid, offset, datasize);
    m_fw_recv[imu]->set_thread(m_Ntask, destid, offset, datasize);

    imu = 1;
    for (int ith_t = 0; ith_t < m_Ntask_t; ++ith_t) {
      for (int ith_z = 0; ith_z < m_Ntask_z; ++ith_z) {
        int itask    = ith_z + ith_t * m_Ntask_z;
        int isite_cp = itask * m_Mz * m_Mt * m_Nx2;
        destid[itask]   = itask;
        offset[itask]   = sizeof(double) * Nvcd2 * isite_cp;
        datasize[itask] = sizeof(double) * Nvcd2 * m_Mz * m_Mt * m_Nx2;
      }
    }
    m_bw_send[imu]->set_thread(m_Ntask, destid, offset, datasize);
    m_fw_send[imu]->set_thread(m_Ntask, destid, offset, datasize);
    m_bw_recv[imu]->set_thread(m_Ntask, destid, offset, datasize);
    m_fw_recv[imu]->set_thread(m_Ntask, destid, offset, datasize);

    imu = 2;
    for (int ith_t = 0; ith_t < m_Ntask_t; ++ith_t) {
      for (int ith_z = 0; ith_z < m_Ntask_z; ++ith_z) {
        int itask = ith_z + m_Ntask_z * ith_t;
        destid[itask]      = -1;
        offset_up[itask]   = 0;
        offset_lw[itask]   = 0;
        datasize_up[itask] = 0;
        datasize_lw[itask] = 0;
        if (ith_z == 0) {
          destid[itask]      = (m_Ntask_z - 1) + ith_t * m_Ntask_z;
          offset_lw[itask]   = sizeof(double) * Nvcd2 * ith_t * m_Mt * m_Nx2 * m_Ny;
          datasize_lw[itask] = sizeof(double) * Nvcd2 * m_Mt * m_Nx2 * m_Ny;
        }
        if (ith_z == m_Ntask_z - 1) {
          destid[itask]      = ith_t * m_Ntask_z;
          offset_up[itask]   = sizeof(double) * Nvcd2 * ith_t * m_Mt * m_Nx2 * m_Ny;
          datasize_up[itask] = sizeof(double) * Nvcd2 * m_Mt * m_Nx2 * m_Ny;
        }
      }
    }
    m_bw_send[imu]->set_thread(m_Ntask, destid, offset_lw, datasize_lw);
    m_bw_recv[imu]->set_thread(m_Ntask, destid, offset_up, datasize_up);
    m_fw_send[imu]->set_thread(m_Ntask, destid, offset_up, datasize_up);
    m_fw_recv[imu]->set_thread(m_Ntask, destid, offset_lw, datasize_lw);

    imu = 3;
    for (int ith_t = 0; ith_t < m_Ntask_t; ++ith_t) {
      for (int ith_z = 0; ith_z < m_Ntask_z; ++ith_z) {
        int itask = ith_z + m_Ntask_z * ith_t;
        destid[itask]      = -1;
        offset_up[itask]   = 0;
        offset_lw[itask]   = 0;
        datasize_up[itask] = 0;
        datasize_lw[itask] = 0;
        if (ith_t == 0) {
          destid[itask]      = ith_z + (m_Ntask_t - 1) * m_Ntask_z;
          offset_lw[itask]   = sizeof(double) * Nvcd2 * ith_z * m_Mz * m_Nx2 * m_Ny;
          datasize_lw[itask] = sizeof(double) * Nvcd2 * m_Mz * m_Nx2 * m_Ny;
        }
        if (ith_t == m_Ntask_t - 1) {
          destid[itask]      = ith_z;
          offset_up[itask]   = sizeof(double) * Nvcd2 * ith_z * m_Mz * m_Nx2 * m_Ny;
          datasize_up[itask] = sizeof(double) * Nvcd2 * m_Mz * m_Nx2 * m_Ny;
        }
      }
    }
    m_bw_send[imu]->set_thread(m_Ntask, destid, offset_lw, datasize_lw);
    m_bw_recv[imu]->set_thread(m_Ntask, destid, offset_up, datasize_up);
    m_fw_send[imu]->set_thread(m_Ntask, destid, offset_up, datasize_up);
    m_fw_recv[imu]->set_thread(m_Ntask, destid, offset_lw, datasize_lw);
  }


//====================================================================
  void Fopr_Wilson_eo::scal_thread(int itask,
                                   double *w, double fac)
  {
    int Nvcd = m_Nvc * m_Nd;
    int Nvxy = Nvcd * m_Nx2 * m_Ny;

    int    isite = m_arg[itask].isite;
    double *wp   = &w[Nvcd * isite];

    for (int it = 0; it < m_Mt; ++it) {
      for (int iz = 0; iz < m_Mz; ++iz) {
        for (int ivxy = 0; ivxy < Nvxy; ++ivxy) {
          int iv = ivxy + Nvxy * (iz + m_Nz * it);
          wp[iv] = fac * wp[iv];
        }
      }
    }
  }


//====================================================================
  void Fopr_Wilson_eo::clear_thread(int itask,
                                    double *v)
  {
    int Nvcd = m_Nvc * m_Nd;
    int Nvxy = Nvcd * m_Nx2 * m_Ny;

    int    isite = m_arg[itask].isite;
    double *wp   = &v[Nvcd * isite];

    for (int it = 0; it < m_Mt; ++it) {
      for (int iz = 0; iz < m_Mz; ++iz) {
        for (int ivxy = 0; ivxy < Nvxy; ++ivxy) {
          int iv = ivxy + Nvxy * (iz + m_Nz * it);
          wp[iv] = 0.0;
        }
      }
    }
  }


//====================================================================
  void Fopr_Wilson_eo::mult_xp1_thread(
    int itask, double *vcp1, const double *v1, int ieo)
  {
    int Nvc2  = 2 * m_Nvc;
    int Nvcd  = m_Nvc * m_Nd;
    int Nvcd2 = Nvcd / 2;

    int id1 = 0;
    int id2 = m_Nvc;
    int id3 = m_Nvc * 2;
    int id4 = m_Nvc * 3;

    int idir = 0;

    int isite    = m_arg[itask].isite;
    int isite_cp = m_arg[itask].isite_cpx;
    int iyzt0    = isite / m_Nx2;

    //  double* w2 = &vcp1[Nvcd2*isite_cp];
    double *w2
      = (double *)m_bw_send[idir]->ptr(sizeof(double) * Nvcd2 * isite_cp);
    const double *w1 = &v1[Nvcd * isite];

    double bc2 = m_boundary2[idir];

    int ix  = 0;
    int ibf = 0;

    for (int it = 0; it < m_Mt; ++it) {
      for (int iz = 0; iz < m_Mz; ++iz) {
        for (int iy = 0; iy < m_Ny; ++iy) {
          int iyzt = iy + m_Ny * (iz + m_Nz * it);
          int Leo  = ieo + (1 - 2 * ieo) * m_Leo[iyzt0 + iyzt];
          if (Leo == 1) {
            int is = ix + m_Nx2 * iyzt;
            int in = Nvcd * is;

            int ix1 = Nvc2 * ibf;
            int ix2 = ix1 + m_Nvc;

            for (int ic = 0; ic < m_Nc; ++ic) {
              w2[2 * ic + ix1]     = bc2 * (w1[2 * ic + id1 + in] - w1[2 * ic + 1 + id4 + in]);
              w2[2 * ic + 1 + ix1] = bc2 * (w1[2 * ic + 1 + id1 + in] + w1[2 * ic + id4 + in]);
              w2[2 * ic + ix2]     = bc2 * (w1[2 * ic + id2 + in] - w1[2 * ic + 1 + id3 + in]);
              w2[2 * ic + 1 + ix2] = bc2 * (w1[2 * ic + 1 + id2 + in] + w1[2 * ic + id3 + in]);
            }
            ++ibf;
          }
        }
      }
    }

    m_bw_send[idir]->start_thread(itask);
  }


//====================================================================
  void Fopr_Wilson_eo::mult_xp2_thread(
    int itask, double *v2, const double *vcp2, int ieo)
  {
    int Nvc2  = 2 * m_Nvc;
    int Nvcd  = m_Nvc * m_Nd;
    int Nvcd2 = Nvcd / 2;

    int id1 = 0;
    int id2 = m_Nvc;
    int id3 = m_Nvc * 2;
    int id4 = m_Nvc * 3;

    int idir = 0;

    double wt1r, wt1i, wt2r, wt2i;

    int isite    = m_arg[itask].isite;
    int isite_cp = m_arg[itask].isite_cpx;
    int iyzt0    = isite / m_Nx2;

    double *w2 = &v2[Nvcd * isite];
    //  double* w1 = &vcp2[Nvcd2*isite_cp];
    const double *w1
      = (double *)m_bw_recv[idir]->ptr(sizeof(double) * Nvcd2 * isite_cp);
    const double *u = m_U->ptr(m_Ndf * (isite + ieo * m_Nvol / 2 + idir * m_Nvol));

    m_bw_recv[idir]->wait_thread(itask);

    int ix  = m_Nx2 - 1;
    int ibf = 0;
    for (int it = 0; it < m_Mt; ++it) {
      for (int iz = 0; iz < m_Mz; ++iz) {
        for (int iy = 0; iy < m_Ny; ++iy) {
          int iyzt = iy + m_Ny * (iz + m_Nz * it);
          int Leo  = ieo + (1 - 2 * ieo) * m_Leo[iyzt0 + iyzt];

          if (Leo == 1) {
            int is  = ix + m_Nx2 * iyzt;
            int iv  = Nvcd * is;
            int ig  = m_Ndf * is;
            int ix1 = Nvc2 * ibf;
            int ix2 = ix1 + m_Nvc;

            for (int ic = 0; ic < m_Nc; ++ic) {
              int ic2 = ic * m_Nvc;
              wt1r = mult_uv_r(&u[ic2 + ig], &w1[ix1], m_Nc);
              wt1i = mult_uv_i(&u[ic2 + ig], &w1[ix1], m_Nc);
              wt2r = mult_uv_r(&u[ic2 + ig], &w1[ix2], m_Nc);
              wt2i = mult_uv_i(&u[ic2 + ig], &w1[ix2], m_Nc);
              w2[2 * ic + id1 + iv]     += wt1r;
              w2[2 * ic + 1 + id1 + iv] += wt1i;
              w2[2 * ic + id2 + iv]     += wt2r;
              w2[2 * ic + 1 + id2 + iv] += wt2i;
              w2[2 * ic + id3 + iv]     += wt2i;
              w2[2 * ic + 1 + id3 + iv] += -wt2r;
              w2[2 * ic + id4 + iv]     += wt1i;
              w2[2 * ic + 1 + id4 + iv] += -wt1r;
            }
            ++ibf;
          }
        }
      }
    }
  }


//====================================================================
  void Fopr_Wilson_eo::mult_xpb_thread(
    int itask, double *v2, const double *v1, int ieo)
  {
    int Nvcd = m_Nvc * m_Nd;

    int id1 = 0;
    int id2 = m_Nvc;
    int id3 = m_Nvc * 2;
    int id4 = m_Nvc * 3;

    int idir = 0;

    double vt1[m_Nvc], vt2[m_Nvc];
    double wt1r, wt1i, wt2r, wt2i;

    int isite = m_arg[itask].isite;
    int iyzt0 = isite / m_Nx2;

    double       *w2 = &v2[Nvcd * isite];
    const double *w1 = &v1[Nvcd * isite];
    const double *u  = m_U->ptr(m_Ndf * (isite + ieo * m_Nvol / 2 + idir * m_Nvol));

    for (int it = 0; it < m_Mt; ++it) {
      for (int iz = 0; iz < m_Mz; ++iz) {
        for (int iy = 0; iy < m_Ny; ++iy) {
          int iyzt = iy + m_Ny * (iz + m_Nz * it);
          int Leo  = ieo + (1 - 2 * ieo) * m_Leo[iyzt0 + iyzt];
          for (int ix = 0; ix < m_Nx2 - Leo; ++ix) {
            int is = ix + m_Nx2 * iyzt;
            int iv = Nvcd * is;
            int in = Nvcd * (is + Leo);
            int ig = m_Ndf * is;

            for (int ic = 0; ic < m_Nc; ++ic) {
              vt1[2 * ic]     = w1[2 * ic + id1 + in] - w1[2 * ic + 1 + id4 + in];
              vt1[2 * ic + 1] = w1[2 * ic + 1 + id1 + in] + w1[2 * ic + id4 + in];
              vt2[2 * ic]     = w1[2 * ic + id2 + in] - w1[2 * ic + 1 + id3 + in];
              vt2[2 * ic + 1] = w1[2 * ic + 1 + id2 + in] + w1[2 * ic + id3 + in];
            }

            for (int ic = 0; ic < m_Nc; ++ic) {
              int ic2 = ic * m_Nvc;

              wt1r = mult_uv_r(&u[ic2 + ig], vt1, m_Nc);
              wt1i = mult_uv_i(&u[ic2 + ig], vt1, m_Nc);
              wt2r = mult_uv_r(&u[ic2 + ig], vt2, m_Nc);
              wt2i = mult_uv_i(&u[ic2 + ig], vt2, m_Nc);

              w2[2 * ic + id1 + iv]     += wt1r;
              w2[2 * ic + 1 + id1 + iv] += wt1i;
              w2[2 * ic + id2 + iv]     += wt2r;
              w2[2 * ic + 1 + id2 + iv] += wt2i;
              w2[2 * ic + id3 + iv]     += wt2i;
              w2[2 * ic + 1 + id3 + iv] += -wt2r;
              w2[2 * ic + id4 + iv]     += wt1i;
              w2[2 * ic + 1 + id4 + iv] += -wt1r;
            }
          }
        }
      }
    }
  }


//====================================================================
  void Fopr_Wilson_eo::mult_xm1_thread(
    int itask, double *vcp1, const double *v1, int ieo)
  {
    int Nvc2  = 2 * m_Nvc;
    int Nvcd  = m_Nvc * m_Nd;
    int Nvcd2 = Nvcd / 2;

    int id1 = 0;
    int id2 = m_Nvc;
    int id3 = m_Nvc * 2;
    int id4 = m_Nvc * 3;

    int idir = 0;

    int isite    = m_arg[itask].isite;
    int isite_cp = m_arg[itask].isite_cpx;
    int iyzt0    = isite / m_Nx2;

    //  double* w2 = &vcp1[Nvcd2*isite_cp];
    double *w2
      = (double *)m_fw_send[idir]->ptr(sizeof(double) * Nvcd2 * isite_cp);
    const double *w1 = &v1[Nvcd * isite];
    const double *u  = m_U->ptr(m_Ndf * (isite + (1 - ieo) * m_Nvol / 2 + idir * m_Nvol));

    double vt1[m_Nvc], vt2[m_Nvc];

    int ix  = m_Nx2 - 1;
    int ibf = 0;

    for (int it = 0; it < m_Mt; ++it) {
      for (int iz = 0; iz < m_Mz; ++iz) {
        for (int iy = 0; iy < m_Ny; ++iy) {
          int iyzt = iy + m_Ny * (iz + m_Nz * it);
          int Leo  = ieo + (1 - 2 * ieo) * m_Leo[iyzt0 + iyzt];
          if (Leo == 0) {
            int is = ix + m_Nx2 * iyzt;
            int in = Nvcd * is;
            int ig = m_Ndf * is;

            int ix1 = Nvc2 * ibf;
            int ix2 = ix1 + m_Nvc;

            for (int ic = 0; ic < m_Nc; ++ic) {
              vt1[2 * ic]     = w1[2 * ic + id1 + in] + w1[2 * ic + 1 + id4 + in];
              vt1[2 * ic + 1] = w1[2 * ic + 1 + id1 + in] - w1[2 * ic + id4 + in];
              vt2[2 * ic]     = w1[2 * ic + id2 + in] + w1[2 * ic + 1 + id3 + in];
              vt2[2 * ic + 1] = w1[2 * ic + 1 + id2 + in] - w1[2 * ic + id3 + in];
            }

            for (int ic = 0; ic < m_Nc; ++ic) {
              int icr = 2 * ic;
              w2[icr + ix1]     = mult_udagv_r(&u[icr + ig], vt1, m_Nc);
              w2[icr + 1 + ix1] = mult_udagv_i(&u[icr + ig], vt1, m_Nc);
              w2[icr + ix2]     = mult_udagv_r(&u[icr + ig], vt2, m_Nc);
              w2[icr + 1 + ix2] = mult_udagv_i(&u[icr + ig], vt2, m_Nc);
            }
            ++ibf;
          }
        }
      }
    }

    m_fw_send[idir]->start_thread(itask);
  }


//====================================================================
  void Fopr_Wilson_eo::mult_xm2_thread(
    int itask, double *v2, const double *vcp2, int ieo)
  {
    int Nvc2  = 2 * m_Nvc;
    int Nvcd  = m_Nvc * m_Nd;
    int Nvcd2 = Nvcd / 2;

    int id1 = 0;
    int id2 = m_Nvc;
    int id3 = m_Nvc * 2;
    int id4 = m_Nvc * 3;

    int    idir = 0;
    double bc2  = m_boundary2[idir];

    double wt1r, wt1i, wt2r, wt2i;

    int isite    = m_arg[itask].isite;
    int isite_cp = m_arg[itask].isite_cpx;
    int iyzt0    = isite / m_Nx2;

    double *w2 = &v2[Nvcd * isite];
    //  double* w1 = &vcp2[Nvcd2*isite_cp];
    const double *w1
      = (double *)m_fw_recv[idir]->ptr(sizeof(double) * Nvcd2 * isite_cp);

    m_fw_recv[idir]->wait_thread(itask);

    int ix  = 0;
    int ibf = 0;
    for (int it = 0; it < m_Mt; ++it) {
      for (int iz = 0; iz < m_Mz; ++iz) {
        for (int iy = 0; iy < m_Ny; ++iy) {
          int iyzt = iy + m_Ny * (iz + m_Nz * it);
          int Leo  = ieo + (1 - 2 * ieo) * m_Leo[iyzt0 + iyzt];
          if (Leo == 0) {
            int is = ix + m_Nx2 * iyzt;
            int iv = Nvcd * is;

            int ix1 = Nvc2 * ibf;
            int ix2 = ix1 + m_Nvc;

            for (int ic = 0; ic < m_Nc; ++ic) {
              int icr = 2 * ic;
              int ici = 2 * ic + 1;
              w2[icr + id1 + iv] += bc2 * w1[icr + ix1];
              w2[ici + id1 + iv] += bc2 * w1[ici + ix1];
              w2[icr + id2 + iv] += bc2 * w1[icr + ix2];
              w2[ici + id2 + iv] += bc2 * w1[ici + ix2];
              w2[icr + id3 + iv] += -bc2 * w1[ici + ix2];
              w2[ici + id3 + iv] += +bc2 * w1[icr + ix2];
              w2[icr + id4 + iv] += -bc2 * w1[ici + ix1];
              w2[ici + id4 + iv] += +bc2 * w1[icr + ix1];
            }
            ++ibf;
          }
        }
      }
    }
  }


//====================================================================
  void Fopr_Wilson_eo::mult_xmb_thread(
    int itask, double *v2, const double *v1, int ieo)
  {
    int Nvcd = m_Nvc * m_Nd;

    int id1 = 0;
    int id2 = m_Nvc;
    int id3 = m_Nvc * 2;
    int id4 = m_Nvc * 3;

    int idir = 0;

    double vt1[m_Nvc], vt2[m_Nvc];
    double wt1r, wt1i, wt2r, wt2i;

    int isite = m_arg[itask].isite;
    int iyzt0 = isite / m_Nx2;

    double       *w2 = &v2[Nvcd * isite];
    const double *w1 = &v1[Nvcd * isite];
    const double *u  = m_U->ptr(m_Ndf * (isite + (1 - ieo) * m_Nvol / 2 + idir * m_Nvol));

    for (int it = 0; it < m_Mt; ++it) {
      for (int iz = 0; iz < m_Mz; ++iz) {
        for (int iy = 0; iy < m_Ny; ++iy) {
          int iyzt = iy + m_Ny * (iz + m_Nz * it);
          int Leo  = ieo + (1 - 2 * ieo) * m_Leo[iyzt0 + iyzt];
          int Meo  = 1 - Leo;
          for (int ix = Meo; ix < m_Nx2; ++ix) {
            int is = ix + m_Nx2 * iyzt;
            int iv = Nvcd * is;
            int in = Nvcd * (is - Meo);
            int ig = m_Ndf * (is - Meo);

            for (int ic = 0; ic < m_Nc; ++ic) {
              vt1[2 * ic]     = w1[2 * ic + id1 + in] + w1[2 * ic + 1 + id4 + in];
              vt1[2 * ic + 1] = w1[2 * ic + 1 + id1 + in] - w1[2 * ic + id4 + in];
              vt2[2 * ic]     = w1[2 * ic + id2 + in] + w1[2 * ic + 1 + id3 + in];
              vt2[2 * ic + 1] = w1[2 * ic + 1 + id2 + in] - w1[2 * ic + id3 + in];
            }

            for (int ic = 0; ic < m_Nc; ++ic) {
              int ic2 = 2 * ic;

              wt1r = mult_udagv_r(&u[ic2 + ig], vt1, m_Nc);
              wt1i = mult_udagv_i(&u[ic2 + ig], vt1, m_Nc);
              wt2r = mult_udagv_r(&u[ic2 + ig], vt2, m_Nc);
              wt2i = mult_udagv_i(&u[ic2 + ig], vt2, m_Nc);

              w2[2 * ic + id1 + iv]     += wt1r;
              w2[2 * ic + 1 + id1 + iv] += wt1i;
              w2[2 * ic + id2 + iv]     += wt2r;
              w2[2 * ic + 1 + id2 + iv] += wt2i;
              w2[2 * ic + id3 + iv]     += -wt2i;
              w2[2 * ic + 1 + id3 + iv] += +wt2r;
              w2[2 * ic + id4 + iv]     += -wt1i;
              w2[2 * ic + 1 + id4 + iv] += +wt1r;
            }
          }
        }
      }
    }
  }


//====================================================================
  void Fopr_Wilson_eo::mult_yp1_thread(
    int itask, double *vcp1, const double *v1, int ieo)
  {
    int Nvc2  = 2 * m_Nvc;
    int Nvcd  = m_Nvc * m_Nd;
    int Nvcd2 = Nvcd / 2;

    int id1 = 0;
    int id2 = m_Nvc;
    int id3 = m_Nvc * 2;
    int id4 = m_Nvc * 3;

    int idir = 1;

    int isite    = m_arg[itask].isite;
    int isite_cp = m_arg[itask].isite_cpy;

    //  double* w2 = &vcp1[Nvcd2*isite_cp];
    double *w2
      = (double *)m_bw_send[idir]->ptr(sizeof(double) * Nvcd2 * isite_cp);
    const double *w1 = &v1[Nvcd * isite];

    double bc2 = m_boundary2[idir];

    int iy = 0;

    for (int it = 0; it < m_Mt; ++it) {
      for (int iz = 0; iz < m_Mz; ++iz) {
        for (int ix = 0; ix < m_Nx2; ++ix) {
          int is  = ix + m_Nx2 * (iy + m_Ny * (iz + m_Nz * it));
          int is2 = ix + m_Nx2 * (iz + m_Mz * it);
          int in  = Nvcd * is;
          int ix1 = Nvc2 * is2;
          int ix2 = ix1 + m_Nvc;

          for (int ic = 0; ic < m_Nc; ++ic) {
            w2[2 * ic + ix1]     = bc2 * (w1[2 * ic + id1 + in] + w1[2 * ic + id4 + in]);
            w2[2 * ic + 1 + ix1] = bc2 * (w1[2 * ic + 1 + id1 + in] + w1[2 * ic + 1 + id4 + in]);
            w2[2 * ic + ix2]     = bc2 * (w1[2 * ic + id2 + in] - w1[2 * ic + id3 + in]);
            w2[2 * ic + 1 + ix2] = bc2 * (w1[2 * ic + 1 + id2 + in] - w1[2 * ic + 1 + id3 + in]);
          }
        }
      }
    }

    m_bw_send[idir]->start_thread(itask);
  }


//====================================================================
  void Fopr_Wilson_eo::mult_yp2_thread(
    int itask, double *v2, const double *vcp2, int ieo)
  {
    int Nvc2  = 2 * m_Nvc;
    int Nvcd  = m_Nvc * m_Nd;
    int Nvcd2 = Nvcd / 2;

    int id1 = 0;
    int id2 = m_Nvc;
    int id3 = m_Nvc * 2;
    int id4 = m_Nvc * 3;

    int idir = 1;

    double wt1r, wt1i, wt2r, wt2i;

    int isite    = m_arg[itask].isite;
    int isite_cp = m_arg[itask].isite_cpy;

    double *w2 = &v2[Nvcd * isite];
    //  double* w1 = &vcp2[Nvcd2*isite_cp];
    const double *w1
      = (double *)m_bw_recv[idir]->ptr(sizeof(double) * Nvcd2 * isite_cp);
    const double *u = m_U->ptr(m_Ndf * (isite + ieo * m_Nvol / 2 + idir * m_Nvol));

    m_bw_recv[idir]->wait_thread(itask);

    int iy = m_Ny - 1;
    for (int it = 0; it < m_Mt; ++it) {
      for (int iz = 0; iz < m_Mz; ++iz) {
        for (int ix = 0; ix < m_Nx2; ++ix) {
          int is  = ix + m_Nx2 * (iy + m_Ny * (iz + m_Nz * it));
          int is2 = ix + m_Nx2 * (iz + m_Mz * it);
          int iv  = Nvcd * is;
          int ig  = m_Ndf * is;
          int ix1 = Nvc2 * is2;
          int ix2 = ix1 + m_Nvc;

          for (int ic = 0; ic < m_Nc; ++ic) {
            int ic2 = ic * m_Nvc;

            wt1r = mult_uv_r(&u[ic2 + ig], &w1[ix1], m_Nc);
            wt1i = mult_uv_i(&u[ic2 + ig], &w1[ix1], m_Nc);
            wt2r = mult_uv_r(&u[ic2 + ig], &w1[ix2], m_Nc);
            wt2i = mult_uv_i(&u[ic2 + ig], &w1[ix2], m_Nc);

            w2[2 * ic + id1 + iv]     += wt1r;
            w2[2 * ic + 1 + id1 + iv] += wt1i;
            w2[2 * ic + id2 + iv]     += wt2r;
            w2[2 * ic + 1 + id2 + iv] += wt2i;
            w2[2 * ic + id3 + iv]     += -wt2r;
            w2[2 * ic + 1 + id3 + iv] += -wt2i;
            w2[2 * ic + id4 + iv]     += wt1r;
            w2[2 * ic + 1 + id4 + iv] += wt1i;
          }
        }
      }
    }
  }


//====================================================================
  void Fopr_Wilson_eo::mult_ypb_thread(
    int itask, double *v2, const double *v1, int ieo)
  {
    int Nvcd = m_Nvc * m_Nd;

    int id1 = 0;
    int id2 = m_Nvc;
    int id3 = m_Nvc * 2;
    int id4 = m_Nvc * 3;

    int idir = 1;

    double vt1[m_Nvc], vt2[m_Nvc];
    double wt1r, wt1i, wt2r, wt2i;

    int isite = m_arg[itask].isite;

    double       *w2 = &v2[Nvcd * isite];
    const double *w1 = &v1[Nvcd * isite];
    const double *u  = m_U->ptr(m_Ndf * (isite + ieo * m_Nvol / 2 + idir * m_Nvol));

    for (int it = 0; it < m_Mt; ++it) {
      for (int iz = 0; iz < m_Mz; ++iz) {
        for (int iy = 0; iy < m_Ny - 1; ++iy) {
          for (int ix = 0; ix < m_Nx2; ++ix) {
            int is = ix + m_Nx2 * (iy + m_Ny * (iz + m_Nz * it));
            int iv = Nvcd * is;
            int in = Nvcd * (is + m_Nx2);
            int ig = m_Ndf * is;

            for (int ic = 0; ic < m_Nc; ++ic) {
              vt1[2 * ic]     = w1[2 * ic + id1 + in] + w1[2 * ic + id4 + in];
              vt1[2 * ic + 1] = w1[2 * ic + 1 + id1 + in] + w1[2 * ic + 1 + id4 + in];
              vt2[2 * ic]     = w1[2 * ic + id2 + in] - w1[2 * ic + id3 + in];
              vt2[2 * ic + 1] = w1[2 * ic + 1 + id2 + in] - w1[2 * ic + 1 + id3 + in];
            }

            for (int ic = 0; ic < m_Nc; ++ic) {
              int ic2 = ic * m_Nvc;

              wt1r = mult_uv_r(&u[ic2 + ig], vt1, m_Nc);
              wt1i = mult_uv_i(&u[ic2 + ig], vt1, m_Nc);
              wt2r = mult_uv_r(&u[ic2 + ig], vt2, m_Nc);
              wt2i = mult_uv_i(&u[ic2 + ig], vt2, m_Nc);

              w2[2 * ic + id1 + iv]     += wt1r;
              w2[2 * ic + 1 + id1 + iv] += wt1i;
              w2[2 * ic + id2 + iv]     += wt2r;
              w2[2 * ic + 1 + id2 + iv] += wt2i;
              w2[2 * ic + id3 + iv]     += -wt2r;
              w2[2 * ic + 1 + id3 + iv] += -wt2i;
              w2[2 * ic + id4 + iv]     += wt1r;
              w2[2 * ic + 1 + id4 + iv] += wt1i;
            }
          }
        }
      }
    }
  }


//====================================================================
  void Fopr_Wilson_eo::mult_ym1_thread(
    int itask, double *vcp1, const double *v1, int ieo)
  {
    int Nvc2  = 2 * m_Nvc;
    int Nvcd  = m_Nvc * m_Nd;
    int Nvcd2 = Nvcd / 2;

    int id1 = 0;
    int id2 = m_Nvc;
    int id3 = m_Nvc * 2;
    int id4 = m_Nvc * 3;

    int idir = 1;

    int isite    = m_arg[itask].isite;
    int isite_cp = m_arg[itask].isite_cpy;

    //  double* w2 = &vcp1[Nvcd2*isite_cp];
    double *w2
      = (double *)m_fw_send[idir]->ptr(sizeof(double) * Nvcd2 * isite_cp);
    const double *w1 = &v1[Nvcd * isite];
    const double *u  = m_U->ptr(m_Ndf * (isite + (1 - ieo) * m_Nvol / 2 + idir * m_Nvol));

    double vt1[m_Nvc], vt2[m_Nvc];

    int iy = m_Ny - 1;

    for (int it = 0; it < m_Mt; ++it) {
      for (int iz = 0; iz < m_Mz; ++iz) {
        for (int ix = 0; ix < m_Nx2; ++ix) {
          int is  = ix + m_Nx2 * (iy + m_Ny * (iz + m_Nz * it));
          int is2 = ix + m_Nx2 * (iz + m_Mz * it);
          int in  = Nvcd * is;
          int ig  = m_Ndf * is;
          int ix1 = Nvc2 * is2;
          int ix2 = ix1 + m_Nvc;

          for (int ic = 0; ic < m_Nc; ++ic) {
            vt1[2 * ic]     = w1[2 * ic + id1 + in] - w1[2 * ic + id4 + in];
            vt1[2 * ic + 1] = w1[2 * ic + 1 + id1 + in] - w1[2 * ic + 1 + id4 + in];
            vt2[2 * ic]     = w1[2 * ic + id2 + in] + w1[2 * ic + id3 + in];
            vt2[2 * ic + 1] = w1[2 * ic + 1 + id2 + in] + w1[2 * ic + 1 + id3 + in];
          }

          for (int ic = 0; ic < m_Nc; ++ic) {
            int icr = 2 * ic;
            w2[icr + ix1]     = mult_udagv_r(&u[icr + ig], vt1, m_Nc);
            w2[icr + 1 + ix1] = mult_udagv_i(&u[icr + ig], vt1, m_Nc);
            w2[icr + ix2]     = mult_udagv_r(&u[icr + ig], vt2, m_Nc);
            w2[icr + 1 + ix2] = mult_udagv_i(&u[icr + ig], vt2, m_Nc);
          }
        }
      }
    }

    m_fw_send[idir]->start_thread(itask);
  }


//====================================================================
  void Fopr_Wilson_eo::mult_ym2_thread(
    int itask, double *v2, const double *vcp2, int ieo)
  {
    int Nvc2  = 2 * m_Nvc;
    int Nvcd  = m_Nvc * m_Nd;
    int Nvcd2 = Nvcd / 2;

    int id1 = 0;
    int id2 = m_Nvc;
    int id3 = m_Nvc * 2;
    int id4 = m_Nvc * 3;

    int    idir = 1;
    double bc2  = m_boundary2[idir];

    double wt1r, wt1i, wt2r, wt2i;

    int isite    = m_arg[itask].isite;
    int isite_cp = m_arg[itask].isite_cpy;

    double *w2 = &v2[Nvcd * isite];
    //  double* w1 = &vcp2[Nvcd2*isite_cp];
    const double *w1
      = (double *)m_fw_recv[idir]->ptr(sizeof(double) * Nvcd2 * isite_cp);

    m_fw_recv[idir]->wait_thread(itask);

    int iy = 0;
    for (int it = 0; it < m_Mt; ++it) {
      for (int iz = 0; iz < m_Mz; ++iz) {
        for (int ix = 0; ix < m_Nx2; ++ix) {
          int is  = ix + m_Nx2 * (iy + m_Ny * (iz + m_Nz * it));
          int is2 = ix + m_Nx2 * (iz + m_Mz * it);
          int iv  = Nvcd * is;
          int ix1 = Nvc2 * is2;
          int ix2 = ix1 + m_Nvc;

          for (int ic = 0; ic < m_Nc; ++ic) {
            int icr = 2 * ic;
            int ici = 2 * ic + 1;
            w2[icr + id1 + iv] += bc2 * w1[icr + ix1];
            w2[ici + id1 + iv] += bc2 * w1[ici + ix1];
            w2[icr + id2 + iv] += bc2 * w1[icr + ix2];
            w2[ici + id2 + iv] += bc2 * w1[ici + ix2];
            w2[icr + id3 + iv] += bc2 * w1[icr + ix2];
            w2[ici + id3 + iv] += bc2 * w1[ici + ix2];
            w2[icr + id4 + iv] += -bc2 * w1[icr + ix1];
            w2[ici + id4 + iv] += -bc2 * w1[ici + ix1];
          }
        }
      }
    }
  }


//====================================================================
  void Fopr_Wilson_eo::mult_ymb_thread(
    int itask, double *v2, const double *v1, int ieo)
  {
    int Nvcd = m_Nvc * m_Nd;

    int id1 = 0;
    int id2 = m_Nvc;
    int id3 = m_Nvc * 2;
    int id4 = m_Nvc * 3;

    int idir = 1;

    double vt1[m_Nvc], vt2[m_Nvc];
    double wt1r, wt1i, wt2r, wt2i;

    int isite = m_arg[itask].isite;

    double       *w2 = &v2[Nvcd * isite];
    const double *w1 = &v1[Nvcd * isite];
    const double *u  = m_U->ptr(m_Ndf * (isite + (1 - ieo) * m_Nvol / 2 + idir * m_Nvol));

    for (int it = 0; it < m_Mt; ++it) {
      for (int iz = 0; iz < m_Mz; ++iz) {
        for (int iy = 1; iy < m_Ny; ++iy) {
          for (int ix = 0; ix < m_Nx2; ++ix) {
            int is = ix + m_Nx2 * (iy + m_Ny * (iz + m_Nz * it));
            int iv = Nvcd * is;
            int in = Nvcd * (is - m_Nx2);
            int ig = m_Ndf * (is - m_Nx2);

            for (int ic = 0; ic < m_Nc; ++ic) {
              vt1[2 * ic]     = w1[2 * ic + id1 + in] - w1[2 * ic + id4 + in];
              vt1[2 * ic + 1] = w1[2 * ic + 1 + id1 + in] - w1[2 * ic + 1 + id4 + in];
              vt2[2 * ic]     = w1[2 * ic + id2 + in] + w1[2 * ic + id3 + in];
              vt2[2 * ic + 1] = w1[2 * ic + 1 + id2 + in] + w1[2 * ic + 1 + id3 + in];
            }

            for (int ic = 0; ic < m_Nc; ++ic) {
              int ic2 = 2 * ic;
              wt1r = mult_udagv_r(&u[ic2 + ig], vt1, m_Nc);
              wt1i = mult_udagv_i(&u[ic2 + ig], vt1, m_Nc);
              wt2r = mult_udagv_r(&u[ic2 + ig], vt2, m_Nc);
              wt2i = mult_udagv_i(&u[ic2 + ig], vt2, m_Nc);

              w2[ic2 + id1 + iv]     += wt1r;
              w2[ic2 + 1 + id1 + iv] += wt1i;
              w2[ic2 + id2 + iv]     += wt2r;
              w2[ic2 + 1 + id2 + iv] += wt2i;
              w2[ic2 + id3 + iv]     += wt2r;
              w2[ic2 + 1 + id3 + iv] += wt2i;
              w2[ic2 + id4 + iv]     += -wt1r;
              w2[ic2 + 1 + id4 + iv] += -wt1i;
            }
          }
        }
      }
    }
  }


//====================================================================
  void Fopr_Wilson_eo::mult_zp1_thread(
    int itask, double *vcp1, const double *v1, int ieo)
  {
    int Nvc2  = 2 * m_Nvc;
    int Nvcd  = m_Nvc * m_Nd;
    int Nvcd2 = Nvcd / 2;

    int id1 = 0;
    int id2 = m_Nvc;
    int id3 = m_Nvc * 2;
    int id4 = m_Nvc * 3;

    int idir = 2;

    int isite    = m_arg[itask].isite;
    int isite_cp = m_arg[itask].isite_cpz;

    //  double* w2 = &vcp1[Nvcd2*isite_cp];
    double *w2
      = (double *)m_bw_send[idir]->ptr(sizeof(double) * Nvcd2 * isite_cp);
    const double *w1 = &v1[Nvcd * isite];

    double bc2 = m_boundary2[idir];

    if (m_arg[itask].kz0 == 1) {
      int Nxy = m_Nx2 * m_Ny;
      int iz  = 0;
      for (int it = 0; it < m_Mt; ++it) {
        for (int ixy = 0; ixy < Nxy; ++ixy) {
          int is  = ixy + Nxy * (iz + m_Nz * it);
          int is2 = ixy + Nxy * it;

          int in  = Nvcd * is;
          int ix1 = Nvc2 * is2;
          int ix2 = ix1 + m_Nvc;

          for (int ic = 0; ic < m_Nc; ++ic) {
            w2[2 * ic + ix1]     = bc2 * (w1[2 * ic + id1 + in] - w1[2 * ic + 1 + id3 + in]);
            w2[2 * ic + 1 + ix1] = bc2 * (w1[2 * ic + 1 + id1 + in] + w1[2 * ic + id3 + in]);
            w2[2 * ic + ix2]     = bc2 * (w1[2 * ic + id2 + in] + w1[2 * ic + 1 + id4 + in]);
            w2[2 * ic + 1 + ix2] = bc2 * (w1[2 * ic + 1 + id2 + in] - w1[2 * ic + id4 + in]);
          }
        }
      }
    }

    m_bw_send[idir]->start_thread(itask);
  }


//====================================================================
  void Fopr_Wilson_eo::mult_zp2_thread(
    int itask, double *v2, const double *vcp2, int ieo)
  {
    int Nvc2  = 2 * m_Nvc;
    int Nvcd  = m_Nvc * m_Nd;
    int Nvcd2 = Nvcd / 2;

    int id1 = 0;
    int id2 = m_Nvc;
    int id3 = m_Nvc * 2;
    int id4 = m_Nvc * 3;

    int idir = 2;

    double wt1r, wt1i, wt2r, wt2i;

    int isite    = m_arg[itask].isite;
    int isite_cp = m_arg[itask].isite_cpz;

    double *w2 = &v2[Nvcd * isite];
    //  double* w1 = &vcp2[Nvcd2*isite_cp];
    const double *w1
      = (double *)m_bw_recv[idir]->ptr(sizeof(double) * Nvcd2 * isite_cp);
    const double *u = m_U->ptr(m_Ndf * (isite + ieo * m_Nvol / 2 + idir * m_Nvol));

    m_bw_recv[idir]->wait_thread(itask);

    if (m_arg[itask].kz1 == 1) {
      int Nxy = m_Nx2 * m_Ny;
      int iz  = m_Mz - 1;
      for (int it = 0; it < m_Mt; ++it) {
        for (int ixy = 0; ixy < Nxy; ++ixy) {
          int is  = ixy + Nxy * (iz + m_Nz * it);
          int is2 = ixy + Nxy * it;
          int iv  = Nvcd * is;
          int ig  = m_Ndf * is;
          int ix1 = Nvc2 * is2;
          int ix2 = ix1 + m_Nvc;

          for (int ic = 0; ic < m_Nc; ++ic) {
            int ic2 = ic * m_Nvc;

            wt1r = mult_uv_r(&u[ic2 + ig], &w1[ix1], m_Nc);
            wt1i = mult_uv_i(&u[ic2 + ig], &w1[ix1], m_Nc);
            wt2r = mult_uv_r(&u[ic2 + ig], &w1[ix2], m_Nc);
            wt2i = mult_uv_i(&u[ic2 + ig], &w1[ix2], m_Nc);

            w2[2 * ic + id1 + iv]     += wt1r;
            w2[2 * ic + 1 + id1 + iv] += wt1i;
            w2[2 * ic + id2 + iv]     += wt2r;
            w2[2 * ic + 1 + id2 + iv] += wt2i;
            w2[2 * ic + id3 + iv]     += wt1i;
            w2[2 * ic + 1 + id3 + iv] += -wt1r;
            w2[2 * ic + id4 + iv]     += -wt2i;
            w2[2 * ic + 1 + id4 + iv] += wt2r;
          }
        }
      }
    }
  }


//====================================================================
  void Fopr_Wilson_eo::mult_zpb_thread(
    int itask, double *v2, const double *v1, int ieo)
  {
    int Nvcd = m_Nvc * m_Nd;

    int id1 = 0;
    int id2 = m_Nvc;
    int id3 = m_Nvc * 2;
    int id4 = m_Nvc * 3;

    int idir = 2;

    double vt1[m_Nvc], vt2[m_Nvc];
    double wt1r, wt1i, wt2r, wt2i;

    int isite = m_arg[itask].isite;

    double       *w2 = &v2[Nvcd * isite];
    const double *w1 = &v1[Nvcd * isite];
    const double *u  = m_U->ptr(m_Ndf * (isite + ieo * m_Nvol / 2 + idir * m_Nvol));

    int kz1 = m_arg[itask].kz1;
    int Nxy = m_Nx2 * m_Ny;

    for (int it = 0; it < m_Mt; ++it) {
      for (int iz = 0; iz < m_Mz - kz1; ++iz) {
        for (int ixy = 0; ixy < Nxy; ++ixy) {
          int is = ixy + Nxy * (iz + m_Nz * it);
          int iv = Nvcd * is;
          int in = Nvcd * (is + Nxy);
          int ig = m_Ndf * is;

          for (int ic = 0; ic < m_Nc; ++ic) {
            vt1[2 * ic]     = w1[2 * ic + id1 + in] - w1[2 * ic + 1 + id3 + in];
            vt1[2 * ic + 1] = w1[2 * ic + 1 + id1 + in] + w1[2 * ic + id3 + in];
            vt2[2 * ic]     = w1[2 * ic + id2 + in] + w1[2 * ic + 1 + id4 + in];
            vt2[2 * ic + 1] = w1[2 * ic + 1 + id2 + in] - w1[2 * ic + id4 + in];
          }

          for (int ic = 0; ic < m_Nc; ++ic) {
            int ic2 = ic * m_Nvc;

            wt1r = mult_uv_r(&u[ic2 + ig], vt1, m_Nc);
            wt1i = mult_uv_i(&u[ic2 + ig], vt1, m_Nc);
            wt2r = mult_uv_r(&u[ic2 + ig], vt2, m_Nc);
            wt2i = mult_uv_i(&u[ic2 + ig], vt2, m_Nc);

            w2[2 * ic + id1 + iv]     += wt1r;
            w2[2 * ic + 1 + id1 + iv] += wt1i;
            w2[2 * ic + id2 + iv]     += wt2r;
            w2[2 * ic + 1 + id2 + iv] += wt2i;
            w2[2 * ic + id3 + iv]     += wt1i;
            w2[2 * ic + 1 + id3 + iv] += -wt1r;
            w2[2 * ic + id4 + iv]     += -wt2i;
            w2[2 * ic + 1 + id4 + iv] += wt2r;
          }
        }
      }
    }
  }


//====================================================================
  void Fopr_Wilson_eo::mult_zm1_thread(
    int itask, double *vcp1, const double *v1, int ieo)
  {
    int Nvc2  = 2 * m_Nvc;
    int Nvcd  = m_Nvc * m_Nd;
    int Nvcd2 = Nvcd / 2;

    int id1 = 0;
    int id2 = m_Nvc;
    int id3 = m_Nvc * 2;
    int id4 = m_Nvc * 3;

    int idir = 2;

    int isite    = m_arg[itask].isite;
    int isite_cp = m_arg[itask].isite_cpz;

    //  double* w2 = &vcp1[Nvcd2*isite_cp];
    double *w2
      = (double *)m_fw_send[idir]->ptr(sizeof(double) * Nvcd2 * isite_cp);
    const double *w1 = &v1[Nvcd * isite];
    const double *u  = m_U->ptr(m_Ndf * (isite + (1 - ieo) * m_Nvol / 2 + idir * m_Nvol));

    double vt1[m_Nvc], vt2[m_Nvc];

    if (m_arg[itask].kz1 == 1) {
      int Nxy = m_Nx2 * m_Ny;
      int iz  = m_Mz - 1;
      for (int it = 0; it < m_Mt; ++it) {
        for (int ixy = 0; ixy < Nxy; ++ixy) {
          int is  = ixy + Nxy * (iz + m_Nz * it);
          int is2 = ixy + Nxy * it;
          int in  = Nvcd * is;
          int ig  = m_Ndf * is;
          int ix1 = Nvc2 * is2;
          int ix2 = ix1 + m_Nvc;

          for (int ic = 0; ic < m_Nc; ++ic) {
            vt1[2 * ic]     = w1[2 * ic + id1 + in] + w1[2 * ic + 1 + id3 + in];
            vt1[2 * ic + 1] = w1[2 * ic + 1 + id1 + in] - w1[2 * ic + id3 + in];
            vt2[2 * ic]     = w1[2 * ic + id2 + in] - w1[2 * ic + 1 + id4 + in];
            vt2[2 * ic + 1] = w1[2 * ic + 1 + id2 + in] + w1[2 * ic + id4 + in];
          }

          for (int ic = 0; ic < m_Nc; ++ic) {
            int icr = 2 * ic;
            w2[icr + ix1]     = mult_udagv_r(&u[icr + ig], vt1, m_Nc);
            w2[icr + 1 + ix1] = mult_udagv_i(&u[icr + ig], vt1, m_Nc);
            w2[icr + ix2]     = mult_udagv_r(&u[icr + ig], vt2, m_Nc);
            w2[icr + 1 + ix2] = mult_udagv_i(&u[icr + ig], vt2, m_Nc);
          }
        }
      }
    }

    m_fw_send[idir]->start_thread(itask);
  }


//====================================================================
  void Fopr_Wilson_eo::mult_zm2_thread(
    int itask, double *v2, const double *vcp2, int ieo)
  {
    int Nvc2  = 2 * m_Nvc;
    int Nvcd  = m_Nvc * m_Nd;
    int Nvcd2 = Nvcd / 2;

    int id1 = 0;
    int id2 = m_Nvc;
    int id3 = m_Nvc * 2;
    int id4 = m_Nvc * 3;

    int    idir = 2;
    double bc2  = m_boundary2[idir];

    double wt1r, wt1i, wt2r, wt2i;

    int isite    = m_arg[itask].isite;
    int isite_cp = m_arg[itask].isite_cpz;

    double *w2 = &v2[Nvcd * isite];
    //  double* w1 = &vcp2[Nvcd2*isite_cp];
    const double *w1
      = (double *)m_fw_recv[idir]->ptr(sizeof(double) * Nvcd2 * isite_cp);

    m_fw_recv[idir]->wait_thread(itask);

    if (m_arg[itask].kz0 == 1) {
      int Nxy = m_Nx2 * m_Ny;

      int iz = 0;
      for (int it = 0; it < m_Mt; ++it) {
        for (int ixy = 0; ixy < Nxy; ++ixy) {
          int is  = ixy + Nxy * (iz + m_Nz * it);
          int is2 = ixy + Nxy * it;
          int iv  = Nvcd * is;
          int ix1 = Nvc2 * is2;
          int ix2 = ix1 + m_Nvc;

          for (int ic = 0; ic < m_Nc; ++ic) {
            int icr = 2 * ic;
            int ici = 2 * ic + 1;
            w2[icr + id1 + iv] += bc2 * w1[icr + ix1];
            w2[ici + id1 + iv] += bc2 * w1[ici + ix1];
            w2[icr + id2 + iv] += bc2 * w1[icr + ix2];
            w2[ici + id2 + iv] += bc2 * w1[ici + ix2];
            w2[icr + id3 + iv] += -bc2 * w1[ici + ix1];
            w2[ici + id3 + iv] += bc2 * w1[icr + ix1];
            w2[icr + id4 + iv] += bc2 * w1[ici + ix2];
            w2[ici + id4 + iv] += -bc2 * w1[icr + ix2];
          }
        }
      }
    }
  }


//====================================================================
  void Fopr_Wilson_eo::mult_zmb_thread(
    int itask, double *v2, const double *v1, int ieo)
  {
    int Nvcd = m_Nvc * m_Nd;

    int id1 = 0;
    int id2 = m_Nvc;
    int id3 = m_Nvc * 2;
    int id4 = m_Nvc * 3;

    int idir = 2;

    double vt1[m_Nvc], vt2[m_Nvc];
    double wt1r, wt1i, wt2r, wt2i;

    int isite = m_arg[itask].isite;

    double       *w2 = &v2[Nvcd * isite];
    const double *w1 = &v1[Nvcd * isite];
    const double *u  = m_U->ptr(m_Ndf * (isite + (1 - ieo) * m_Nvol / 2 + idir * m_Nvol));

    int kz0 = m_arg[itask].kz0;
    int Nxy = m_Nx2 * m_Ny;

    for (int it = 0; it < m_Mt; ++it) {
      for (int iz = kz0; iz < m_Mz; ++iz) {
        for (int ixy = 0; ixy < Nxy; ++ixy) {
          int is = ixy + Nxy * (iz + m_Nz * it);
          int iv = Nvcd * is;
          int in = Nvcd * (is - Nxy);
          int ig = m_Ndf * (is - Nxy);

          for (int ic = 0; ic < m_Nc; ++ic) {
            vt1[2 * ic]     = w1[2 * ic + id1 + in] + w1[2 * ic + 1 + id3 + in];
            vt1[2 * ic + 1] = w1[2 * ic + 1 + id1 + in] - w1[2 * ic + id3 + in];
            vt2[2 * ic]     = w1[2 * ic + id2 + in] - w1[2 * ic + 1 + id4 + in];
            vt2[2 * ic + 1] = w1[2 * ic + 1 + id2 + in] + w1[2 * ic + id4 + in];
          }

          for (int ic = 0; ic < m_Nc; ++ic) {
            int ic2 = 2 * ic;
            wt1r = mult_udagv_r(&u[ic2 + ig], vt1, m_Nc);
            wt1i = mult_udagv_i(&u[ic2 + ig], vt1, m_Nc);
            wt2r = mult_udagv_r(&u[ic2 + ig], vt2, m_Nc);
            wt2i = mult_udagv_i(&u[ic2 + ig], vt2, m_Nc);

            w2[ic2 + id1 + iv]     += wt1r;
            w2[ic2 + 1 + id1 + iv] += wt1i;
            w2[ic2 + id2 + iv]     += wt2r;
            w2[ic2 + 1 + id2 + iv] += wt2i;
            w2[ic2 + id3 + iv]     += -wt1i;
            w2[ic2 + 1 + id3 + iv] += wt1r;
            w2[ic2 + id4 + iv]     += wt2i;
            w2[ic2 + 1 + id4 + iv] += -wt2r;
          }
        }
      }
    }
  }


//====================================================================
  void Fopr_Wilson_eo::mult_tp1_dirac_thread(
    int itask, double *vcp1, const double *v1, int ieo)
  {
    int Nvc2  = 2 * m_Nvc;
    int Nvcd  = m_Nvc * m_Nd;
    int Nvcd2 = Nvcd / 2;

    int id1 = 0;
    int id2 = m_Nvc;
    int id3 = m_Nvc * 2;
    int id4 = m_Nvc * 3;

    int idir = 3;

    int isite    = m_arg[itask].isite;
    int isite_cp = m_arg[itask].isite_cpt;

    //  double* w2 = &vcp1[Nvcd2*isite_cp];
    double *w2
      = (double *)m_bw_send[idir]->ptr(sizeof(double) * Nvcd2 * isite_cp);
    const double *w1 = &v1[Nvcd * isite];

    double bc2 = m_boundary2[idir];

    if (m_arg[itask].kt0 == 1) {
      int Nxy = m_Nx2 * m_Ny;
      int it  = 0;
      for (int iz = 0; iz < m_Mz; ++iz) {
        for (int ixy = 0; ixy < Nxy; ++ixy) {
          int is  = ixy + Nxy * (iz + m_Nz * it);
          int is2 = ixy + Nxy * iz;

          int in  = Nvcd * is;
          int ix1 = Nvc2 * is2;
          int ix2 = ix1 + m_Nvc;

          for (int ic = 0; ic < m_Nc; ++ic) {
            w2[2 * ic + ix1]     = 2.0 * bc2 * w1[2 * ic + id3 + in];
            w2[2 * ic + 1 + ix1] = 2.0 * bc2 * w1[2 * ic + 1 + id3 + in];
            w2[2 * ic + ix2]     = 2.0 * bc2 * w1[2 * ic + id4 + in];
            w2[2 * ic + 1 + ix2] = 2.0 * bc2 * w1[2 * ic + 1 + id4 + in];
          }
        }
      }
    }

    m_bw_send[idir]->start_thread(itask);
  }


//====================================================================
  void Fopr_Wilson_eo::mult_tp2_dirac_thread(
    int itask, double *v2, const double *vcp2, int ieo)
  {
    int Nvc2  = 2 * m_Nvc;
    int Nvcd  = m_Nvc * m_Nd;
    int Nvcd2 = Nvcd / 2;

    int id1 = 0;
    int id2 = m_Nvc;
    int id3 = m_Nvc * 2;
    int id4 = m_Nvc * 3;

    int idir = 3;

    double wt1r, wt1i, wt2r, wt2i;

    int isite    = m_arg[itask].isite;
    int isite_cp = m_arg[itask].isite_cpt;

    double *w2 = &v2[Nvcd * isite];
    //  double* w1 = &vcp2[Nvcd2*isite_cp];
    const double *w1
      = (double *)m_bw_recv[idir]->ptr(sizeof(double) * Nvcd2 * isite_cp);
    const double *u = m_U->ptr(m_Ndf * (isite + ieo * m_Nvol / 2 + idir * m_Nvol));

    m_bw_recv[idir]->wait_thread(itask);

    if (m_arg[itask].kt1 == 1) {
      int Nxy = m_Nx2 * m_Ny;
      int it  = m_Mt - 1;
      for (int iz = 0; iz < m_Mz; ++iz) {
        for (int ixy = 0; ixy < Nxy; ++ixy) {
          int is  = ixy + Nxy * (iz + m_Nz * it);
          int is2 = ixy + Nxy * iz;
          int iv  = Nvcd * is;
          int ig  = m_Ndf * is;
          int ix1 = Nvc2 * is2;
          int ix2 = ix1 + m_Nvc;

          for (int ic = 0; ic < m_Nc; ++ic) {
            int ic2 = ic * m_Nvc;

            wt1r = mult_uv_r(&u[ic2 + ig], &w1[ix1], m_Nc);
            wt1i = mult_uv_i(&u[ic2 + ig], &w1[ix1], m_Nc);
            wt2r = mult_uv_r(&u[ic2 + ig], &w1[ix2], m_Nc);
            wt2i = mult_uv_i(&u[ic2 + ig], &w1[ix2], m_Nc);

            w2[2 * ic + id3 + iv]     += wt1r;
            w2[2 * ic + 1 + id3 + iv] += wt1i;
            w2[2 * ic + id4 + iv]     += wt2r;
            w2[2 * ic + 1 + id4 + iv] += wt2i;
          }
        }
      }
    }
  }


//====================================================================
  void Fopr_Wilson_eo::mult_tpb_dirac_thread(
    int itask, double *v2, const double *v1, int ieo)
  {
    int Nvcd = m_Nvc * m_Nd;

    int id1 = 0;
    int id2 = m_Nvc;
    int id3 = m_Nvc * 2;
    int id4 = m_Nvc * 3;

    int idir = 3;

    double vt1[m_Nvc], vt2[m_Nvc];
    double wt1r, wt1i, wt2r, wt2i;

    int isite = m_arg[itask].isite;

    double       *w2 = &v2[Nvcd * isite];
    const double *w1 = &v1[Nvcd * isite];
    const double *u  = m_U->ptr(m_Ndf * (isite + ieo * m_Nvol / 2 + idir * m_Nvol));

    int kt1  = m_arg[itask].kt1;
    int Nxy  = m_Nx2 * m_Ny;
    int Nxyz = Nxy * m_Nz;

    for (int it = 0; it < m_Mt - kt1; ++it) {
      for (int iz = 0; iz < m_Mz; ++iz) {
        for (int ixy = 0; ixy < Nxy; ++ixy) {
          int is = ixy + Nxy * (iz + m_Nz * it);
          int iv = Nvcd * is;
          int in = Nvcd * (is + Nxyz);
          int ig = m_Ndf * is;

          for (int ic = 0; ic < m_Nc; ++ic) {
            vt1[2 * ic]     = 2.0 * w1[2 * ic + id3 + in];
            vt1[2 * ic + 1] = 2.0 * w1[2 * ic + 1 + id3 + in];
            vt2[2 * ic]     = 2.0 * w1[2 * ic + id4 + in];
            vt2[2 * ic + 1] = 2.0 * w1[2 * ic + 1 + id4 + in];
          }

          for (int ic = 0; ic < m_Nc; ++ic) {
            int ic2 = ic * m_Nvc;

            wt1r = mult_uv_r(&u[ic2 + ig], vt1, m_Nc);
            wt1i = mult_uv_i(&u[ic2 + ig], vt1, m_Nc);
            wt2r = mult_uv_r(&u[ic2 + ig], vt2, m_Nc);
            wt2i = mult_uv_i(&u[ic2 + ig], vt2, m_Nc);

            w2[2 * ic + id3 + iv]     += wt1r;
            w2[2 * ic + 1 + id3 + iv] += wt1i;
            w2[2 * ic + id4 + iv]     += wt2r;
            w2[2 * ic + 1 + id4 + iv] += wt2i;
          }
        }
      }
    }
  }


//====================================================================
  void Fopr_Wilson_eo::mult_tm1_dirac_thread(
    int itask, double *vcp1, const double *v1, int ieo)
  {
    int Nvc2  = 2 * m_Nvc;
    int Nvcd  = m_Nvc * m_Nd;
    int Nvcd2 = Nvcd / 2;

    int id1 = 0;
    int id2 = m_Nvc;
    int id3 = m_Nvc * 2;
    int id4 = m_Nvc * 3;

    int idir = 3;

    int isite    = m_arg[itask].isite;
    int isite_cp = m_arg[itask].isite_cpt;

    //  double* w2 = &vcp1[Nvcd2*isite_cp];
    double *w2
      = (double *)m_fw_send[idir]->ptr(sizeof(double) * Nvcd2 * isite_cp);
    const double *w1 = &v1[Nvcd * isite];
    const double *u  = m_U->ptr(m_Ndf * (isite + (1 - ieo) * m_Nvol / 2 + idir * m_Nvol));

    double vt1[m_Nvc], vt2[m_Nvc];

    if (m_arg[itask].kt1 == 1) {
      int Nxy = m_Nx2 * m_Ny;
      int it  = m_Mt - 1;
      for (int iz = 0; iz < m_Mz; ++iz) {
        for (int ixy = 0; ixy < Nxy; ++ixy) {
          int is  = ixy + Nxy * (iz + m_Nz * it);
          int is2 = ixy + Nxy * iz;
          int in  = Nvcd * is;
          int ig  = m_Ndf * is;
          int ix1 = Nvc2 * is2;
          int ix2 = ix1 + m_Nvc;

          for (int ic = 0; ic < m_Nc; ++ic) {
            vt1[2 * ic]     = 2.0 * w1[2 * ic + id1 + in];
            vt1[2 * ic + 1] = 2.0 * w1[2 * ic + 1 + id1 + in];
            vt2[2 * ic]     = 2.0 * w1[2 * ic + id2 + in];
            vt2[2 * ic + 1] = 2.0 * w1[2 * ic + 1 + id2 + in];
          }

          for (int ic = 0; ic < m_Nc; ++ic) {
            int icr = 2 * ic;
            w2[icr + ix1]     = mult_udagv_r(&u[icr + ig], vt1, m_Nc);
            w2[icr + 1 + ix1] = mult_udagv_i(&u[icr + ig], vt1, m_Nc);
            w2[icr + ix2]     = mult_udagv_r(&u[icr + ig], vt2, m_Nc);
            w2[icr + 1 + ix2] = mult_udagv_i(&u[icr + ig], vt2, m_Nc);
          }
        }
      }
    }

    m_fw_send[idir]->start_thread(itask);
  }


//====================================================================
  void Fopr_Wilson_eo::mult_tm2_dirac_thread(
    int itask, double *v2, const double *vcp2, int ieo)
  {
    int Nvc2  = 2 * m_Nvc;
    int Nvcd  = m_Nvc * m_Nd;
    int Nvcd2 = Nvcd / 2;

    int id1 = 0;
    int id2 = m_Nvc;
    int id3 = m_Nvc * 2;
    int id4 = m_Nvc * 3;

    int    idir = 3;
    double bc2  = m_boundary2[idir];

    double wt1r, wt1i, wt2r, wt2i;

    int isite    = m_arg[itask].isite;
    int isite_cp = m_arg[itask].isite_cpt;

    double *w2 = &v2[Nvcd * isite];
    //  double* w1 = &vcp2[Nvcd2*isite_cp];
    const double *w1
      = (double *)m_fw_recv[idir]->ptr(sizeof(double) * Nvcd2 * isite_cp);

    m_fw_recv[idir]->wait_thread(itask);

    if (m_arg[itask].kt0 == 1) {
      int Nxy = m_Nx2 * m_Ny;
      int it  = 0;
      for (int iz = 0; iz < m_Mz; ++iz) {
        for (int ixy = 0; ixy < Nxy; ++ixy) {
          int is  = ixy + Nxy * (iz + m_Nz * it);
          int is2 = ixy + Nxy * iz;
          int iv  = Nvcd * is;
          int ix1 = Nvc2 * is2;
          int ix2 = ix1 + m_Nvc;

          for (int ic = 0; ic < m_Nc; ++ic) {
            int icr = 2 * ic;
            int ici = 2 * ic + 1;
            w2[icr + id1 + iv] += bc2 * w1[icr + ix1];
            w2[ici + id1 + iv] += bc2 * w1[ici + ix1];
            w2[icr + id2 + iv] += bc2 * w1[icr + ix2];
            w2[ici + id2 + iv] += bc2 * w1[ici + ix2];
          }
        }
      }
    }
  }


//====================================================================
  void Fopr_Wilson_eo::mult_tmb_dirac_thread(
    int itask, double *v2, const double *v1, int ieo)
  {
    int Nvcd = m_Nvc * m_Nd;

    int id1 = 0;
    int id2 = m_Nvc;
    int id3 = m_Nvc * 2;
    int id4 = m_Nvc * 3;

    int idir = 3;

    double vt1[m_Nvc], vt2[m_Nvc];
    double wt1r, wt1i, wt2r, wt2i;

    int isite = m_arg[itask].isite;

    double       *w2 = &v2[Nvcd * isite];
    const double *w1 = &v1[Nvcd * isite];
    const double *u  = m_U->ptr(m_Ndf * (isite + (1 - ieo) * m_Nvol / 2 + idir * m_Nvol));

    int kt0  = m_arg[itask].kt0;
    int Nxy  = m_Nx2 * m_Ny;
    int Nxyz = Nxy * m_Nz;

    for (int it = kt0; it < m_Mt; ++it) {
      for (int iz = 0; iz < m_Mz; ++iz) {
        for (int ixy = 0; ixy < Nxy; ++ixy) {
          int is = ixy + Nxy * (iz + m_Nz * it);
          int iv = Nvcd * is;
          int in = Nvcd * (is - Nxyz);
          int ig = m_Ndf * (is - Nxyz);

          for (int ic = 0; ic < m_Nc; ++ic) {
            vt1[2 * ic]     = 2.0 * w1[2 * ic + id1 + in];
            vt1[2 * ic + 1] = 2.0 * w1[2 * ic + 1 + id1 + in];
            vt2[2 * ic]     = 2.0 * w1[2 * ic + id2 + in];
            vt2[2 * ic + 1] = 2.0 * w1[2 * ic + 1 + id2 + in];
          }

          for (int ic = 0; ic < m_Nc; ++ic) {
            int ic2 = 2 * ic;
            wt1r = mult_udagv_r(&u[ic2 + ig], vt1, m_Nc);
            wt1i = mult_udagv_i(&u[ic2 + ig], vt1, m_Nc);
            wt2r = mult_udagv_r(&u[ic2 + ig], vt2, m_Nc);
            wt2i = mult_udagv_i(&u[ic2 + ig], vt2, m_Nc);

            w2[ic2 + id1 + iv]     += wt1r;
            w2[ic2 + 1 + id1 + iv] += wt1i;
            w2[ic2 + id2 + iv]     += wt2r;
            w2[ic2 + 1 + id2 + iv] += wt2i;
          }
        }
      }
    }
  }


//====================================================================
  void Fopr_Wilson_eo::mult_tp1_chiral_thread(
    int itask, double *vcp1, const double *v1, int ieo)
  {
    int Nvc2  = 2 * m_Nvc;
    int Nvcd  = m_Nvc * m_Nd;
    int Nvcd2 = Nvcd / 2;

    int id1 = 0;
    int id2 = m_Nvc;
    int id3 = m_Nvc * 2;
    int id4 = m_Nvc * 3;

    int idir = 3;

    int isite    = m_arg[itask].isite;
    int isite_cp = m_arg[itask].isite_cpt;

    //  double* w2 = &vcp1[Nvcd2*isite_cp];
    double *w2
      = (double *)m_bw_send[idir]->ptr(sizeof(double) * Nvcd2 * isite_cp);
    const double *w1 = &v1[Nvcd * isite];

    double bc2 = m_boundary2[idir];

    if (m_arg[itask].kt0 == 1) {
      int Nxy = m_Nx2 * m_Ny;
      int it  = 0;
      for (int iz = 0; iz < m_Mz; ++iz) {
        for (int ixy = 0; ixy < Nxy; ++ixy) {
          int is  = ixy + Nxy * (iz + m_Nz * it);
          int is2 = ixy + Nxy * iz;

          int in  = Nvcd * is;
          int ix1 = Nvc2 * is2;
          int ix2 = ix1 + m_Nvc;

          for (int ic = 0; ic < m_Nc; ++ic) {
            w2[2 * ic + ix1]     = bc2 * (w1[2 * ic + id1 + in] + w1[2 * ic + id3 + in]);
            w2[2 * ic + 1 + ix1] = bc2 * (w1[2 * ic + 1 + id1 + in] + w1[2 * ic + 1 + id3 + in]);
            w2[2 * ic + ix2]     = bc2 * (w1[2 * ic + id2 + in] + w1[2 * ic + id4 + in]);
            w2[2 * ic + 1 + ix2] = bc2 * (w1[2 * ic + 1 + id2 + in] + w1[2 * ic + 1 + id4 + in]);
          }
        }
      }
    }

    m_bw_send[idir]->start_thread(itask);
  }


//====================================================================
  void Fopr_Wilson_eo::mult_tp2_chiral_thread(
    int itask, double *v2, const double *vcp2, int ieo)
  {
    int Nvc2  = 2 * m_Nvc;
    int Nvcd  = m_Nvc * m_Nd;
    int Nvcd2 = Nvcd / 2;

    int id1 = 0;
    int id2 = m_Nvc;
    int id3 = m_Nvc * 2;
    int id4 = m_Nvc * 3;

    int idir = 3;

    double wt1r, wt1i, wt2r, wt2i;

    int isite    = m_arg[itask].isite;
    int isite_cp = m_arg[itask].isite_cpt;

    double *w2 = &v2[Nvcd * isite];
    //  double* w1 = &vcp2[Nvcd2*isite_cp];
    const double *w1
      = (double *)m_bw_recv[idir]->ptr(sizeof(double) * Nvcd2 * isite_cp);
    const double *u = m_U->ptr(m_Ndf * (isite + ieo * m_Nvol / 2 + idir * m_Nvol));

    m_bw_recv[idir]->wait_thread(itask);

    if (m_arg[itask].kt1 == 1) {
      int Nxy = m_Nx2 * m_Ny;
      int it  = m_Mt - 1;
      for (int iz = 0; iz < m_Mz; ++iz) {
        for (int ixy = 0; ixy < Nxy; ++ixy) {
          int is  = ixy + Nxy * (iz + m_Nz * it);
          int is2 = ixy + Nxy * iz;
          int iv  = Nvcd * is;
          int ig  = m_Ndf * is;
          int ix1 = Nvc2 * is2;
          int ix2 = ix1 + m_Nvc;

          for (int ic = 0; ic < m_Nc; ++ic) {
            int ic2 = ic * m_Nvc;

            wt1r = mult_uv_r(&u[ic2 + ig], &w1[ix1], m_Nc);
            wt1i = mult_uv_i(&u[ic2 + ig], &w1[ix1], m_Nc);
            wt2r = mult_uv_r(&u[ic2 + ig], &w1[ix2], m_Nc);
            wt2i = mult_uv_i(&u[ic2 + ig], &w1[ix2], m_Nc);

            w2[2 * ic + id1 + iv]     += wt1r;
            w2[2 * ic + 1 + id1 + iv] += wt1i;
            w2[2 * ic + id2 + iv]     += wt2r;
            w2[2 * ic + 1 + id2 + iv] += wt2i;
            w2[2 * ic + id3 + iv]     += wt1r;
            w2[2 * ic + 1 + id3 + iv] += wt1i;
            w2[2 * ic + id4 + iv]     += wt2r;
            w2[2 * ic + 1 + id4 + iv] += wt2i;
          }
        }
      }
    }
  }


//====================================================================
  void Fopr_Wilson_eo::mult_tpb_chiral_thread(
    int itask, double *v2, const double *v1, int ieo)
  {
    int Nvcd = m_Nvc * m_Nd;

    int id1 = 0;
    int id2 = m_Nvc;
    int id3 = m_Nvc * 2;
    int id4 = m_Nvc * 3;

    int idir = 3;

    double vt1[m_Nvc], vt2[m_Nvc];
    double wt1r, wt1i, wt2r, wt2i;

    int isite = m_arg[itask].isite;

    double       *w2 = &v2[Nvcd * isite];
    const double *w1 = &v1[Nvcd * isite];
    const double *u  = m_U->ptr(m_Ndf * (isite + ieo * m_Nvol / 2 + idir * m_Nvol));

    int kt1  = m_arg[itask].kt1;
    int Nxy  = m_Nx2 * m_Ny;
    int Nxyz = Nxy * m_Nz;

    for (int it = 0; it < m_Mt - kt1; ++it) {
      for (int iz = 0; iz < m_Mz; ++iz) {
        for (int ixy = 0; ixy < Nxy; ++ixy) {
          int is = ixy + Nxy * (iz + m_Nz * it);
          int iv = Nvcd * is;
          int in = Nvcd * (is + Nxyz);
          int ig = m_Ndf * is;

          for (int ic = 0; ic < m_Nc; ++ic) {
            vt1[2 * ic]     = w1[2 * ic + id1 + in] + w1[2 * ic + id3 + in];
            vt1[2 * ic + 1] = w1[2 * ic + 1 + id1 + in] + w1[2 * ic + 1 + id3 + in];
            vt2[2 * ic]     = w1[2 * ic + id2 + in] + w1[2 * ic + id4 + in];
            vt2[2 * ic + 1] = w1[2 * ic + 1 + id2 + in] + w1[2 * ic + 1 + id4 + in];
          }

          for (int ic = 0; ic < m_Nc; ++ic) {
            int ic2 = ic * m_Nvc;

            wt1r = mult_uv_r(&u[ic2 + ig], vt1, m_Nc);
            wt1i = mult_uv_i(&u[ic2 + ig], vt1, m_Nc);
            wt2r = mult_uv_r(&u[ic2 + ig], vt2, m_Nc);
            wt2i = mult_uv_i(&u[ic2 + ig], vt2, m_Nc);

            w2[2 * ic + id1 + iv]     += wt1r;
            w2[2 * ic + 1 + id1 + iv] += wt1i;
            w2[2 * ic + id2 + iv]     += wt2r;
            w2[2 * ic + 1 + id2 + iv] += wt2i;
            w2[2 * ic + id3 + iv]     += wt1r;
            w2[2 * ic + 1 + id3 + iv] += wt1i;
            w2[2 * ic + id4 + iv]     += wt2r;
            w2[2 * ic + 1 + id4 + iv] += wt2i;
          }
        }
      }
    }
  }


//====================================================================
  void Fopr_Wilson_eo::mult_tm1_chiral_thread(
    int itask, double *vcp1, const double *v1, int ieo)
  {
    int Nvc2  = 2 * m_Nvc;
    int Nvcd  = m_Nvc * m_Nd;
    int Nvcd2 = Nvcd / 2;

    int id1 = 0;
    int id2 = m_Nvc;
    int id3 = m_Nvc * 2;
    int id4 = m_Nvc * 3;

    int idir = 3;

    int isite    = m_arg[itask].isite;
    int isite_cp = m_arg[itask].isite_cpt;

    //  double* w2 = &vcp1[Nvcd2*isite_cp];
    double *w2
      = (double *)m_fw_send[idir]->ptr(sizeof(double) * Nvcd2 * isite_cp);
    const double *w1 = &v1[Nvcd * isite];
    const double *u  = m_U->ptr(m_Ndf * (isite + (1 - ieo) * m_Nvol / 2 + idir * m_Nvol));

    double vt1[m_Nvc], vt2[m_Nvc];

    if (m_arg[itask].kt1 == 1) {
      int Nxy = m_Nx2 * m_Ny;
      int it  = m_Mt - 1;
      for (int iz = 0; iz < m_Mz; ++iz) {
        for (int ixy = 0; ixy < Nxy; ++ixy) {
          int is  = ixy + Nxy * (iz + m_Nz * it);
          int is2 = ixy + Nxy * iz;
          int in  = Nvcd * is;
          int ig  = m_Ndf * is;
          int ix1 = Nvc2 * is2;
          int ix2 = ix1 + m_Nvc;

          for (int ic = 0; ic < m_Nc; ++ic) {
            vt1[2 * ic]     = w1[2 * ic + id1 + in] - w1[2 * ic + id3 + in];
            vt1[2 * ic + 1] = w1[2 * ic + 1 + id1 + in] - w1[2 * ic + 1 + id3 + in];
            vt2[2 * ic]     = w1[2 * ic + id2 + in] - w1[2 * ic + id4 + in];
            vt2[2 * ic + 1] = w1[2 * ic + 1 + id2 + in] - w1[2 * ic + 1 + id4 + in];
          }

          for (int ic = 0; ic < m_Nc; ++ic) {
            int icr = 2 * ic;
            w2[icr + ix1]     = mult_udagv_r(&u[icr + ig], vt1, m_Nc);
            w2[icr + 1 + ix1] = mult_udagv_i(&u[icr + ig], vt1, m_Nc);
            w2[icr + ix2]     = mult_udagv_r(&u[icr + ig], vt2, m_Nc);
            w2[icr + 1 + ix2] = mult_udagv_i(&u[icr + ig], vt2, m_Nc);
          }
        }
      }
    }

    m_fw_send[idir]->start_thread(itask);
  }


//====================================================================
  void Fopr_Wilson_eo::mult_tm2_chiral_thread(
    int itask, double *v2, const double *vcp2, int ieo)
  {
    int Nvc2  = 2 * m_Nvc;
    int Nvcd  = m_Nvc * m_Nd;
    int Nvcd2 = Nvcd / 2;

    int id1 = 0;
    int id2 = m_Nvc;
    int id3 = m_Nvc * 2;
    int id4 = m_Nvc * 3;

    int    idir = 3;
    double bc2  = m_boundary2[idir];

    double wt1r, wt1i, wt2r, wt2i;

    int isite    = m_arg[itask].isite;
    int isite_cp = m_arg[itask].isite_cpt;

    double *w2 = &v2[Nvcd * isite];
    //  double* w1 = &vcp2[Nvcd2*isite_cp];
    const double *w1
      = (double *)m_fw_recv[idir]->ptr(sizeof(double) * Nvcd2 * isite_cp);

    m_fw_recv[idir]->wait_thread(itask);

    if (m_arg[itask].kt0 == 1) {
      int Nxy = m_Nx2 * m_Ny;
      int it  = 0;
      for (int iz = 0; iz < m_Mz; ++iz) {
        for (int ixy = 0; ixy < Nxy; ++ixy) {
          int is  = ixy + Nxy * (iz + m_Nz * it);
          int is2 = ixy + Nxy * iz;
          int iv  = Nvcd * is;
          int ix1 = Nvc2 * is2;
          int ix2 = ix1 + m_Nvc;

          for (int ic = 0; ic < m_Nc; ++ic) {
            int icr = 2 * ic;
            int ici = 2 * ic + 1;
            w2[icr + id1 + iv] += bc2 * w1[icr + ix1];
            w2[ici + id1 + iv] += bc2 * w1[ici + ix1];
            w2[icr + id2 + iv] += bc2 * w1[icr + ix2];
            w2[ici + id2 + iv] += bc2 * w1[ici + ix2];
            w2[icr + id3 + iv] -= bc2 * w1[icr + ix1];
            w2[ici + id3 + iv] -= bc2 * w1[ici + ix1];
            w2[icr + id4 + iv] -= bc2 * w1[icr + ix2];
            w2[ici + id4 + iv] -= bc2 * w1[ici + ix2];
          }
        }
      }
    }
  }


//====================================================================
  void Fopr_Wilson_eo::mult_tmb_chiral_thread(
    int itask, double *v2, const double *v1, int ieo)
  {
    int Nvcd = m_Nvc * m_Nd;

    int id1 = 0;
    int id2 = m_Nvc;
    int id3 = m_Nvc * 2;
    int id4 = m_Nvc * 3;

    int idir = 3;

    double vt1[m_Nvc], vt2[m_Nvc];
    double wt1r, wt1i, wt2r, wt2i;

    int isite = m_arg[itask].isite;

    double       *w2 = &v2[Nvcd * isite];
    const double *w1 = &v1[Nvcd * isite];
    const double *u  = m_U->ptr(m_Ndf * (isite + (1 - ieo) * m_Nvol / 2 + idir * m_Nvol));

    int kt0  = m_arg[itask].kt0;
    int Nxy  = m_Nx2 * m_Ny;
    int Nxyz = Nxy * m_Nz;

    for (int it = kt0; it < m_Mt; ++it) {
      for (int iz = 0; iz < m_Mz; ++iz) {
        for (int ixy = 0; ixy < Nxy; ++ixy) {
          int is = ixy + Nxy * (iz + m_Nz * it);
          int iv = Nvcd * is;
          int in = Nvcd * (is - Nxyz);
          int ig = m_Ndf * (is - Nxyz);

          for (int ic = 0; ic < m_Nc; ++ic) {
            vt1[2 * ic]     = w1[2 * ic + id1 + in] - w1[2 * ic + id3 + in];
            vt1[2 * ic + 1] = w1[2 * ic + 1 + id1 + in] - w1[2 * ic + 1 + id3 + in];
            vt2[2 * ic]     = w1[2 * ic + id2 + in] - w1[2 * ic + id4 + in];
            vt2[2 * ic + 1] = w1[2 * ic + 1 + id2 + in] - w1[2 * ic + 1 + id4 + in];
          }

          for (int ic = 0; ic < m_Nc; ++ic) {
            int ic2 = 2 * ic;
            wt1r = mult_udagv_r(&u[ic2 + ig], vt1, m_Nc);
            wt1i = mult_udagv_i(&u[ic2 + ig], vt1, m_Nc);
            wt2r = mult_udagv_r(&u[ic2 + ig], vt2, m_Nc);
            wt2i = mult_udagv_i(&u[ic2 + ig], vt2, m_Nc);

            w2[ic2 + id1 + iv]     += wt1r;
            w2[ic2 + 1 + id1 + iv] += wt1i;
            w2[ic2 + id2 + iv]     += wt2r;
            w2[ic2 + 1 + id2 + iv] += wt2i;
            w2[ic2 + id3 + iv]     -= wt1r;
            w2[ic2 + 1 + id3 + iv] -= wt1i;
            w2[ic2 + id4 + iv]     -= wt2r;
            w2[ic2 + 1 + id4 + iv] -= wt2i;
          }
        }
      }
    }
  }


//====================================================================
  void Fopr_Wilson_eo::gm5_dirac_thread(
    int itask, double *v2, const double *v1)
  {
    int Nvcd = m_Nvc * m_Nd;
    int Nxy  = m_Nx2 * m_Ny;

    int id1 = 0;
    int id2 = m_Nvc;
    int id3 = m_Nvc * 2;
    int id4 = m_Nvc * 3;

    int          isite = m_arg[itask].isite;
    double       *w2   = &v2[Nvcd * isite];
    const double *w1   = &v1[Nvcd * isite];

    for (int it = 0; it < m_Mt; ++it) {
      for (int iz = 0; iz < m_Mz; ++iz) {
        for (int ixy = 0; ixy < Nxy; ++ixy) {
          int iv = Nvcd * (ixy + Nxy * (iz + m_Nz * it));
          for (int ivc = 0; ivc < m_Nvc; ++ivc) {
            w2[ivc + id1 + iv] = w1[ivc + id3 + iv];
            w2[ivc + id2 + iv] = w1[ivc + id4 + iv];
            w2[ivc + id3 + iv] = w1[ivc + id1 + iv];
            w2[ivc + id4 + iv] = w1[ivc + id2 + iv];
          }
        }
      }
    }
  }


//====================================================================
  void Fopr_Wilson_eo::gm5_chiral_thread(
    int itask, double *v2, const double *v1)
  {
    int Nvcd = m_Nvc * m_Nd;
    int Nxy  = m_Nx2 * m_Ny;

    int id1 = 0;
    int id2 = m_Nvc;
    int id3 = m_Nvc * 2;
    int id4 = m_Nvc * 3;

    int          isite = m_arg[itask].isite;
    double       *w2   = &v2[Nvcd * isite];
    const double *w1   = &v1[Nvcd * isite];

    for (int it = 0; it < m_Mt; ++it) {
      for (int iz = 0; iz < m_Mz; ++iz) {
        for (int ixy = 0; ixy < Nxy; ++ixy) {
          int iv = Nvcd * (ixy + Nxy * (iz + m_Nz * it));
          for (int ivc = 0; ivc < m_Nvc; ++ivc) {
            w2[ivc + id1 + iv] = w1[ivc + id1 + iv];
            w2[ivc + id2 + iv] = w1[ivc + id2 + iv];
            w2[ivc + id3 + iv] = -w1[ivc + id3 + iv];
            w2[ivc + id4 + iv] = -w1[ivc + id4 + iv];
          }
        }
      }
    }
  }


//====================================================================
  void Fopr_Wilson_eo::gm5_dirac_thread(int itask,
                                        double *v1)
  {
    int Nvcd = m_Nvc * m_Nd;
    int Nxy  = m_Nx2 * m_Ny;

    int id1 = 0;
    int id2 = m_Nvc;
    int id3 = m_Nvc * 2;
    int id4 = m_Nvc * 3;

    int    isite = m_arg[itask].isite;
    double *w1   = &v1[Nvcd * isite];

    for (int it = 0; it < m_Mt; ++it) {
      for (int iz = 0; iz < m_Mz; ++iz) {
        for (int ixy = 0; ixy < Nxy; ++ixy) {
          int iv = Nvcd * (ixy + Nxy * (iz + m_Nz * it));
          for (int ivc = 0; ivc < m_Nvc; ++ivc) {
            double wt1 = w1[ivc + id1 + iv];
            double wt2 = w1[ivc + id2 + iv];
            w1[ivc + id1 + iv] = w1[ivc + id3 + iv];
            w1[ivc + id2 + iv] = w1[ivc + id4 + iv];
            w1[ivc + id3 + iv] = wt1;
            w1[ivc + id4 + iv] = wt2;
          }
        }
      }
    }
  }


//====================================================================
  void Fopr_Wilson_eo::gm5_chiral_thread(int itask,
                                         double *v1)
  {
    int Nvcd = m_Nvc * m_Nd;
    int Nxy  = m_Nx2 * m_Ny;

    int id1 = 0;
    int id2 = m_Nvc;
    int id3 = m_Nvc * 2;
    int id4 = m_Nvc * 3;

    int    isite = m_arg[itask].isite;
    double *w1   = &v1[Nvcd * isite];

    for (int it = 0; it < m_Mt; ++it) {
      for (int iz = 0; iz < m_Mz; ++iz) {
        for (int ixy = 0; ixy < Nxy; ++ixy) {
          int iv = Nvcd * (ixy + Nxy * (iz + m_Nz * it));
          for (int ivc = 0; ivc < m_Nvc; ++ivc) {
            w1[ivc + id3 + iv] = -w1[ivc + id3 + iv];
            w1[ivc + id4 + iv] = -w1[ivc + id4 + iv];
          }
        }
      }
    }
  }


//====================================================================
}
//============================================================END=====
#endif
