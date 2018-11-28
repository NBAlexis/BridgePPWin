#include "BridgeLib_Private.h"
#if USE_IMP

/*!
        @file    $Id:: fopr_Wilson_impl_thread.cpp #$

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 00:49:46 #$

        @version $LastChangedRevision: 1561 $
*/

#include <string>
using std::string;

#include "fopr_Wilson_impl.h"

#include "ResourceManager/threadManager_OpenMP.h"

namespace Imp {
#if defined USE_GROUP_SU3
#include "fopr_Wilson_impl_SU3.inc"
#elif defined USE_GROUP_SU2
#include "fopr_Wilson_impl_SU2.inc"
#elif defined USE_GROUP_SU_N
#include "fopr_Wilson_impl_SU_N.inc"
#endif

// const std::string Fopr_Wilson::class_name = "Fopr_Wilson";

//====================================================================
  void Fopr_Wilson::setup_thread()
  {
    m_Nthread = ThreadManager_OpenMP::get_num_threads_available();

    // The following setup corresponds to uniform division of volume.
    if (m_Nthread <= m_Nt) {
      m_Ntask_t = m_Nthread;
    } else if (m_Nthread <= m_Nz * m_Nt) {
      m_Ntask_t = m_Nt;
    } else {
      vout.crucial(m_vl, "Error at %s: Too large Nthread = %d\n", class_name.c_str(), m_Nthread);
      exit(EXIT_FAILURE);
    }

    m_Ntask_z = m_Nthread / m_Ntask_t;

    if (m_Ntask_z * m_Ntask_t != m_Nthread) {
      vout.crucial(m_vl, "Error at %s:  Nz = %d and Nt = %d do not match Nthread = %d\n",
                   class_name.c_str(), m_Nz, m_Nt, m_Nthread);
      exit(EXIT_FAILURE);
    }

    m_Ntask = m_Ntask_t * m_Ntask_z;
    m_Mz    = m_Nz / m_Ntask_z;
    m_Mt    = m_Nt / m_Ntask_t;

    if (m_Mz * m_Ntask_z != m_Nz) {
      vout.crucial(m_vl, "Error at %s:  Mz = %d and Ntask_z = %d do not match Nz = %d\n",
                   class_name.c_str(), m_Mz, m_Ntask_z, m_Nz);
      exit(EXIT_FAILURE);
    }

    if (m_Mt * m_Ntask_t != m_Nt) {
      vout.crucial(m_vl, "Error at %s:  Mt = %d and Ntask_t = %d do not match Nt = %d\n",
                   class_name.c_str(), m_Mt, m_Ntask_t, m_Nt);
      exit(EXIT_FAILURE);
    }

    // The following setup is not monotonic division, and requires
    // barrier at the beginning and end of mult (D and gamma5).
    //                                    [H.Matsufuru 22 Oct 2013]
    //  if(m_Nthread >= 64){
    //    m_Ntask_z = 8;
    //  }else if(m_Nthread >= 16){
    //    m_Ntask_z = 4;
    //  }else if(m_Nthread >= 4){
    //    m_Ntask_z = 2;
    //  }else{
    //    m_Ntask_z = 1;
    //  }
    //  m_Ntask_t = m_Nthread/m_Ntask_z;
    //  m_Ntask = m_Ntask_t * m_Ntask_z;
    //  m_Mz = m_Nz/m_Ntask_z;
    //  m_Mt = m_Nt/m_Ntask_t;

    vout.general(m_vl, "  Nthread = %d\n", m_Nthread);
    vout.general(m_vl, "  Ntask   = %d\n", m_Ntask);
    vout.general(m_vl, "  Ntask_z = %d  Ntask_t = %d\n", m_Ntask_z, m_Ntask_t);
    vout.general(m_vl, "  Mz      = %d  Mt      = %d\n", m_Mz, m_Mt);

    // setup of arguments
    int Nxy = m_Nx * m_Ny;
    m_arg.resize(m_Ntask);
    for (int ith_t = 0; ith_t < m_Ntask_t; ++ith_t) {
      for (int ith_z = 0; ith_z < m_Ntask_z; ++ith_z) {
        int itask = ith_z + m_Ntask_z * ith_t;

        m_arg[itask].isite = (ith_z * m_Mz + ith_t * (m_Nz * m_Mt)) * Nxy;

        m_arg[itask].kt0 = 0;
        m_arg[itask].kt1 = 0;
        m_arg[itask].kz0 = 0;
        m_arg[itask].kz1 = 0;
        if (ith_t == 0) m_arg[itask].kt0 = 1;
        if (ith_z == 0) m_arg[itask].kz0 = 1;
        if (ith_t == m_Ntask_t - 1) m_arg[itask].kt1 = 1;
        if (ith_z == m_Ntask_z - 1) m_arg[itask].kz1 = 1;

        m_arg[itask].isite_cpx = itask * m_Mz * m_Mt * m_Ny;
        m_arg[itask].isite_cpy = itask * m_Mz * m_Mt * m_Nx;
        m_arg[itask].isite_cpz = ith_t * m_Mt * Nxy;
        m_arg[itask].isite_cpt = ith_z * m_Mz * Nxy;
      }
    }
  }


//====================================================================
  void Fopr_Wilson::daypx_thread(int itask,
                                 double *v2, double fac, const double *v1)
  {
    int Nvcd = m_Nvc * m_Nd;
    int Nvxy = Nvcd * m_Nx * m_Ny;

    int          isite = m_arg[itask].isite;
    double       *w2   = &v2[Nvcd * isite];
    const double *w1   = &v1[Nvcd * isite];

    for (int it = 0; it < m_Mt; ++it) {
      for (int iz = 0; iz < m_Mz; ++iz) {
        for (int ivxy = 0; ivxy < Nvxy; ++ivxy) {
          int iv = ivxy + Nvxy * (iz + m_Nz * it);
          w2[iv] = fac * w2[iv] + w1[iv];
        }
      }
    }
  }


//====================================================================
  void Fopr_Wilson::clear_thread(int itask,
                                 double *v2)
  {
    int Nvcd = m_Nvc * m_Nd;
    int Nvxy = Nvcd * m_Nx * m_Ny;

    int    isite = m_arg[itask].isite;
    double *w2   = &v2[Nvcd * isite];

    for (int it = 0; it < m_Mt; ++it) {
      for (int iz = 0; iz < m_Mz; ++iz) {
        for (int ivxy = 0; ivxy < Nvxy; ++ivxy) {
          int iv = ivxy + Nvxy * (iz + m_Nz * it);
          w2[iv] = 0.0;
        }
      }
    }
  }


//====================================================================
  void Fopr_Wilson::mult_xp1_thread(int itask,
                                    double *vcp1, const double *v1)
  {
    int Nvc2  = 2 * m_Nvc;
    int Nvcd  = m_Nvc * m_Nd;
    int Nvcd2 = Nvcd / 2;

    int id1 = 0;
    int id2 = m_Nvc;
    int id3 = m_Nvc * 2;
    int id4 = m_Nvc * 3;

    int isite    = m_arg[itask].isite;
    int isite_cp = m_arg[itask].isite_cpx;

    double       *w2 = &vcp1[Nvcd2 * isite_cp];
    const double *w1 = &v1[Nvcd * isite];

    int    idir = 0;
    double bc2  = m_boundary2[idir];

    int ix = 0;

    for (int it = 0; it < m_Mt; ++it) {
      for (int iz = 0; iz < m_Mz; ++iz) {
        for (int iy = 0; iy < m_Ny; ++iy) {
          int is  = ix + m_Nx * (iy + m_Ny * (iz + m_Nz * it));
          int is2 = iy + m_Ny * (iz + m_Mz * it);
          int in  = Nvcd * is;
          int ix1 = Nvc2 * is2;
          int ix2 = ix1 + m_Nvc;

          for (int ic = 0; ic < m_Nc; ++ic) {
            w2[2 * ic + ix1]     = bc2 * (w1[2 * ic + id1 + in] - w1[2 * ic + 1 + id4 + in]);
            w2[2 * ic + 1 + ix1] = bc2 * (w1[2 * ic + 1 + id1 + in] + w1[2 * ic + id4 + in]);
            w2[2 * ic + ix2]     = bc2 * (w1[2 * ic + id2 + in] - w1[2 * ic + 1 + id3 + in]);
            w2[2 * ic + 1 + ix2] = bc2 * (w1[2 * ic + 1 + id2 + in] + w1[2 * ic + id3 + in]);
          }
        }
      }
    }
  }


//====================================================================
  void Fopr_Wilson::mult_xp2_thread(int itask,
                                    double *v2, const double *vcp2)
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

    double       *w2 = &v2[Nvcd * isite];
    const double *w1 = &vcp2[Nvcd2 * isite_cp];
    const double *u  = m_U->ptr(m_Ndf * (isite + idir * m_Nvol));

    int ix = m_Nx - 1;

    for (int it = 0; it < m_Mt; ++it) {
      for (int iz = 0; iz < m_Mz; ++iz) {
        for (int iy = 0; iy < m_Ny; ++iy) {
          int is  = ix + m_Nx * (iy + m_Ny * (iz + m_Nz * it));
          int is2 = iy + m_Ny * (iz + m_Mz * it);
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
            w2[2 * ic + id3 + iv]     += wt2i;
            w2[2 * ic + 1 + id3 + iv] += -wt2r;
            w2[2 * ic + id4 + iv]     += wt1i;
            w2[2 * ic + 1 + id4 + iv] += -wt1r;
          }
        }
      }
    }
  }


//====================================================================
  void Fopr_Wilson::mult_xpb_thread(int itask,
                                    double *v2, const double *v1)
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

    double       *w2 = &v2[Nvcd * isite];
    const double *w1 = &v1[Nvcd * isite];
    const double *u  = m_U->ptr(m_Ndf * (isite + idir * m_Nvol));

    for (int it = 0; it < m_Mt; ++it) {
      for (int iz = 0; iz < m_Mz; ++iz) {
        for (int iy = 0; iy < m_Ny; ++iy) {
          for (int ix = 0; ix < m_Nx - 1; ++ix) {
            int is = ix + m_Nx * (iy + m_Ny * (iz + m_Nz * it));
            int iv = Nvcd * is;
            int in = Nvcd * (is + 1);
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
  void Fopr_Wilson::mult_xm1_thread(int itask,
                                    double *vcp1, const double *v1)
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

    double       *w2 = &vcp1[Nvcd2 * isite_cp];
    const double *w1 = &v1[Nvcd * isite];
    const double *u  = m_U->ptr(m_Ndf * (isite + idir * m_Nvol));

    double vt1[m_Nvc], vt2[m_Nvc];

    int ix = m_Nx - 1;

    for (int it = 0; it < m_Mt; ++it) {
      for (int iz = 0; iz < m_Mz; ++iz) {
        for (int iy = 0; iy < m_Ny; ++iy) {
          int is  = ix + m_Nx * (iy + m_Ny * (iz + m_Nz * it));
          int is2 = iy + m_Ny * (iz + m_Mz * it);
          int in  = Nvcd * is;
          int ig  = m_Ndf * is;
          int ix1 = Nvc2 * is2;
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
        }
      }
    }
  }


//====================================================================
  void Fopr_Wilson::mult_xm2_thread(int itask,
                                    double *v2, const double *vcp2)
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

    // double wt1r, wt1i, wt2r, wt2i;

    int isite    = m_arg[itask].isite;
    int isite_cp = m_arg[itask].isite_cpx;

    double       *w2 = &v2[Nvcd * isite];
    const double *w1 = &vcp2[Nvcd2 * isite_cp];

    int ix = 0;

    for (int it = 0; it < m_Mt; ++it) {
      for (int iz = 0; iz < m_Mz; ++iz) {
        for (int iy = 0; iy < m_Ny; ++iy) {
          int is  = ix + m_Nx * (iy + m_Ny * (iz + m_Nz * it));
          int is2 = iy + m_Ny * (iz + m_Mz * it);
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
            w2[icr + id3 + iv] += -bc2 * w1[ici + ix2];
            w2[ici + id3 + iv] += +bc2 * w1[icr + ix2];
            w2[icr + id4 + iv] += -bc2 * w1[ici + ix1];
            w2[ici + id4 + iv] += +bc2 * w1[icr + ix1];
          }
        }
      }
    }
  }


//====================================================================
  void Fopr_Wilson::mult_xmb_thread(int itask,
                                    double *v2, const double *v1)
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

    double       *w2 = &v2[Nvcd * isite];
    const double *w1 = &v1[Nvcd * isite];
    const double *u  = m_U->ptr(m_Ndf * (isite + idir * m_Nvol));

    for (int it = 0; it < m_Mt; ++it) {
      for (int iz = 0; iz < m_Mz; ++iz) {
        for (int iy = 0; iy < m_Ny; ++iy) {
          for (int ix = 1; ix < m_Nx; ++ix) {
            int is = ix + m_Nx * (iy + m_Ny * (iz + m_Nz * it));
            int iv = Nvcd * is;
            int in = Nvcd * (is - 1);
            int ig = m_Ndf * (is - 1);

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
  void Fopr_Wilson::mult_yp1_thread(int itask,
                                    double *vcp1, const double *v1)
  {
    int Nvc2  = 2 * m_Nvc;
    int Nvcd  = m_Nvc * m_Nd;
    int Nvcd2 = Nvcd / 2;

    int id1 = 0;
    int id2 = m_Nvc;
    int id3 = m_Nvc * 2;
    int id4 = m_Nvc * 3;

    int isite    = m_arg[itask].isite;
    int isite_cp = m_arg[itask].isite_cpy;

    double       *w2 = &vcp1[Nvcd2 * isite_cp];
    const double *w1 = &v1[Nvcd * isite];

    int    idir = 1;
    double bc2  = m_boundary2[idir];

    int iy = 0;

    for (int it = 0; it < m_Mt; ++it) {
      for (int iz = 0; iz < m_Mz; ++iz) {
        for (int ix = 0; ix < m_Nx; ++ix) {
          int is  = ix + m_Nx * (iy + m_Ny * (iz + m_Nz * it));
          int is2 = ix + m_Nx * (iz + m_Mz * it);
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
  }


//====================================================================
  void Fopr_Wilson::mult_yp2_thread(int itask,
                                    double *v2, const double *vcp2)
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

    double       *w2 = &v2[Nvcd * isite];
    const double *w1 = &vcp2[Nvcd2 * isite_cp];
    const double *u  = m_U->ptr(m_Ndf * (isite + idir * m_Nvol));

    int iy = m_Ny - 1;

    for (int it = 0; it < m_Mt; ++it) {
      for (int iz = 0; iz < m_Mz; ++iz) {
        for (int ix = 0; ix < m_Nx; ++ix) {
          int is  = ix + m_Nx * (iy + m_Ny * (iz + m_Nz * it));
          int is2 = ix + m_Nx * (iz + m_Mz * it);
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
  void Fopr_Wilson::mult_ypb_thread(int itask,
                                    double *v2, const double *v1)
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
    const double *u  = m_U->ptr(m_Ndf * (isite + idir * m_Nvol));

    for (int it = 0; it < m_Mt; ++it) {
      for (int iz = 0; iz < m_Mz; ++iz) {
        for (int iy = 0; iy < m_Ny - 1; ++iy) {
          for (int ix = 0; ix < m_Nx; ++ix) {
            int is = ix + m_Nx * (iy + m_Ny * (iz + m_Nz * it));
            int iv = Nvcd * is;
            int in = Nvcd * (is + m_Nx);
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
  void Fopr_Wilson::mult_ym1_thread(int itask,
                                    double *vcp1, const double *v1)
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

    double       *w2 = &vcp1[Nvcd2 * isite_cp];
    const double *w1 = &v1[Nvcd * isite];
    const double *u  = m_U->ptr(m_Ndf * (isite + idir * m_Nvol));

    double vt1[m_Nvc], vt2[m_Nvc];

    int iy = m_Ny - 1;

    for (int it = 0; it < m_Mt; ++it) {
      for (int iz = 0; iz < m_Mz; ++iz) {
        for (int ix = 0; ix < m_Nx; ++ix) {
          int is  = ix + m_Nx * (iy + m_Ny * (iz + m_Nz * it));
          int is2 = ix + m_Nx * (iz + m_Mz * it);
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
  }


//====================================================================
  void Fopr_Wilson::mult_ym2_thread(int itask,
                                    double *v2, const double *vcp2)
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

    // double wt1r, wt1i, wt2r, wt2i;

    int isite    = m_arg[itask].isite;
    int isite_cp = m_arg[itask].isite_cpy;

    double       *w2 = &v2[Nvcd * isite];
    const double *w1 = &vcp2[Nvcd2 * isite_cp];

    int iy = 0;

    for (int it = 0; it < m_Mt; ++it) {
      for (int iz = 0; iz < m_Mz; ++iz) {
        for (int ix = 0; ix < m_Nx; ++ix) {
          int is  = ix + m_Nx * (iy + m_Ny * (iz + m_Nz * it));
          int is2 = ix + m_Nx * (iz + m_Mz * it);
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
  void Fopr_Wilson::mult_ymb_thread(int itask,
                                    double *v2, const double *v1)
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
    const double *u  = m_U->ptr(m_Ndf * (isite + idir * m_Nvol));

    for (int it = 0; it < m_Mt; ++it) {
      for (int iz = 0; iz < m_Mz; ++iz) {
        for (int iy = 1; iy < m_Ny; ++iy) {
          for (int ix = 0; ix < m_Nx; ++ix) {
            int is = ix + m_Nx * (iy + m_Ny * (iz + m_Nz * it));
            int iv = Nvcd * is;
            int in = Nvcd * (is - m_Nx);
            int ig = m_Ndf * (is - m_Nx);

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
  void Fopr_Wilson::mult_zp1_thread(int itask,
                                    double *vcp1, const double *v1)
  {
    int Nvc2  = 2 * m_Nvc;
    int Nvcd  = m_Nvc * m_Nd;
    int Nvcd2 = Nvcd / 2;

    int id1 = 0;
    int id2 = m_Nvc;
    int id3 = m_Nvc * 2;
    int id4 = m_Nvc * 3;

    int isite    = m_arg[itask].isite;
    int isite_cp = m_arg[itask].isite_cpz;

    double       *w2 = &vcp1[Nvcd2 * isite_cp];
    const double *w1 = &v1[Nvcd * isite];

    int    idir = 2;
    double bc2  = m_boundary2[idir];

    if (m_arg[itask].kz0 == 1) {
      int Nxy = m_Nx * m_Ny;
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
  }


//====================================================================
  void Fopr_Wilson::mult_zp2_thread(int itask,
                                    double *v2, const double *vcp2)
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

    double       *w2 = &v2[Nvcd * isite];
    const double *w1 = &vcp2[Nvcd2 * isite_cp];
    const double *u  = m_U->ptr(m_Ndf * (isite + idir * m_Nvol));

    if (m_arg[itask].kz1 == 1) {
      int Nxy = m_Nx * m_Ny;
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
  void Fopr_Wilson::mult_zpb_thread(int itask,
                                    double *v2, const double *v1)
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
    const double *u  = m_U->ptr(m_Ndf * (isite + idir * m_Nvol));

    int kz1 = m_arg[itask].kz1;
    int Nxy = m_Nx * m_Ny;

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
  void Fopr_Wilson::mult_zm1_thread(int itask,
                                    double *vcp1, const double *v1)
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

    double       *w2 = &vcp1[Nvcd2 * isite_cp];
    const double *w1 = &v1[Nvcd * isite];
    const double *u  = m_U->ptr(m_Ndf * (isite + idir * m_Nvol));

    double vt1[m_Nvc], vt2[m_Nvc];

    if (m_arg[itask].kz1 == 1) {
      int Nxy = m_Nx * m_Ny;
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
  }


//====================================================================
  void Fopr_Wilson::mult_zm2_thread(int itask,
                                    double *v2, const double *vcp2)
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

    // double wt1r, wt1i, wt2r, wt2i;

    int isite    = m_arg[itask].isite;
    int isite_cp = m_arg[itask].isite_cpz;

    double       *w2 = &v2[Nvcd * isite];
    const double *w1 = &vcp2[Nvcd2 * isite_cp];

    if (m_arg[itask].kz0 == 1) {
      int Nxy = m_Nx * m_Ny;

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
  void Fopr_Wilson::mult_zmb_thread(int itask,
                                    double *v2, const double *v1)
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
    const double *u  = m_U->ptr(m_Ndf * (isite + idir * m_Nvol));

    int kz0 = m_arg[itask].kz0;
    int Nxy = m_Nx * m_Ny;

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
  void Fopr_Wilson::mult_tp1_dirac_thread(int itask,
                                          double *vcp1, const double *v1)
  {
    int Nvc2  = 2 * m_Nvc;
    int Nvcd  = m_Nvc * m_Nd;
    int Nvcd2 = Nvcd / 2;

    // int id1 = 0;
    // int id2 = m_Nvc;
    int id3 = m_Nvc * 2;
    int id4 = m_Nvc * 3;

    int isite    = m_arg[itask].isite;
    int isite_cp = m_arg[itask].isite_cpt;

    double       *w2 = &vcp1[Nvcd2 * isite_cp];
    const double *w1 = &v1[Nvcd * isite];

    int    idir = 3;
    double bc2  = m_boundary2[idir];

    if (m_arg[itask].kt0 == 1) {
      int Nxy = m_Nx * m_Ny;
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
  }


//====================================================================
  void Fopr_Wilson::mult_tp2_dirac_thread(int itask,
                                          double *v2, const double *vcp2)
  {
    int Nvc2  = 2 * m_Nvc;
    int Nvcd  = m_Nvc * m_Nd;
    int Nvcd2 = Nvcd / 2;

    // int id1 = 0;
    // int id2 = m_Nvc;
    int id3 = m_Nvc * 2;
    int id4 = m_Nvc * 3;

    int idir = 3;

    double wt1r, wt1i, wt2r, wt2i;

    int isite    = m_arg[itask].isite;
    int isite_cp = m_arg[itask].isite_cpt;

    double       *w2 = &v2[Nvcd * isite];
    const double *w1 = &vcp2[Nvcd2 * isite_cp];
    const double *u  = m_U->ptr(m_Ndf * (isite + idir * m_Nvol));

    if (m_arg[itask].kt1 == 1) {
      int Nxy = m_Nx * m_Ny;
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
  void Fopr_Wilson::mult_tpb_dirac_thread(int itask,
                                          double *v2, const double *v1)
  {
    int Nvcd = m_Nvc * m_Nd;

    // int id1 = 0;
    // int id2 = m_Nvc;
    int id3 = m_Nvc * 2;
    int id4 = m_Nvc * 3;

    int idir = 3;

    double vt1[m_Nvc], vt2[m_Nvc];
    double wt1r, wt1i, wt2r, wt2i;

    int isite = m_arg[itask].isite;

    double       *w2 = &v2[Nvcd * isite];
    const double *w1 = &v1[Nvcd * isite];
    const double *u  = m_U->ptr(m_Ndf * (isite + idir * m_Nvol));

    int kt1  = m_arg[itask].kt1;
    int Nxy  = m_Nx * m_Ny;
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
  void Fopr_Wilson::mult_tm1_dirac_thread(int itask,
                                          double *vcp1, const double *v1)
  {
    int Nvc2  = 2 * m_Nvc;
    int Nvcd  = m_Nvc * m_Nd;
    int Nvcd2 = Nvcd / 2;

    int id1 = 0;
    int id2 = m_Nvc;
    // int id3 = m_Nvc * 2;
    // int id4 = m_Nvc * 3;

    int idir = 3;

    int isite    = m_arg[itask].isite;
    int isite_cp = m_arg[itask].isite_cpt;

    double       *w2 = &vcp1[Nvcd2 * isite_cp];
    const double *w1 = &v1[Nvcd * isite];
    const double *u  = m_U->ptr(m_Ndf * (isite + idir * m_Nvol));

    double vt1[m_Nvc], vt2[m_Nvc];

    if (m_arg[itask].kt1 == 1) {
      int Nxy = m_Nx * m_Ny;
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
  }


//====================================================================
  void Fopr_Wilson::mult_tm2_dirac_thread(int itask,
                                          double *v2, const double *vcp2)
  {
    int Nvc2  = 2 * m_Nvc;
    int Nvcd  = m_Nvc * m_Nd;
    int Nvcd2 = Nvcd / 2;

    int id1 = 0;
    int id2 = m_Nvc;
    // int id3 = m_Nvc * 2;
    // int id4 = m_Nvc * 3;

    int    idir = 3;
    double bc2  = m_boundary2[idir];

    // double wt1r, wt1i, wt2r, wt2i;

    int isite    = m_arg[itask].isite;
    int isite_cp = m_arg[itask].isite_cpt;

    double       *w2 = &v2[Nvcd * isite];
    const double *w1 = &vcp2[Nvcd2 * isite_cp];

    if (m_arg[itask].kt0 == 1) {
      int Nxy = m_Nx * m_Ny;
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
  void Fopr_Wilson::mult_tmb_dirac_thread(int itask,
                                          double *v2, const double *v1)
  {
    int Nvcd = m_Nvc * m_Nd;

    int id1 = 0;
    int id2 = m_Nvc;
    // int id3 = m_Nvc * 2;
    // int id4 = m_Nvc * 3;

    int idir = 3;

    double vt1[m_Nvc], vt2[m_Nvc];
    double wt1r, wt1i, wt2r, wt2i;

    int isite = m_arg[itask].isite;

    double       *w2 = &v2[Nvcd * isite];
    const double *w1 = &v1[Nvcd * isite];
    const double *u  = m_U->ptr(m_Ndf * (isite + idir * m_Nvol));

    int kt0  = m_arg[itask].kt0;
    int Nxy  = m_Nx * m_Ny;
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
  void Fopr_Wilson::mult_tp1_chiral_thread(int itask,
                                           double *vcp1, const double *v1)
  {
    int Nvc2  = 2 * m_Nvc;
    int Nvcd  = m_Nvc * m_Nd;
    int Nvcd2 = Nvcd / 2;

    int id1 = 0;
    int id2 = m_Nvc;
    int id3 = m_Nvc * 2;
    int id4 = m_Nvc * 3;

    int isite    = m_arg[itask].isite;
    int isite_cp = m_arg[itask].isite_cpt;

    double       *w2 = &vcp1[Nvcd2 * isite_cp];
    const double *w1 = &v1[Nvcd * isite];

    int    idir = 3;
    double bc2  = m_boundary2[idir];

    if (m_arg[itask].kt0 == 1) {
      int Nxy = m_Nx * m_Ny;
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
  }


//====================================================================
  void Fopr_Wilson::mult_tp2_chiral_thread(int itask,
                                           double *v2, const double *vcp2)
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

    double       *w2 = &v2[Nvcd * isite];
    const double *w1 = &vcp2[Nvcd2 * isite_cp];
    const double *u  = m_U->ptr(m_Ndf * (isite + idir * m_Nvol));

    if (m_arg[itask].kt1 == 1) {
      int Nxy = m_Nx * m_Ny;
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
  void Fopr_Wilson::mult_tpb_chiral_thread(int itask,
                                           double *v2, const double *v1)
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
    const double *u  = m_U->ptr(m_Ndf * (isite + idir * m_Nvol));

    int kt1  = m_arg[itask].kt1;
    int Nxy  = m_Nx * m_Ny;
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
  void Fopr_Wilson::mult_tm1_chiral_thread(int itask,
                                           double *vcp1, const double *v1)
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

    double       *w2 = &vcp1[Nvcd2 * isite_cp];
    const double *w1 = &v1[Nvcd * isite];
    const double *u  = m_U->ptr(m_Ndf * (isite + idir * m_Nvol));

    double vt1[m_Nvc], vt2[m_Nvc];

    if (m_arg[itask].kt1 == 1) {
      int Nxy = m_Nx * m_Ny;
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
  }


//====================================================================
  void Fopr_Wilson::mult_tm2_chiral_thread(int itask,
                                           double *v2, const double *vcp2)
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

    // double wt1r, wt1i, wt2r, wt2i;

    int isite    = m_arg[itask].isite;
    int isite_cp = m_arg[itask].isite_cpt;

    double       *w2 = &v2[Nvcd * isite];
    const double *w1 = &vcp2[Nvcd2 * isite_cp];

    if (m_arg[itask].kt0 == 1) {
      int Nxy = m_Nx * m_Ny;
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
  void Fopr_Wilson::mult_tmb_chiral_thread(int itask,
                                           double *v2, const double *v1)
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
    const double *u  = m_U->ptr(m_Ndf * (isite + idir * m_Nvol));

    int kt0  = m_arg[itask].kt0;
    int Nxy  = m_Nx * m_Ny;
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
  void Fopr_Wilson::gm5_dirac_thread(int itask,
                                     double *v2, const double *v1)
  {
    int Nvcd = m_Nvc * m_Nd;
    int Nxy  = m_Nx * m_Ny;

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
  void Fopr_Wilson::gm5_chiral_thread(int itask,
                                      double *v2, const double *v1)
  {
    int Nvcd = m_Nvc * m_Nd;
    int Nxy  = m_Nx * m_Ny;

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
}
//============================================================END=====

#endif
