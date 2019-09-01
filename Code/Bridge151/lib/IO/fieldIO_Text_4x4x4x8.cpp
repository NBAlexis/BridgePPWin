/*!
        @file    fieldIO_Text_4x4x4x8.cpp

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate: 2014-10-09 17:20:36 #$

        @version $LastChangedRevision: 1929 $
*/
#include "BridgeLib_Private.h"
#include "fieldIO_Text_4x4x4x8.h"

#include <fstream>

const std::string FieldIO_Text_4x4x4x8::class_name = "FieldIO_Text_4x4x4x8";

//====================================================================
void FieldIO_Text_4x4x4x8::read_file(Field *v, const std::string filename)
{
  const int nin_field = v->nin();
  const int nex_field = v->nex();

  const int Nc   = CommonParameters::Nc();
  const int Ndim = CommonParameters::Ndim();
  const int NinG = 2 * Nc * Nc;

  const long_t Lvol = CommonParameters::Lvol();

  const int NxF = 4;
  const int NyF = 4;
  const int NzF = 4;
  const int NtF = 8;

  const int   LvolF = NxF * NyF * NzF * NtF;
  const Field v_in(NinG, LvolF, Ndim);

  const int nin_file = NinG;
  const int nex_file = Ndim;

  assert(nin_file == nin_field);
  assert(nex_file == nex_field);


  vout.detailed(m_vl, "%s: file format: nin=%d, nex=%d, Lvol=%ld\n", __func__, nin_file, nex_file, Lvol);
  vout.detailed(m_vl, "%s: field format: nin=%d, nex=%d, Lvol=%ld\n", __func__, nin_field, nex_field, v->nvol());

  // temporary field holding the whole space-time data.
  //  Field vtmp;
  Field vtmp(nin_field, LvolF, nex_field);

  if (Communicator::is_primary()) {
    vout.detailed(m_vl, "reading field data from %s\n", filename.c_str());

    //    vtmp.reset(nin_field, Lvol, nex_field);

    std::fstream config(filename.c_str(), std::ios::in);
    if (!config.is_open()) {
      vout.crucial(m_vl, "Error at %s: file open error: %s may not exist.\n", __func__, filename.c_str());
      exit(EXIT_FAILURE);
    }

    // int    s, t;
    double val;
    for (int isite = 0; isite < LvolF; ++isite) {
      for (int j = 0; j < nex_file; ++j) {
        for (int i = 0; i < nin_file; ++i) {
          config >> val;

          if (!config.good()) {
            if (config.eof()) {
              // file short.
              vout.crucial(m_vl, "Error at %s: file size too small.\n", __func__);
              exit(EXIT_FAILURE);
            }
            if (config.fail()) {
              // invalid data.
              vout.crucial(m_vl, "Error at %s: invalid data.\n", __func__);
              exit(EXIT_FAILURE);
            }
            // other error.
            vout.crucial(m_vl, "Error at %s: io error.\n", __func__);
            exit(EXIT_FAILURE);
          }

          //  m_format->file_to_field(s, t, i, j);
          vtmp.set(i, isite, j, val);
        }
      }
    }

    int dummy;
    config >> dummy;

    if (!config.eof()) {
      // file size mismatch
      vout.crucial(m_vl, "Warning at %s: file size larger than expected.\n", __func__);
      // warning only
    }

    config.close();
  }
  Communicator::sync();

  vout.detailed(m_vl, "read successful\n");

  // copy read configuration to all the nodes.
  const int count = NinG * LvolF * Ndim;
  Communicator::broadcast(count, (double *)vtmp.ptr(0), 0);

  Index_lex idx_local;
  Index_lex indexF(NxF, NyF, NzF, NtF);

  const int Nx = CommonParameters::Nx();
  const int Ny = CommonParameters::Ny();
  const int Nz = CommonParameters::Nz();
  const int Nt = CommonParameters::Nt();

  const int ipe_x = Communicator::ipe(0);
  const int ipe_y = Communicator::ipe(1);
  const int ipe_z = Communicator::ipe(2);
  const int ipe_t = Communicator::ipe(3);

  for (int j = 0; j < Ndim; ++j) {
    for (int t = 0; t < Nt; ++t) {
      int t_global = (t + ipe_t * Nt) % NtF;

      for (int z = 0; z < Nz; ++z) {
        int z_global = (z + ipe_z * Nz) % NzF;

        for (int y = 0; y < Ny; ++y) {
          int y_global = (y + ipe_y * Ny) % NyF;

          for (int x = 0; x < Nx; ++x) {
            int x_global = (x + ipe_x * Nx) % NxF;

            int local_site  = idx_local.site(x, y, z, t);
            int global_site = indexF.site(x_global, y_global, z_global, t_global);

            for (int i = 0; i < NinG; ++i) {
              v->set(i, local_site, j, vtmp.cmp(i, global_site, j));
            }
          }
        }
      }
    }
  }
}


//====================================================================
void FieldIO_Text_4x4x4x8::write_file(Field *v, const std::string filename)
{
  vout.crucial(m_vl, "Warning at %s: no write method is defined.\n",
               class_name.c_str());
}


//====================================================================
//============================================================END=====
