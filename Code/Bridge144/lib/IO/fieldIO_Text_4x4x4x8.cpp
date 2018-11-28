#include "BridgeLib_Private.h"

/*!
        @file    $Id: fieldIO_Text_4x4x4x8.cpp #$

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate: 2014-10-09 17:20:36 #$

        @version $LastChangedRevision: 1561 $
*/

#include "fieldIO_Text_4x4x4x8.h"

#include <fstream>

#include "bridgeIO.h"
using Bridge::vout;

const std::string FieldIO_Text_4x4x4x8::class_name = "FieldIO_Text_4x4x4x8";

//====================================================================
void FieldIO_Text_4x4x4x8::read_file(Field *v, string filename)
{
  int nin_field = v->nin();
  int nex_field = v->nex();

  int   Nc    = CommonParameters::Nc();
  int   Ndim  = CommonParameters::Ndim();
  int   NinG  = 2 * Nc * Nc;
  int   NxF   = 4;
  int   NyF   = 4;
  int   NzF   = 4;
  int   NtF   = 8;
  int   LvolF = NxF * NyF * NzF * NtF;
  Field v_in(NinG, LvolF, Ndim);

  int nin_file = NinG;
  int nex_file = Ndim;

  assert(nin_file == nin_field);
  assert(nex_file == nex_field);

  int Lvol = CommonParameters::Lvol();

  vout.detailed(m_vl, "%s: file format: nin=%d, nex=%d, Lvol=%d\n", __func__, nin_file, nex_file, Lvol);
  vout.detailed(m_vl, "%s: field format: nin=%d, nex=%d, Lvol=%d\n", __func__, nin_field, nex_field, v->nvol());

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
  int count = NinG * LvolF * Ndim;
  Communicator::broadcast(count, (double *)vtmp.ptr(0), 0);

  Index_lex idx_f;
  Index_lex indexF(NxF, NyF, NzF, NtF);

  int Nx = CommonParameters::Nx();
  int Ny = CommonParameters::Ny();
  int Nz = CommonParameters::Nz();
  int Nt = CommonParameters::Nt();

  int ipex = Communicator::ipe(0);
  int ipey = Communicator::ipe(1);
  int ipez = Communicator::ipe(2);
  int ipet = Communicator::ipe(3);

  for (int j = 0; j < Ndim; ++j) {
    for (int t = 0; t < Nt; ++t) {
      int t2 = (t + ipet * Nt) % NtF;
      for (int z = 0; z < Nz; ++z) {
        int z2 = (z + ipez * Nz) % NzF;
        for (int y = 0; y < Ny; ++y) {
          int y2 = (y + ipey * Ny) % NyF;
          for (int x = 0; x < Nx; ++x) {
            int x2    = (x + ipex * Nx) % NxF;
            int lsite = idx_f.site(x, y, z, t);
            int gsite = indexF.site(x2, y2, z2, t2);

            for (int i = 0; i < NinG; ++i) {
              v->set(i, lsite, j, vtmp.cmp(i, gsite, j));
            }
          }
        }
      }
    }
  }
}


//====================================================================
void FieldIO_Text_4x4x4x8::write_file(Field *, std::string )
{
  vout.crucial(m_vl, "Warning at %s: no write method is defined.\n",
               class_name.c_str());
}


//====================================================================
//============================================================END=====
