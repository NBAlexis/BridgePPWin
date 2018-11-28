/*!
        @file    $Id: fieldIO_Text.cpp #$

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate: 2014-10-09 17:20:36 #$

        @version $LastChangedRevision: 1561 $
*/

#include "fieldIO_Text.h"

#include <fstream>

#include "bridgeIO.h"
using Bridge::vout;

const std::string FieldIO_Text::class_name = "FieldIO_Text";

//====================================================================
void FieldIO_Text::read_file(Field *v, string filename)
{
  int nin_field = v->nin();
  int nex_field = v->nex();

  int nin_file = m_format->nin();
  int nex_file = m_format->nex();

  if ((nin_file == 0) || (nex_file == 0)) {
    nin_file = nin_field;
    nex_file = nex_field;
  }

  int Lvol = CommonParameters::Lvol();

  vout.detailed(m_vl, "%s: file format: nin=%d, nex=%d, Lvol=%d\n", __func__, nin_file, nex_file, Lvol);
  vout.detailed(m_vl, "%s: field format: nin=%d, nex=%d, Lvol=%d\n", __func__, nin_field, nex_field, v->nvol());

  // temporary field holding the whole space-time data.
  Field vtmp;

  if (Communicator::is_primary()) {
    vout.detailed(m_vl, "reading field data from %s\n", filename.c_str());

    vtmp.reset(nin_field, Lvol, nex_field);

    std::fstream config(filename.c_str(), std::ios::in);
    if (!config.is_open()) {
      vout.crucial(m_vl, "Error at %s: file open error: %s may not exist.\n", __func__, filename.c_str());
      exit(EXIT_FAILURE);
    }

    int    s, t;
    double val;
    for (int j = 0; j < nex_file; ++j) {
      for (int isite = 0; isite < Lvol; ++isite) {
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

          m_format->file_to_field(s, t, i, j);
          vtmp.set(s, isite, t, val);
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

  FieldIO::deliver(v, &vtmp);

  vout.detailed(m_vl, "read successful\n");
}


//====================================================================
void FieldIO_Text::write_file(Field *v, std::string filename)
{
  if (v == 0) {
    vout.crucial(m_vl, "%s: field empty.\n", __func__);
    return;
  }

  int nin_field = v->nin();
  int nex_field = v->nex();

  int nin_file = m_format->nin();
  int nex_file = m_format->nex();

  if ((nin_file == 0) || (nex_file == 0)) {
    nin_file = nin_field;
    nex_file = nex_field;
  }

  int Lvol = CommonParameters::Lvol();

  vout.detailed(m_vl, "%s: file format: nin=%d, nex=%d, Lvol=%d\n", __func__, nin_file, nex_file, Lvol);
  vout.detailed(m_vl, "%s: field format: nin=%d, nex=%d, Lvol=%d\n", __func__, nin_field, nex_field, v->nvol());

  // temporary field holding the whole space-time data.
  Field vtmp;
  if (Communicator::is_primary()) {
    vtmp.reset(nin_field, Lvol, nex_field);
  }

  FieldIO::gather(&vtmp, v);

  if (Communicator::is_primary()) {
    vout.detailed(m_vl, "writing field data to %s\n", filename.c_str());

    std::fstream config(filename.c_str(), std::ios::out);
    if (!config.is_open()) {
      vout.crucial(m_vl, "Error at %s: file open error: %s\n", __func__, filename.c_str());
      exit(EXIT_FAILURE);
    }

    config.setf(std::ios_base::scientific, std::ios_base::floatfield);
    config.precision(14);

    int    s, t;
    double val;
    for (int j = 0; j < nex_file; ++j) {
      for (int isite = 0; isite < Lvol; ++isite) {
        for (int i = 0; i < nin_file; ++i) {
          m_format->file_to_field(s, t, i, j);

          val = vtmp.cmp(s, isite, t);
          config << val << std::endl;

          if (config.fail()) {
            vout.crucial(m_vl, "Error at %s: write failed.\n", __func__);
            exit(EXIT_FAILURE);
          }
        }
      }
    }

    config.close();
  }

  vout.detailed(m_vl, "write successful\n");
}


//====================================================================
//============================================================END=====
