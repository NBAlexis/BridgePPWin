/*!
        @file    $Id: fieldIO_Binary_Distributed.cpp #$

        @brief

        @author  Tatsumi Aoyama (aoym)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate: 2013-01-22 13:51:53 +0900 #$

        @version $LastChangedRevision: 1561 $
*/

#include "fieldIO_Binary_Distributed.h"
#include "Tools/file_utils.h"

const std::string FieldIO_Binary_Distributed::class_name = "FieldIO_Binary_Distributed";

//====================================================================
void FieldIO_Binary_Distributed::read_file(Field *u, string filename_base)
{
  static const char _function_name[] = "FieldIO_Binary_Distributed::read_file";
  string            filename
    = FileUtils::generate_filename("%s.x%dy%dz%dt%d",
                                   filename_base.c_str(),
                                   Communicator::ipe(0),
                                   Communicator::ipe(1),
                                   Communicator::ipe(2),
                                   Communicator::ipe(3));

  int nin_field = u->nin();
  int nex_field = u->nex();

  int nin_file = m_format->nin();
  int nex_file = m_format->nex();

  if ((nin_file == 0) || (nex_file == 0)) {
    nin_file = nin_field;
    nex_file = nex_field;
  }

  int nvol = CommonParameters::Nvol();

  std::ifstream infile(filename.c_str(), std::ios::in | std::ios::binary);
  if (!infile) {
    vout.crucial(m_vl, "%s: file open failed.\n", _function_name);
    exit(EXIT_FAILURE);
  }

  bool do_swap = (is_bigendian() == false);
  if (do_swap) {
    vout.detailed(m_vl, "host endian is little: byte swap performed.\n");
  }

  // read block size
  const int block_size = nin_file;
  char      buf[sizeof(double) * block_size];

  for (int j = 0; j < nex_file; ++j) {
    for (int isite = 0; isite < nvol; ++isite) {
      // read 1 block
      infile.read(buf, sizeof(double) * block_size);

      if (!infile) {
        if (infile.eof()) {  // short file
          vout.crucial(m_vl, "%s: file size too small.\n", _function_name);
        } else {
          vout.crucial(m_vl, "%s: io error.\n", _function_name);
        }

        exit(EXIT_FAILURE);
      }

      if (do_swap) {
        byte_swap(buf, sizeof(double), block_size);
      }

      double *ptr = (double *)buf;

      for (int i = 0; i < nin_file; ++i) {
        int s, t;
        m_format->file_to_field(s, t, i, j);

        u->set(s, isite, t, ptr[i]);
      }
    }
  }

  infile.close();
  vout.detailed(m_vl, "read succeeded.\n");
}


//====================================================================
void FieldIO_Binary_Distributed::write_file(Field *u, string filename_base)
{
  static const char _function_name[] = "FieldIO_Binary_Distributed::write_file";
  string            filename
    = FileUtils::generate_filename("%s.x%dy%dz%dt%d",
                                   filename_base.c_str(),
                                   Communicator::ipe(0),
                                   Communicator::ipe(1),
                                   Communicator::ipe(2),
                                   Communicator::ipe(3));

  int nin_field = u->nin();
  int nex_field = u->nex();

  int nin_file = m_format->nin();
  int nex_file = m_format->nex();

  if ((nin_file == 0) || (nex_file == 0)) {
    nin_file = nin_field;
    nex_file = nex_field;
  }

  int nvol = CommonParameters::Nvol();

  bool do_swap = (is_bigendian() == false);
  if (do_swap) {
    vout.detailed(m_vl, "host endian is little: byte swap performed.\n");
  }

  // write block size
  const int block_size = nin_file;
  char      buf[sizeof(double) * block_size];

  std::ofstream outfile(filename.c_str(), std::ios::out | std::ios::binary);
  if (!outfile) {
    vout.crucial(m_vl, "%s: file open failed\n", _function_name);
    exit(EXIT_FAILURE);
  }

  for (int j = 0; j < nex_file; ++j) {
    for (int isite = 0; isite < nvol; ++isite) {
      double *ptr = (double *)buf;

      for (int i = 0; i < nin_file; ++i) {
        int s, t;
        m_format->file_to_field(s, t, i, j);

        ptr[i] = u->cmp(s, isite, t);
      }

      if (do_swap) {
        byte_swap(buf, sizeof(double), block_size);
      }

      outfile.write(buf, sizeof(double) * block_size);

      if (!outfile) {  // error
        vout.crucial(m_vl, "%s: io error.\n", _function_name);
        exit(EXIT_FAILURE);
      }
    }
  }

  outfile.close();
  vout.detailed(m_vl, "write succeeded.\n");
}


//====================================================================
//============================================================END=====
