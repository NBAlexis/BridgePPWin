/*!
        @file    fieldIO_LIME_Parallel.cpp

        @brief

        @author  Tatsumi Aoyama (aoym)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
 */
#include "BridgeLib_Private.h"
// this code only makes sense in MPI environment.
#ifdef USE_MPI

#include "fieldIO_LIME_Parallel.h"
#include <list>

using Bridge::vout;

const std::string FieldIO_LIME_Parallel::class_name = "FieldIO_LIME_Parallel";

//====================================================================
// private definitions

namespace {
// ILDG metadata
  const char ildg_metadata_template[] =
    "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
    "<ildgFormat xmlns=\"http://www.lqcd.org/ildg\"\n"
    "            xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"\n"
    "            xsi:schemaLocation=\"http://www.lqcd.org/ildg http://www.lqcd.org/ildg/filefmt.xsd\">\n"
    "  <version> 1.0 </version>\n"
    "  <field> su3gauge </field>\n"
    "  <precision> %zu </precision>\n"
    "  <lx> %d </lx> <ly> %d </ly> <lz> %d </lz> <lt> %d </lt>\n"
    "</ildgFormat>";

// LIME header format

#define LIME_MAGIC    ((uint32_t)0x456789ab)
#define MB_MASK       ((uint16_t)0x8000)
#define ME_MASK       ((uint16_t)0x4000)

  struct LIME_header
  {
    uint32_t magic;
    uint16_t version;
    uint16_t bitfield;
    uint64_t length;
    char     type[128];
  };

// local header info
  struct LIME_record_info
  {
    uint64_t offset;
    uint64_t length;
    char     type[128];
  };

  typedef std::list<LIME_record_info>    LIME_message_info;
  typedef std::list<LIME_message_info>   LIME_file_info;

  int traverse(MPI_File& fh, LIME_file_info& file_info);
  int read_lime_header(MPI_File& fh, LIME_header& header);
  int read_lime_content(MPI_File& fh, const MPI_Offset offset, char *buf, const size_t length);
  int find_record_offset(const LIME_file_info& file_info, const char *type, MPI_Offset& pos);

  int report_file_info(const LIME_file_info& file_info);

  size_t write_lime_header(MPI_File& fh, const char *type, const size_t length, const uint16_t flag);
  size_t write_lime_record(MPI_File& fh, const char *type, const char *content, const size_t length, const uint16_t flag);

  Bridge::VerboseLevel vl;
}

//====================================================================
void FieldIO_LIME_Parallel::read_file(Field *u, const std::string filename)
{
  static const char _function_name[] = "FieldIO_LIME_Parallel::read_file";

  initialize();

  MPI_File fh;
  int      ret;

  const int nin_file = m_format->nin();
  const int nex_file = m_format->nex();

  // Field::element_type *buf = new Field::element_type [m_nvol*m_nvector];
  double *buf = new double [nin_file * nex_file * m_nvol];
  if (!buf) {
    vout.crucial(m_vl, "Error at %s: allocate buffer failed.\n", _function_name);
    exit(EXIT_FAILURE);
  }

  ret = MPI_File_open(Communicator_impl::world(), const_cast<char *>(filename.c_str()), MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
  if (ret) {
    vout.crucial(m_vl, "Error at %s: MPI_File_open failed.\n", _function_name);
    exit(EXIT_FAILURE);
  }

  MPI_Offset pos = 0;

  // process records
  if (Communicator::is_primary()) {
    LIME_file_info file_info;

    traverse(fh, file_info);

    report_file_info(file_info);

    if (!find_record_offset(file_info, "ildg-binary-data", pos)) {
      vout.crucial(m_vl, "Error at %s: binary data record not found.\n", _function_name);
      exit(EXIT_FAILURE);
    }
  }

  Communicator::Base::broadcast(sizeof(off_t), &pos, 0);

  ret = MPI_File_set_view(fh, pos, m_type_vector, m_type_tiled, const_cast<char *>("native"), MPI_INFO_NULL);
  if (ret) {
    vout.crucial(m_vl, "Error at %s: MPI_File_set_view failed.\n", _function_name);
    exit(EXIT_FAILURE);
  }

  ret = MPI_File_read_all(fh, (void *)buf, m_nvol * nex_file, m_type_vector, MPI_STATUS_IGNORE);
  if (ret) {
    vout.crucial(m_vl, "Error at %s: MPI_File_read_all failed.\n", _function_name);
    exit(EXIT_FAILURE);
  }

  ret = MPI_File_close(&fh);
  if (ret) {
    vout.crucial(m_vl, "Error at %s: MPI_File_close failed.\n", _function_name);
    exit(EXIT_FAILURE);
  }

  if (!is_bigendian()) {
//    convert_endian(buf, sizeof(Field::element_type), m_nvol*m_nvector);
    convert_endian(buf, sizeof(double), m_nvol * nin_file * nex_file);
  }

  // unpack buffer
  double *p = buf;

  for (int j = 0; j < nex_file; ++j) {
    for (int isite = 0; isite < m_nvol; ++isite) {
      for (int i = 0; i < nin_file; ++i) {
        int s, t;
        m_format->file_to_field(s, t, i, j);

        u->set(s, isite, t, *p++);
      }
    }
  }

  delete [] buf;

  finalize();
}


//====================================================================
void FieldIO_LIME_Parallel::write_file(Field *u, const std::string filename)
{
  static const char _function_name[] = "FieldIO_LIME_Parallel::write_file";

  initialize();

  const int nin_file = m_format->nin();
  const int nex_file = m_format->nex();

  //  Field::element_type *buf = new Field::element_type [m_nvol*m_nvector];
  double *buf = new double [nin_file * nex_file * m_nvol];
  if (!buf) {
    vout.crucial(m_vl, "Error at %s: allocate buffer failed.\n", _function_name);
    exit(EXIT_FAILURE);
  }

  size_t data_length = sizeof(double) * nin_file * nex_file * CommonParameters::Lvol();

  // pack buffer
  double *p = buf;

  for (int j = 0; j < nex_file; ++j) {
    for (int isite = 0; isite < m_nvol; ++isite) {
      for (int i = 0; i < nin_file; ++i) {
        int s, t;
        m_format->file_to_field(s, t, i, j);

        *p++ = u->cmp(s, isite, t);
      }
    }
  }

  if (!is_bigendian()) {
    //    convert_endian(buf, sizeof(Field::element_type), m_nvol*m_nvector);
    convert_endian(buf, sizeof(double), nin_file * nex_file * m_nvol);
  }

  MPI_File fh;
  int      ret;

  ret = MPI_File_open(Communicator_impl::world(), const_cast<char *>(filename.c_str()), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
  if (ret) {
    vout.crucial(m_vl, "Error at %s: MPI_File_open failed.\n", _function_name);
    exit(EXIT_FAILURE);
  }

  MPI_Offset pos = 0;

  if (Communicator::is_primary()) {
    // metadata
    char metadata[2048];
    sprintf(metadata, ildg_metadata_template,
            sizeof(double) * 8, /* bit width */
            CommonParameters::Lx(),
            CommonParameters::Ly(),
            CommonParameters::Lz(),
            CommonParameters::Lt());

    pos += write_lime_record(fh, "ildg-format", metadata, strlen(metadata), MB_MASK);

    // content: write header
    pos += write_lime_header(fh, "ildg-binary-data", data_length, ME_MASK);
  }

  Communicator::Base::broadcast(sizeof(off_t), &pos, 0);

  // content: write data

  ret = MPI_File_set_view(fh, pos, m_type_vector, m_type_tiled, const_cast<char *>("native"), MPI_INFO_NULL);
  if (ret) {
    vout.crucial(m_vl, "Error at %s: MPI_File_set_view failed.\n", _function_name);
    exit(EXIT_FAILURE);
  }

  ret = MPI_File_write_all(fh, (void *)buf, m_nvol * nex_file, m_type_vector, MPI_STATUS_IGNORE);
  if (ret) {
    vout.crucial(m_vl, "Error at %s: MPI_File_write_all failed.\n", _function_name);
    exit(EXIT_FAILURE);
  }

  // content: padding if needed
  if (data_length % 8 > 0) {
    size_t padding_size = (8 - data_length % 8) % 8;

    const char blank[8] = "";
    ret = MPI_File_write_at(fh, static_cast<MPI_Offset>(pos + data_length), const_cast<char *>(blank), static_cast<int>(padding_size), MPI_BYTE, MPI_STATUS_IGNORE);
    if (ret) {
      vout.crucial(vl, "Error at %s: write padding failed.\n", _function_name);
      exit(EXIT_FAILURE);
    }

    vout.general(m_vl, "%s: padding %lu bytes added.\n", _function_name, padding_size);
  }

  ret = MPI_File_close(&fh);
  if (ret) {
    vout.crucial(m_vl, "Error at %s: MPI_File_close failed.\n", _function_name);
    exit(EXIT_FAILURE);
  }

  delete [] buf;

  finalize();
}


//====================================================================
int FieldIO_LIME_Parallel::initialize()
{
  static const char _function_name[] = "FieldIO_LIME_Parallel::initialize";

  // store verbose level to private parameter.
  vl = m_vl;

  if (m_is_initialized) return EXIT_SUCCESS;

  const int nin_file = m_format->nin();
  const int nex_file = m_format->nex();

  if ((nin_file == 0) || (nex_file == 0)) {
    vout.crucial(m_vl, "%s: incompatible file format.\n", _function_name);
    exit(EXIT_FAILURE);
  }


  const int ndim = CommonParameters::Ndim();

  int *global_dims = new int[ndim];
  global_dims[0] = CommonParameters::Lx();
  global_dims[1] = CommonParameters::Ly();
  global_dims[2] = CommonParameters::Lz();
  global_dims[3] = CommonParameters::Lt();

  int *local_dims = new int[ndim];
  local_dims[0] = CommonParameters::Nx();
  local_dims[1] = CommonParameters::Ny();
  local_dims[2] = CommonParameters::Nz();
  local_dims[3] = CommonParameters::Nt();

  m_nvol = 1;
  for (int i = 0; i < ndim; ++i) {
    m_nvol *= local_dims[i];
  }

  int *grid_pos = new int[ndim];
  for (int i = 0; i < ndim; ++i) {
    grid_pos[i] = Communicator::ipe(i);
  }

  int *starts = new int[ndim];
  for (int i = 0; i < ndim; ++i) {
    starts[i] = local_dims[i] * grid_pos[i];
  }

  int ret = 0;

//  MPI_Datatype m_type_vector;
//  ret = MPI_Type_contiguous(sizeof(Field::element_type)*nin_file, MPI_BYTE, &m_type_vector);
  ret = MPI_Type_contiguous(sizeof(double) * nin_file, MPI_BYTE, &m_type_vector);
  if (ret) {
    vout.general(m_vl, "%s: MPI_Type_Contiguous failed.\n", _function_name);
    exit(EXIT_FAILURE);
  }

  ret = MPI_Type_commit(&m_type_vector);
  if (ret) {
    vout.general(m_vl, "%s: MPI_Type_commit failed.\n", _function_name);
    exit(EXIT_FAILURE);
  }

//  MPI_Datatype m_type_tiled;
  ret = MPI_Type_create_subarray(ndim, global_dims, local_dims, starts, MPI_ORDER_FORTRAN, m_type_vector, &m_type_tiled);
  if (ret) {
    vout.general(m_vl, "%s: MPI_Type_create_subarray failed.\n", _function_name);
    exit(EXIT_FAILURE);
  }

  ret = MPI_Type_commit(&m_type_tiled);
  if (ret) {
    vout.general(m_vl, "%s: MPI_Type_commit failed.\n", _function_name);
    exit(EXIT_FAILURE);
  }

  m_is_initialized = true;

  delete [] starts;
  delete [] grid_pos;
  delete [] local_dims;
  delete [] global_dims;

  vout.detailed(m_vl, "FieldIO_LIME_Parallel via MPI I/O initialize done.\n");

  return EXIT_SUCCESS;
}


//====================================================================
int FieldIO_LIME_Parallel::finalize()
{
  static const char _function_name[] = "FieldIO_LIME_Parallel::finalize";

  if (!m_is_initialized) return EXIT_SUCCESS;

  int ret;

  ret = MPI_Type_free(&m_type_tiled);
  if (ret) {
    vout.general(m_vl, "%s: MPI_Type_free failed.\n", _function_name);
    exit(EXIT_FAILURE);
  }

  ret = MPI_Type_free(&m_type_vector);
  if (ret) {
    vout.general(m_vl, "%s: MPI_Type_free failed.\n", _function_name);
    exit(EXIT_FAILURE);
  }

  m_is_initialized = false;

  vout.detailed(m_vl, "%s via MPI I/O finalize done.\n", class_name.c_str());

  return EXIT_SUCCESS;
}


//====================================================================
namespace {
//--------------------------------------------------------------------
  int read_lime_header(MPI_File& fh, LIME_header& header)
  {
    MPI_Status status;

    int ret = MPI_File_read(fh, (void *)&header, sizeof(LIME_header), MPI_BYTE, &status);

    int count;

    MPI_Get_count(&status, MPI_BYTE, &count);
    if (count != sizeof(LIME_header)) {  // data length short. end of file.
      return 1;
    }

    if (ret) {
      vout.crucial(vl, "%s: io error.\n", __func__);
      return -1;
    }

    if (FieldIO::is_bigendian() == false) {
      FieldIO::convert_endian(&header.magic, 4, 1);
      FieldIO::convert_endian(&header.version, 2, 1);
      FieldIO::convert_endian(&header.bitfield, 2, 1);
      FieldIO::convert_endian(&header.length, 8, 1);
    }

    if (header.magic != LIME_MAGIC) {
      vout.crucial(vl, "not lime header.\n");
      return -1;
    }

    return 0;
  }


//====================================================================
  int traverse(MPI_File& fh, LIME_file_info& file_info)
  {
    // go to the beginning of the file
    MPI_File_seek(fh, 0, MPI_SEEK_SET);

    MPI_Offset pos = 0;

    LIME_message_info message_info;

    while (true)
    {
      LIME_header header;
      int         stat = read_lime_header(fh, header);

//     if (stat == -1) {  // bad input
//       return -1;
//     } else if (stat == 1) {  // end of file
//       break;
//     } else if (stat != 0) {
//       // unknown status.
//       return -1;
//     }

//      vout.detailed(vl, "read_lime_header: stat = %d, pos = %lu\n", stat, pos);

      if (stat != 0) {
        break;
      }

      // read ok
      pos += sizeof(LIME_header);

      LIME_record_info record_info;

      memcpy((void *)&record_info, (void *)&header, sizeof(LIME_record_info));
      record_info.offset = pos;

      // padding (0-7)
      size_t padding_size = (8 - header.length % 8) % 8;

      // seek to next record
      // MPI_File_seek(fh, header.length + padding_size, MPI_SEEK_CUR);
      // pos += header.length + padding_size;

      // *** workaround for openmpi-2.x
      pos += header.length + padding_size;
      MPI_File_seek(fh, pos, MPI_SEEK_SET);

      // store record info
      if ((header.bitfield & MB_MASK) == MB_MASK) {
        message_info.clear();
      }

      message_info.push_back(record_info);

      if ((header.bitfield & ME_MASK) == ME_MASK) {
        file_info.push_back(message_info);
//      message_info.clear();
      }
    }

    return 0;
  }


//====================================================================
  int find_record_offset(const LIME_file_info& file_info, const char *type, MPI_Offset& pos)
  {
    bool is_found = false;

    for (LIME_file_info::const_iterator p = file_info.begin(); p != file_info.end(); ++p) {
      for (LIME_message_info::const_iterator q = p->begin(); q != p->end(); ++q) {
        if (strncmp(q->type, type, strlen(type)) == 0) {
          is_found = true;
          pos      = q->offset;
          break;
        }
      }
    }

    return is_found ? 1 : 0;
  }


//====================================================================
  int read_record_content(MPI_File& fh, const LIME_file_info& file_info, const char *type, std::string& content)
  {
    bool             is_found = false;
    LIME_record_info info = LIME_record_info();

    for (LIME_file_info::const_iterator p = file_info.begin(); p != file_info.end(); ++p) {
      for (LIME_message_info::const_iterator q = p->begin(); q != p->end(); ++q) {
        if (strncmp(q->type, type, strlen(type)) == 0) {
          is_found = true;
          info     = *q;
          break;
        }
      }
    }

    if (!is_found) {
      return 0;
    }

    MPI_Status status;
    char       *buf = new char [info.length + 1];
    MPI_File_read_at(fh, info.offset, buf, static_cast<int>(info.length), MPI_BYTE, &status);

    int count;
    MPI_Get_count(&status, MPI_BYTE, &count);

    if (count != info.length) {
      vout.crucial(vl, "Error at %s: read error. content length mismatch.\n", __func__);
      exit(EXIT_FAILURE);
    }

    content = std::string(buf);

    return 1;
  }


//====================================================================
  int report_file_info(const LIME_file_info& file_info)
  {
    Bridge::VerboseLevel vlo = vl;

    vl = Bridge::DETAILED;

    for (LIME_file_info::const_iterator p = file_info.begin(); p != file_info.end(); ++p) {
      vout.detailed(vl, "Message:\n");

      for (LIME_message_info::const_iterator q = p->begin(); q != p->end(); ++q) {
        vout.detailed(vl, "\tRecord:\n");
        vout.detailed(vl, "\t\toffset = %lu\n", q->offset);
        vout.detailed(vl, "\t\tsize   = %lu\n", q->length);
        vout.detailed(vl, "\t\ttype   = %s\n", q->type);
      }
    }

    vl = vlo;

    return 0;
  }


//====================================================================
  size_t write_lime_header(MPI_File& fh, const char *type, const size_t length, const uint16_t flag)
  {
    LIME_header header;

    memset(&header, 0, sizeof(LIME_header));

    header.magic    = LIME_MAGIC;
    header.version  = (uint16_t)1;
    header.bitfield = flag;
    strncpy(header.type, type, 128);
    header.length = length;

    if (FieldIO::is_bigendian() == false) {
      FieldIO::convert_endian(&header.magic, 4, 1);
      FieldIO::convert_endian(&header.version, 2, 1);
      FieldIO::convert_endian(&header.bitfield, 2, 1);
      FieldIO::convert_endian(&header.length, 8, 1);
    }

    MPI_Status status;
    int        ret = MPI_File_write(fh, (void *)&header, sizeof(LIME_header), MPI_BYTE, &status);

    if (ret) {
      vout.crucial(vl, "%s: write header failed.\n", __func__);
      return 0;
    }

    return sizeof(LIME_header);  // length written.
  }


//====================================================================
  size_t write_lime_record(MPI_File& fh, const char *type, const char *content, const size_t length, const uint16_t flag)
  {
    const char blank[8] = "";

    if (write_lime_header(fh, type, length, flag) == 0) {
      return 0;
    }

    const size_t padding_size = (8 - length % 8) % 8;

    MPI_Status status;
    int        ret = MPI_File_write(fh, const_cast<char *>(content), static_cast<int>(length), MPI_BYTE, &status);
    if (ret) {
      vout.crucial(vl, "%s: write content failed.\n", __func__);
      return 0;
    }

    if (padding_size > 0) {
      ret = MPI_File_write(fh, const_cast<char *>(blank), static_cast<int>(padding_size), MPI_BYTE, &status);
      if (ret) {
        vout.crucial(vl, "%s: write padding failed.\n", __func__);
        return 0;
      }
    }

    return sizeof(LIME_header) + length + padding_size;
  }


//--------------------------------------------------------------------
} // unnamed namespace
#endif  /* USE_MPI */

//====================================================================
//============================================================END=====
