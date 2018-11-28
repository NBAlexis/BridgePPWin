#include "BridgeLib_Private.h"

/*!
        @file    $Id: fieldIO_LIME.cpp #$

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate: 2015-03-24 18:19:19 #$

        @version $LastChangedRevision: 1571 $
*/

#include <fstream>
#include <string.h>

#include "fieldIO_LIME.h"
#include "io_format_gauge.h"

#if USE_ILDG_METADATA
#include "ildg_metadata.h"
#endif

#include "bridgeIO.h"
using Bridge::vout;


#ifdef USE_LIMELIB

//typedef unsigned short int  n_uint16_t;
//typedef unsigned int        n_uint32_t;
//typedef unsigned long int   n_uint64_t;

#ifdef off_t
#undef off_t
#endif
#define off_t    n_uint64_t

static const off_t block_size = 64;

const std::string FieldIO_LIME::class_name = "FieldIO_LIME";

#if USE_ILDG_METADATA

//================================================================
void FieldIO_LIME::check_metadata(const ILDG_Format::Params *params)
{
  // check if config is as expected.
  int Lx = CommonParameters::Lx();
  int Ly = CommonParameters::Ly();
  int Lz = CommonParameters::Lz();
  int Lt = CommonParameters::Lt();

  if (!((params->Lx == Lx) &&
        (params->Ly == Ly) &&
        (params->Lz == Lz) &&
        (params->Lt == Lt))) {
    vout.crucial(m_vl, "Error at %s: lattice size mismatch. config=(%d,%d,%d,%d), expected=(%d,%d,%d,%d).\n",
                 class_name.c_str(),
                 params->Lx, params->Ly, params->Lz, params->Lt,
                 Lx, Ly, Lz, Lt);
    exit(EXIT_FAILURE);
  }

  if (params->prec == 32) {
    vout.detailed(m_vl, "ildg-format: single precision. extend to double\n");
  } else {
    vout.detailed(m_vl, "ildg-format: double precision\n");
  }
}


//================================================================
void FieldIO_LIME::load_metadata(LimeReader *reader, ILDG_Format::Params *params)
{
  off_t nbytes = limeReaderBytes(reader);

  vout.detailed(m_vl, "limeReaderBytes: %lu bytes to read.\n", nbytes);

  off_t nbytes_alloc = nbytes + ((nbytes % block_size) ? (block_size - (nbytes % block_size)) : 0);

  char *buf = (char *)malloc(nbytes_alloc);
  if (buf == 0) {
    vout.crucial(m_vl, "Error at %s: malloc failed.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }
  vout.detailed(m_vl, "allocated %lu bytes\n", nbytes_alloc);

  int status = limeReaderReadData(buf, &nbytes, reader);
  if (status != LIME_SUCCESS) {
    vout.crucial(m_vl, "Error at %s: limeReaderReadData failed.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }

  // dump Metadata
  vout.detailed(m_vl, "%s\n", buf);

  // process ILDG Metadata
  ILDG_Format::Metadata md;
  md.read_from_buffer(buf).extract(params);

  free(buf);
}


// endif of #if USE_ILDG_METADATA
#endif

//================================================================
void FieldIO_LIME::load_lfn(LimeReader *reader)
{
  off_t nbytes = limeReaderBytes(reader);

  vout.detailed(m_vl, "limeReaderBytes: %lu bytes to read.\n", nbytes);

  off_t nbytes_alloc = (nbytes + 1) + (block_size - ((nbytes + 1) % block_size));

  char *buf = (char *)malloc(nbytes_alloc);
  if (!buf) {
    vout.crucial(m_vl, "Error at %s: malloc failed.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }

  vout.detailed(m_vl, "allocated %lu bytes\n", nbytes_alloc);

  int status = limeReaderReadData(buf, &nbytes, reader);
  if (status != LIME_SUCCESS) {
    vout.crucial(m_vl, "Error at %s: limeReaderReadData failed.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }

  vout.detailed(m_vl, "limeReaderReadData: %lu bytes read.\n", nbytes);

  buf[nbytes] = '\0';

  vout.detailed(m_vl, "ildg-data-lfn: %s\n", buf);

  free(buf);
}


//====================================================================
void FieldIO_LIME::load_data(LimeReader *reader, Field *v)
{
  off_t word_size = 8;  //XXX assume 64bit precision

  off_t nbytes = limeReaderBytes(reader);

  vout.detailed(m_vl, "ildg-binary-data: limeReaderBytes: %lu bytes to read.\n", nbytes);

  // allocate memory for whole config data.
  char *buf = (char *)malloc(nbytes);
  if (!buf) {
    vout.crucial(m_vl, "Error at %s: malloc failed.", __func__);
    exit(EXIT_FAILURE);
  }

  int status = limeReaderReadData(buf, &nbytes, reader);
  if (status != LIME_SUCCESS) {
    vout.crucial(m_vl, "Error at %s: malloc failed.", __func__);
    exit(EXIT_FAILURE);
  }

  // adjust byteorder
  if (!FieldIO::is_bigendian()) {
    byte_swap(buf, word_size, nbytes / word_size);
  }

  // reorder and store
  int nin_file = m_format->nin();
  int nex_file = m_format->nex();

  if ((nin_file == 0) || (nex_file == 0)) {
    nin_file = v->nin();
    nex_file = v->nex();
  }

  int lvol = v->nvol();

  double *p = (double *)buf;

  for (int j = 0; j < nex_file; ++j) {
    for (int isite = 0; isite < lvol; ++isite) {
      for (int i = 0; i < nin_file; ++i) {
        int s, t;
        m_format->file_to_field(s, t, i, j);
        v->set(s, isite, t, *p++);
      }
    }
  }

  free(buf);
}


//====================================================================
void FieldIO_LIME::process_file(Field *v, std::string filename)
{
  FILE *fp = fopen(filename.c_str(), "r");

  if (!fp) {
    vout.crucial(m_vl, "Error at %s: fopen failed.", __func__);
    exit(EXIT_FAILURE);
  }

  LimeReader *reader = limeCreateReader(fp);
  if (!reader) {
    vout.crucial(m_vl, "Error at %s: limeCreateReader failed.", __func__);
    exit(EXIT_FAILURE);
  }

#if USE_ILDG_METADATA
  ILDG_Format::Params params;
#endif

  // scan file for records.
  for ( ; ; ) {
    int status = limeReaderNextRecord(reader);
    if (status != LIME_SUCCESS) break;

    const char *t = limeReaderType(reader);

    if (strcmp(t, "ildg-format") == 0) {
#if USE_ILDG_METADATA
      load_metadata(reader, &params);
      check_metadata(&params);
#endif
    } else if (strcmp(t, "ildg-binary-data") == 0) {
      load_data(reader, v);
    } else if (strcmp(t, "ildg-data-lfn") == 0) {
      load_lfn(reader);
    } else {
      vout.detailed(m_vl, "%s: known record %s.\n", __func__, t);
    }
  } // end loop over records.

  limeDestroyReader(reader);

  fclose(fp);
}


//====================================================================
void FieldIO_LIME::read_file(Field *v, std::string filename)
{
  int nin_field = v->nin();
  int nex_field = v->nex();

  int Lvol = CommonParameters::Lvol();

  Field vtmp;

  if (Communicator::is_primary()) {
    vout.detailed(m_vl, "reading gauge configuration from %s", filename.c_str());

    vtmp.reset(nin_field, Lvol, nex_field);

    process_file(&vtmp, filename);
  }

  FieldIO::deliver(v, &vtmp);

  vout.detailed(m_vl, "read successful\n");
}


//================================================================
#if USE_ILDG_METADATA

void FieldIO_LIME::store_metadata(LimeWriter *writer)
{
  // first, write metadata record.
  vout.detailed(m_vl, "%s: write metadata.\n", __func__);

  ILDG_Format::Params params;
  params.Lx   = CommonParameters::Lx();
  params.Ly   = CommonParameters::Ly();
  params.Lz   = CommonParameters::Lz();
  params.Lt   = CommonParameters::Lt();
  params.type = 0;
  params.prec = 8 * sizeof(double);

  ILDG_Format::Metadata md;
  md.store(&params);

  const int buf_size = 4 * 1024;
  char      *buf     = (char *)malloc(buf_size);

  md.write_to_buffer(buf, buf_size);

  off_t nbytes = strlen(buf);

  LimeRecordHeader *h = limeCreateHeader(1 /* MB */, 0 /* ME */, const_cast<char *>("ildg-format"), nbytes);
  limeWriteRecordHeader(h, writer);
  limeDestroyHeader(h);

  limeWriteRecordData(buf, &nbytes, writer);

  free(buf);
}
#endif

//====================================================================
void FieldIO_LIME::store_data(LimeWriter *writer, Field *v,
                              bool mark_begin, bool mark_end)
{
  // second, write binary data.
  vout.detailed(m_vl, "%s: write binary data.\n", __func__);

//  off_t nbytes = sizeof(Field::element_type) * u->size();
  off_t nbytes = sizeof(double) * v->size();

  LimeRecordHeader *h = limeCreateHeader(
    mark_begin ? 1 : 0 /* MB */,
    mark_end ? 1 : 0 /* ME */,
    const_cast<char *>("ildg-binary-data"), nbytes);
  limeWriteRecordHeader(h, writer);
  limeDestroyHeader(h);

  char *buf = (char *)malloc(nbytes);
  if (!buf) {
    vout.crucial(m_vl, "Error at %s: malloc failed.\n", __func__);
    exit(EXIT_FAILURE);
  }

  // reorder and pack to buffer
  int nin_field = v->nin();
  int nex_field = v->nex();

  int nin_file = m_format->nin();
  int nex_file = m_format->nex();

  if ((nin_file == 0) || (nex_file == 0)) {
    nin_file = nin_field;
    nex_file = nex_field;
  }

  int lvol = CommonParameters::Lvol();

  // Field::element_type *p = (Field::element_type *)buf;
  double *p = (double *)buf;

  for (int j = 0; j < nex_file; ++j) {
    for (int isite = 0; isite < lvol; ++isite) {
      for (int i = 0; i < nin_file; ++i) {
        int s, t;
        m_format->file_to_field(s, t, i, j);
        *p++ = v->cmp(s, isite, t);
      }
    }
  }

  // adjust byteorder
  if (!FieldIO::is_bigendian()) {
    // byte_swap(buf, sizeof(Field::element_type), u->size());
    byte_swap(buf, sizeof(double), v->size());
  }

  // store
  limeWriteRecordData(buf, &nbytes, writer);

  vout.detailed(m_vl, "write succeeded.\n");

  free(buf); // added by s.motoki[12.06.05].
}


//====================================================================
void FieldIO_LIME::store_lfn(LimeWriter *writer, std::string lfn_string)
{
  off_t nbytes = lfn_string.size();

  LimeRecordHeader *h = limeCreateHeader(1 /* MB */, 1 /* ME */, const_cast<char *>("ildg-data-lfn"), nbytes);

  limeWriteRecordHeader(h, writer);
  limeDestroyHeader(h);

  // store
  limeWriteRecordData(const_cast<char *>(lfn_string.c_str()), &nbytes, writer);

  vout.detailed(m_vl, "write succeeded.\n");
}


//====================================================================
void FieldIO_LIME::write_file(Field *v, std::string filename)
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

  Field vtmp;
  if (Communicator::is_primary()) {
    vtmp.reset(nin_field, Lvol, nex_field);
  }

  // gather data
  FieldIO::gather(&vtmp, v);

  // reorder

  // dump to file.
  if (Communicator::is_primary()) {
    vout.detailed(m_vl, "writing gauge configuration to %s\n", filename.c_str());

    FILE *fp = fopen(filename.c_str(), "w");
    if (!fp) {
      vout.crucial(m_vl, "Error at %s: cannot open file for write\n", class_name.c_str());
      exit(EXIT_FAILURE);
    }

    LimeWriter *writer = limeCreateWriter(fp);
    if (!writer) {
      vout.crucial(m_vl, "Error at %s: cannot create limeWriter\n", class_name.c_str());
      exit(EXIT_FAILURE);
    }

    // first, write metadata
#ifdef USE_ILDG_METADATA
    store_metadata(writer);
#endif

    // second, write binary data.
    store_data(writer, &vtmp, true, true);

    // if any, write lfn.
//    store_lfn(writer, "lfn://");

    limeDestroyWriter(writer);

    fclose(fp);
  }

  vout.detailed(m_vl, "write succeeded.\n");

  // cleanup.
}


//====================================================================
void FieldIO_LIME::read_file(std::vector<Field *>& vv, const std::string& filename)
{
  if (vv.size() == 0) return;

  int nin_field = vv[0]->nin();
  int nex_field = vv[0]->nex();

  int Lvol = CommonParameters::Lvol();

  Field vtmp;

  if (Communicator::is_primary()) {
    vout.detailed(m_vl, "reading gauge configuration from %s", filename.c_str());

    vtmp.reset(nin_field, Lvol, nex_field);

    FILE *fp = fopen(filename.c_str(), "r");
    if (!fp) {
      vout.crucial(m_vl, "Error at %s: fopen failed.", __func__);
      exit(EXIT_FAILURE);
    }

    LimeReader *reader = limeCreateReader(fp);
    if (!reader) {
      vout.crucial(m_vl, "Error at %s: limeCreateReader failed.", __func__);
      exit(EXIT_FAILURE);
    }

    // primary node reads file and deliver data to the other nodes.

    int idx = 0;  // idx-th field

    // scan file for records.
    for ( ; ; ) {
      int status = limeReaderNextRecord(reader);
      if (status != LIME_SUCCESS) break;

      const char *t = limeReaderType(reader);

      if (strcmp(t, "ildg-format") == 0) {
#if USE_ILDG_METADATA
        ILDG_Format::Params params;
        load_metadata(reader, &params);
        check_metadata(&params);
#endif
      } else if (strcmp(t, "ildg-binary-data") == 0) {
        load_data(reader, &vtmp);

        FieldIO::deliver(vv[idx++], &vtmp);
      } else if (strcmp(t, "ildg-data-lfn") == 0) {
        load_lfn(reader);
      } else {
        vout.detailed(m_vl, "%s: known record %s.\n", __func__, t);
      }
    } // end loop over records.

    limeDestroyReader(reader);

    fclose(fp);
  } else {
    // other nodes wait for data to be delivered from primary node.
    for (int i = 0, n = vv.size(); i < n; ++i) {
      FieldIO::deliver(vv[i], &vtmp);
    }
  }

  vout.detailed(m_vl, "read successful\n");
}


//====================================================================
void FieldIO_LIME::write_file(std::vector<Field *>& vv, const std::string& filename)
{
  if (vv.size() == 0) return;

  int nin_field = vv[0]->nin();
  int nex_field = vv[0]->nex();

  int nin_file = m_format->nin();
  int nex_file = m_format->nex();

  if ((nin_file == 0) || (nex_file == 0)) {
    nin_file = nin_field;
    nex_file = nex_field;
  }

  int Lvol = CommonParameters::Lvol();

  Field vtmp;
  if (Communicator::is_primary()) {
    vtmp.reset(nin_field, Lvol, nex_field);
  }

  FILE       *fp     = NULL;  // only on primary node.
  LimeWriter *writer = NULL;

  // dump to file.
  if (Communicator::is_primary()) {
    vout.detailed(m_vl, "writing gauge configuration to %s\n", filename.c_str());

    fp = fopen(filename.c_str(), "w");
    if (!fp) {
      vout.crucial(m_vl, "Error at %s: cannot open file for write\n", class_name.c_str());
      exit(EXIT_FAILURE);
    }

    writer = limeCreateWriter(fp);
    if (!writer) {
      vout.crucial(m_vl, "Error at %s: cannot create limeWriter\n", class_name.c_str());
      exit(EXIT_FAILURE);
    }

    // first, write metadata
#ifdef USE_ILDG_METADATA
    store_metadata(writer);
#endif
  }

  for (int i = 0, n = vv.size(); i < n; ++i) {
    // gather data
    FieldIO::gather(&vtmp, vv[i]);

    if (Communicator::is_primary()) {
      // reorder

      // second, write binary data.
      store_data(writer, &vtmp, (i == 0), (i == n - 1));
    }
  }

  if (Communicator::is_primary()) {
    // if any, write lfn.
//    store_lfn(writer, "lfn://");

    limeDestroyWriter(writer);

    fclose(fp);
  }

  vout.detailed(m_vl, "write succeeded.\n");

  // cleanup.
}


#else /* USE_LIMELIB */

void FieldIO_LIME::read_file(Field *v, std::string filename) {}
void FieldIO_LIME::write_file(Field *v, std::string filename) {}
#endif /* USE_LIMELIB */

//====================================================================
//============================================================END=====
