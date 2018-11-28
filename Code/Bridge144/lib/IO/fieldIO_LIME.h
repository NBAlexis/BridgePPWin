/*!
        @file    $Id: fieldIO_LIME.h #$

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate: 2014-10-09 17:20:36 #$

        @version $LastChangedRevision: 1561 $
*/

#ifndef FIELDIO_LIME_INCLUDED
#define FIELDIO_LIME_INCLUDED

#include <string>
using std::string;

#ifdef USE_LIMELIB
extern "C" {
#include <lime.h>
}
#endif

#include "fieldIO.h"
//#include "ildg_metadata.h"

#include "bridgeIO.h"

#undef USE_ILDG_METADATA

//! FieldIO_LIME class for file I/O of Field data in LIME format.

/*!
   The file format treated in this class is the same as ILDG
   file format, while not packed to LIME file.
   The endian is big as the definition of ILDG file.
                                        [28 Dec 2011 H.Matsufuru]

   FieldIO_LIME class provides file I/O of Field data in LIME file format.
   This class uses lime library to handle LIME archive.

   File I/O is performed on the primary node (rank 0 usually), and
   the field data is gathered from/scattered to parallel nodes.

   This class is available when USE_LIMELIB is defined; otherwise
   instantiation fails at runtime.

   Handling of Metadata is not supported in the current version.

   N.B. The present implementation much assumes records to be ILDG data.
*/

#ifdef USE_LIMELIB

class FieldIO_LIME : public FieldIO
{
 public:
  static const std::string class_name;

 public:
  FieldIO_LIME(const IO_Format::Format *format) : FieldIO(format)
  {
  }

  void read_file(Field *v, std::string filename);
  void write_file(Field *v, std::string filename);

  void read_file(std::vector<Field *>& vv, const std::string& filename);
  void write_file(std::vector<Field *>& vv, const std::string& filename);

 private:

  void process_file(Field *v, std::string filename);

  void load_data(LimeReader *reader, Field *v);
  void store_data(LimeWriter *writer, Field *v, bool mark_begin, bool mark_end);

#ifdef USE_ILDG_METADATA
  void load_metadata(LimeReader *reader, ILDG_Format::Params *params);
  void check_metadata(const ILDG_Format::Params *params);
  void store_metadata(LimeWriter *writer);
#endif

  void load_lfn(LimeReader *reader);
  void store_lfn(LimeWriter *reader, std::string lfn_string);
};

#else

class FieldIO_LIME : public FieldIO
{
 public:
  FieldIO_LIME(const IO_Format::Format *format) : FieldIO(format)
  {
    Bridge::vout.crucial("Error at FieldIO_LIME: USE_LIMELIB is not defined.\n");

    exit(EXIT_FAILURE);
  }

  void read_file(Field *v, std::string filename);
  void write_file(Field *v, std::string filename);
};
#endif
#endif /* GAUGECONFIG_ILDG_INCLUDED */
