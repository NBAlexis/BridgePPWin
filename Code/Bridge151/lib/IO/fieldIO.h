/*!
        @file    fieldIO.h

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate: 2014-04-12 16:31:14 #$

        @version $LastChangedRevision: 1929 $
*/

#ifndef FIELDIO_INCLUDED
#define FIELDIO_INCLUDED

#include <string>
#include <stdint.h>
#include <vector>

#include "Field/field.h"
#include "Field/index_lex.h"
#include "Parameters/commonParameters.h"
#include "io_format.h"

//! FieldIO class BAPI for file I/O of space-time distributed data.

/*!
   This class BAPI family is used to setup and output the gauge
   configuration.
   This is the base class BAPI, which implements common functions,
   and read/write methods are implemented in subclasses.
   At present, cutting off the gauge field for each node
   and deliver it to the node is implemented in this class BAPI.
   It may be better to separate that to other class BAPI for
   general usage for other field objects.
                                [28 Dec 2011 H.Matsufuru]

   FieldIO class BAPI provides file I/O of Field data.
   Reading and writing gauge configuration from/to file is now
   implemented on top of this construct.

   For scalar data (independent of space-time index), DataIO class BAPI
   hierarchy is provided.

   FieldIO class BAPI is an abstract base class BAPI; various data format
   (text, binary, etc) and I/O scheme (single node or parallel I/O)
   are realised in the respective subclasses.
   This class BAPI also provides common utilities for gathering/distributing
   data over parallel nodes (if any), and converting byte order.
*/

class BAPI FieldIO
{
 public:
  static const std::string class_name;

 private:
  Index_lex idx;

 protected:
  const IO_Format::Format *m_format;

  Bridge::VerboseLevel m_vl;

 public:

  //!< constructor. format specifies data layout on file.
  FieldIO(const IO_Format::Format *format) : m_format(format), m_vl(CommonParameters::Vlevel()) {}
  virtual ~FieldIO() {}

 private:
  // non-copyable
  FieldIO(const FieldIO&);
  FieldIO& operator=(const FieldIO&);

 public:

  virtual void read_file(Field *v, std::string)  = 0; //!< read data from file.
  virtual void write_file(Field *v, std::string) = 0; //!< write data to file.

  virtual void read_file(std::vector<Field *>& vv, const std::string& filename);
  virtual void write_file(std::vector<Field *>& vv, const std::string& filename);

  void deliver(Field *vlocal, Field *vglobal);        //!< distribute data on primary node over parallel nodes.
  void gather(Field *vglobal, Field *vlocal);         //!< gather data on parallel nodes to primary node.

// protected:

  //!< convert byte order. alternative interface.
  static void byte_swap(void *buf, size_t size, size_t nmemb)
  {
    return convert_endian(buf, size, nmemb);
  }

  //!< convert byte order of data, each of whose element has size bytes.
  static void convert_endian(void *buf, size_t size, size_t nmemb);

  //!< check if machine byte order is big-endian.
  static bool is_bigendian();
};
#endif
