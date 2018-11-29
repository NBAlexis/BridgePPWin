/*!
        @file    $Id: fieldIO_Null.h #$

        @brief

        @author  Tatsumi Aoyama (aoym)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate: 2014-04-12 16:31:14 #$

        @version $LastChangedRevision: 1561 $
*/

#ifndef FIELDIO_NULL_INCLUDED
#define FIELDIO_NULL_INCLUDED

#include <string>
using std::string;

#include "fieldIO.h"
#include "Field/index_lex.h"
#include "Field/field.h"

//! FieldIO_Null class for an analogue of /dev/null

/*!
    FieldIO_Null class provides a dummy I/O entry for
    discarding output, analogous to sending to /dev/null.
                               [26 June 2016 T.Aoyama]
*/

class BAPI FieldIO_Null : public FieldIO
{
 public:
  static const std::string class_name;

 public:
  FieldIO_Null(const IO_Format::Format *format) : FieldIO(format)
  {}

  void read_file(Field *v, string filename);
  void write_file(Field *v, string filename);
};
#endif
