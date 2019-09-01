/*!
        @file    fieldIO_None.h

        @brief

        @author  Tatsumi Aoyama (aoym)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate: 2014-04-12 16:31:14 #$

        @version $LastChangedRevision: 1929 $
*/

#ifndef FIELDIO_NONE_INCLUDED
#define FIELDIO_NONE_INCLUDED

#include <string>
using std::string;

#include "fieldIO.h"
#include "Field/index_lex.h"
#include "Field/field.h"

//! FieldIO_None class BAPI for an analogue of /dev/null

/*!
    FieldIO_None class BAPI provides a dummy I/O entry for
    discarding output, analogous to sending to /dev/null.
                               [26 June 2016 T.Aoyama]
*/

class BAPI FieldIO_None : public FieldIO
{
 public:
  static const std::string class_name;

 public:
  FieldIO_None(const IO_Format::Format *format) : FieldIO(format)
  {}

  void read_file(Field *v, const std::string filename) override;
  void write_file(Field *v, const std::string filename) override;
};
#endif
