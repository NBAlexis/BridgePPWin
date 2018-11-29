/*!
        @file    $Id: fieldIO_Text.h #$

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate: 2014-04-12 16:31:14 #$

        @version $LastChangedRevision: 1561 $
*/

#ifndef FIELDIO_TEXT_INCLUDED
#define FIELDIO_TEXT_INCLUDED

#include <string>
using std::string;

#include "fieldIO.h"
#include "Field/index_lex.h"
#include "Field/field.h"

//! FieldIO_Text class for file I/O of Field data in plain text format.

/*!
   The format of the gauge configuration is assumed to
   the same as ILDG standard index order.
                                    [28 Dec 2011 H.Matsufuru]

   FieldIO_Text class provides file I/O of Field data in plain text format.

   File I/O is performed on the primary node (rank 0 usually), and
   the field data is gathered from/scattered to parallel nodes.
*/

class BAPI FieldIO_Text : public FieldIO
{
 public:
  static const std::string class_name;

 public:
  FieldIO_Text(const IO_Format::Format *format) : FieldIO(format)
  {}

  void read_file(Field *v, string filename);
  void write_file(Field *v, string filename);
};
#endif
