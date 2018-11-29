/*!
        @file    $Id: fieldIO_Text_4x4x4x8.h #$

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate: 2014-10-09 17:20:36 #$

        @version $LastChangedRevision: 1561 $
*/

#ifndef FIELDIO_TEXT_4x4x4x8_INCLUDED
#define FIELDIO_TEXT_4x4x4x8_INCLUDED

#include <string>
using std::string;

#include "fieldIO.h"
#include "Field/index_lex.h"
#include "Field/field.h"

//! FieldIO_Text class for file I/O of Field data in plain text format.

/*!
   This class read 4x4x4x8 in plain text format lattice and
   periodically copy it as a larger size of lattice.
   The configuration is read on the rank 0 node and broadcasted.
   No write method is defined.
                                         [11 Jul 2014 H.Matsufuru]
*/

class BAPI FieldIO_Text_4x4x4x8 : public FieldIO
{
 public:
  static const std::string class_name;

 public:
  FieldIO_Text_4x4x4x8(const IO_Format::Format *format) : FieldIO(format)
  {}

  void read_file(Field *v, string filename);
  void write_file(Field *v, string filename);
};
#endif
