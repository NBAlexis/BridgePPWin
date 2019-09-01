/*!
        @file    fieldIO_Text_4x4x4x8.h

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate: 2014-10-09 17:20:36 #$

        @version $LastChangedRevision: 1929 $
*/

#ifndef FIELDIO_TEXT_4x4x4x8_INCLUDED
#define FIELDIO_TEXT_4x4x4x8_INCLUDED

#include "fieldIO.h"

#include "bridgeIO.h"
using Bridge::vout;

//! FieldIO_Text class BAPI for file I/O of Field data in plain text format.

/*!
   This class BAPI read 4x4x4x8 in plain text format lattice and
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

  void read_file(Field *v, const std::string filename) override;
  void write_file(Field *v, const std::string filename) override;
};
#endif
