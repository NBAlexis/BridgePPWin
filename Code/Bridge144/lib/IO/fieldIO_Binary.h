/*!
        @file    $Id: fieldIO_Binary.h #$

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate: 2014-04-12 16:31:14 #$

        @version $LastChangedRevision: 1561 $
*/

#ifndef FIELDIO_BINARY_INCLUDED
#define FIELDIO_BINARY_INCLUDED

#include <string>
using std::string;

#include "fieldIO.h"

#include "Field/index_lex.h"
#include "Field/field.h"

//! FieldIO_Binary class for file I/O of Field data in binary format.

/*!
   Main ingredients of this class were implemented by T.Yoshie:
     cinput, coutput, big_endian, and byte_swap.
   Incorporated to GaugeConfig_Binary class family by H.Matsufuru.

   The file format treated in this class is the same as ILDG
   file format, while not packed to LIME file.
   The endian is big as the definition of ILDG file.
                                        [28 Dec 2011 H.Matsufuru]

   FieldIO_Binary class provides file I/O of Field data in binary format.
   The inferface is defined in the FieldIO base class, and this class
   defines concrete realisation.

   File I/O is performed on the primary node (rank 0 usually), and
   the field data is gathered from/scattered to parallel nodes.
   Extra memory of global field size is temporarily allocated internally
   as workspace, which may restrict operativity for huge lattice.
 */

class FieldIO_Binary : public FieldIO
{
 public:
  static const std::string class_name;

 public:
  FieldIO_Binary(const IO_Format::Format *format) : FieldIO(format)
  {}

  void read_file(Field *v, string filename);
  void write_file(Field *v, string filename);
};
#endif
