/*!
        @file    fieldIO_Binary.h

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate: 2014-04-12 16:31:14 #$

        @version $LastChangedRevision: 1929 $
*/

#ifndef FIELDIO_BINARY_INCLUDED
#define FIELDIO_BINARY_INCLUDED

#include <string>
using std::string;

#include "fieldIO.h"

#include "Field/index_lex.h"
#include "Field/field.h"

//! FieldIO_Binary class BAPI for file I/O of Field data in binary format.

/*!
   Main ingredients of this class BAPI were implemented by T.Yoshie:
     cinput, coutput, big_endian, and byte_swap.
   Incorporated to GaugeConfig_Binary class BAPI family by H.Matsufuru.

   The file format treated in this class BAPI is the same as ILDG
   file format, while not packed to LIME file.
   The endian is big as the definition of ILDG file.
                                        [28 Dec 2011 H.Matsufuru]

   FieldIO_Binary class BAPI provides file I/O of Field data in binary format.
   The inferface is defined in the FieldIO base class BAPI, and this class BAPI
   defines concrete realisation.

   File I/O is performed on the primary node (rank 0 usually), and
   the field data is gathered from/scattered to parallel nodes.
   Extra memory of global field size is temporarily allocated internally
   as workspace, which may restrict operativity for huge lattice.
 */

class BAPI FieldIO_Binary : public FieldIO
{
 public:
  static const std::string class_name;

 public:
  FieldIO_Binary(const IO_Format::Format *format) : FieldIO(format)
  {}

  void read_file(Field *v, const std::string filename) override;
  void write_file(Field *v, const std::string filename) override;
};
#endif
