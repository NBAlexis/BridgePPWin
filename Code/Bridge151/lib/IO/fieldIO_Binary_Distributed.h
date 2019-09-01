/*!
        @file    fieldIO_Binary_Distributed.h

        @brief

        @author  Tatsumi Aoyama (aoym)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate: 2013-01-22 13:51:53 #$

        @version $LastChangedRevision: 1929 $
*/

#ifndef FIELDIO_BINARY_DISTRIBUTED_INCLUDED
#define FIELDIO_BINARY_DISTRIBUTED_INCLUDED

#include <string>
using std::string;

#include "fieldIO.h"
#include "bridgeIO.h"
using Bridge::vout;

//! FieldIO_Binary_Distributed class BAPI for file I/O of Field data in binary format.

/*!
    The file format treated in this class BAPI is the same as ILDG
    file format, while not packed to LIME file.
    The endian is big as the definition of ILDG file.
                                        [28 Dec 2011 H.Matsufuru]

    FieldIO_Binary_Distributed provides file I/O of Field data in binary format.
    File I/O is performed in parallel by fstream.
    The interface is defined in the FieldIO base class BAPI, and this class BAPI
    defines concrete realisation.
                                        [13 Oct 2015 S.Ueda]
 */

class BAPI FieldIO_Binary_Distributed : public FieldIO
{
 public:
  static const std::string class_name;

 public:
  FieldIO_Binary_Distributed(const IO_Format::Format *format) : FieldIO(format) {}

  void read_file(Field *v, const std::string filename) override;
  void write_file(Field *v, const std::string filename) override;
};
#endif
