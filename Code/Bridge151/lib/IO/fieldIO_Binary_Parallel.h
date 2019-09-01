/*!
        @file    fieldIO_Binary_Parallel.h

        @brief

        @author  Tatsumi Aoyama (aoym)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate: 2013-01-22 13:51:53 #$

        @version $LastChangedRevision: 1929 $
*/

#ifndef FIELDIO_BINARY_PARALLEL_INCLUDED
#define FIELDIO_BINARY_PARALLEL_INCLUDED

#include <string>
using std::string;

#ifdef USE_MPI
#include <mpi.h>
#include "Communicator/MPI/communicator_mpi.h"
#endif

#include "fieldIO.h"

#include "bridgeIO.h"
using Bridge::vout;

//! FieldIO_Binary_Parallel class BAPI for file I/O of Field data in binary format using MPI parallel I/O.

/*!
    The file format treated in this class BAPI is the same as ILDG
    file format, while not packed to LIME file.
    The endian is big as the definition of ILDG file.
                                        [28 Dec 2011 H.Matsufuru]

    FieldIO_Binary_Parallel provides file I/O of Field data in binary format.
    File I/O is performed in parallel relying on MPI I/O.
    The interface is defined in the FieldIO base class BAPI, and this class BAPI
    defines concrete realisation.

    Parallel I/O is enabled when USE_MPI option is turned on; otherwise
    this class BAPI is an alias of FieldIO_Binary that provides serial I/O.

    Simultaneous use of BGNET and MPI requires some care:
    prescription by T.Doi was incorporated. [16 Sep 2014 H.Matsufuru]
 */

#ifdef USE_MPI
class BAPI FieldIO_Binary_Parallel : public FieldIO
{
 public:
  static const std::string class_name;

 public:
  FieldIO_Binary_Parallel(const IO_Format::Format *format);
  ~FieldIO_Binary_Parallel();

  void read_file(Field *v, const std::string filename) override;
  void write_file(Field *v, const std::string filename) override;

 private:

  bool m_is_initialized;  //!< check if initialisation is done.

  int m_nvol;
  int m_nin_file;
  int m_nex_file;

  MPI_Datatype m_type_vector;     //!< internal block
  MPI_Datatype m_type_tiled;      //!< subarray of blocks

  int initialize(const Field *v); //!< initialise MPI datatypes for mapping data location to glboal layout.
  int finalize();                 //!< finalise MPI datatypes.
  int clear_layout();             //!< cleanup layout settings for a field element_type.
};

#else

// for Single version, just an alias of FieldIO_Binary.

class BAPI FieldIO_Binary_Parallel : public FieldIO_Binary
{
 public:
  FieldIO_Binary_Parallel(const IO_Format::Format *format) : FieldIO_Binary(format) {}
};
#endif
#endif
