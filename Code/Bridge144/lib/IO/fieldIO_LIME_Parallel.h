/*!
        @file    $Id:: fieldIO_LIME_Parallel.h #$

        @brief

        @author  Tatsumi Aoyama (aoym)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 00:49:46 #$

        @version $LastChangedRevision: 1561 $
 */

#ifndef FIELDIO_LIME_PARALLEL_INCLUDED
#define FIELDIO_LIME_PARALLEL_INCLUDED

#include <string>
using std::string;

#ifdef USE_MPI
#include <mpi.h>
#ifdef USE_BGNET
#include "Communicator/BGNET/communicator_bgnet.h"
#else
#include "Communicator/MPI/communicator_mpi.h"
#endif
#endif

#include "fieldIO.h"
#include "Field/field.h"
#include "bridgeIO.h"

//! FieldIO_LIME_Parallel class for file I/O of Field data in LIME format.

/*!
   The file format treated in this class is the same as ILDG
   file format, while not packed to LIME file.
   The endian is big as the definition of ILDG file.
                                        [28 Dec 2011 H.Matsufuru]

   FieldIO_LIME_Parallel class provides file I/O of Field data in LIME file format.
   This class uses lime library to handle LIME archive.

   File I/O is performed on the primary node (rank 0 usually), and
   the field data is gathered from/scattered to parallel nodes.

   This class is available when USE_LIMELIB is defined; otherwise
   instantiation fails at runtime.

   Handling of Metadata is not supported in the current version.

   N.B. The present implementation much assumes records to be ILDG data.

   Simultaneous use of BGNET and MPI requires some care:
   prescription by T.Doi was incorporated. [16 Sep 2014 H.Matsufuru]
*/

#ifdef USE_MPI

class FieldIO_LIME_Parallel : public FieldIO
{
 public:
  static const std::string class_name;

 public:
  FieldIO_LIME_Parallel(const IO_Format::Format *format) : FieldIO(format), m_is_initialized(false) {}

  void read_file(Field *v, std::string filename);
  void write_file(Field *v, std::string filename);

 private:

  bool m_is_initialized;  //!< check if initialisation is done.

  int m_nvol;

  MPI_Datatype m_type_vector; //!< internal block
  MPI_Datatype m_type_tiled;  //!< subarray of blocks

  int initialize();           //!< initialise MPI datatypes for mapping data location to glboal layout.
  int finalize();             //!< finalize MPI datatypes.
};

#else

// for Single version, just an alias of FieldIO_LIME.

class FieldIO_LIME_Parallel : public FieldIO_LIME
{
 public:
  FieldIO_LIME_Parallel(const IO_Format::Format *format) : FieldIO_LIME(format) {}
};
#endif /* USE_MPI */
#endif /* GAUGECONFIG_ILDG_INCLUDED */
