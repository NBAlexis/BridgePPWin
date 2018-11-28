/*!
        @file    $Id:: fft.h #$

        @brief

        @author  Yusuke Namekawa (namekawa)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 18:35:38 #$

        @version $LastChangedRevision: 1571 $
*/

#ifndef FFT_INCLUDED
#define FFT_INCLUDED

#ifdef USE_FFTWLIB

#ifdef USE_MPI
#include <fftw3-mpi.h>
#ifdef USE_BGNET
#include "Communicator/BGNET/communicator_bgnet.h"
#else
#include "Communicator/MPI/communicator_mpi.h"
#endif
#else
#include <fftw3.h>
#endif

#ifdef USE_OPENMP
#include "ResourceManager/threadManager_OpenMP.h"
#endif

#include "Field/field.h"
#include "Field/index_lex.h"
#include "Parameters/parameters.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

#ifdef USE_FACTORY
#include "factory.h"
#endif


//! Base class for FFT class family.

/*!
    This class implements a base class of Fast Fourier Transformation.
                               [7 Sep 2015 Y.Namekawa]
 */

class FFT
{
 public:
  FFT() {}

  virtual ~FFT() {}

 private:
  // non-copyable
  FFT(const FFT&);
  FFT& operator=(const FFT&);

 public:
  virtual void fft(Field& field) = 0;  // field is overwritten
  virtual void fft(Field& field_out, const Field& field_in) = 0;

  virtual void set_parameters(const Parameters&) = 0;
  virtual void set_parameters(const std::string str_fft_direction)       = 0;
  virtual void set_parameter_verboselevel(const Bridge::VerboseLevel vl) = 0;

#ifdef USE_FACTORY
 public:
  typedef FFT *(*ProductCreator)();
  typedef FactoryTemplate<FFT, ProductCreator>   Factory;

  static FFT *New(const IdentifierType& subtype)
  {
    ProductCreator p = Factory::Find(subtype);

    return p ? (*p)() : 0;
  }
#endif
};
//- #endif of #ifdef USE_FFTWLIB
#endif
#endif
