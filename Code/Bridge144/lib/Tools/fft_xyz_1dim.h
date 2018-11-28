/*!
        @file    $Id:: fft_xyz_1dim.h #$

        @brief

        @author  Yusuke Namekawa (namekawa)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 00:49:46 #$

        @version $LastChangedRevision: 1561 $
*/

#ifndef FFT_XYZ_1DIM_INCLUDED
#define FFT_XYZ_1DIM_INCLUDED

#ifdef USE_FFTWLIB

#include "fft.h"

//! Fast Fourier Transformation in x,y,z directions.

/*!
   This class implements Fast Fourier Transformation
   in x,y,z directions 1 dim by 1 dim.
   NB. FFTW supports parallelization only in 1 direction.
                                   [06 Jun 2015 Y.Namekawa]
 */

class FFT_xyz_1dim : public FFT
{
 public:
  static const std::string class_name;

 protected:
  Bridge::VerboseLevel m_vl;

 private:
  bool m_is_forward;

  Index_lex m_index;

  fftw_complex *m_in;
  fftw_complex *m_out;

  fftw_plan m_plan;

#ifdef USE_MPI
  ptrdiff_t m_Nsize_in_p, m_start_in_p;
  ptrdiff_t m_Nsize_out_p, m_start_out_p;
#endif


 public:
  FFT_xyz_1dim() {}

  ~FFT_xyz_1dim() {}

  void fft(Field& field);  // field is overwritten
  void fft(Field& field_out, const Field& field_in);

  void set_parameters(const Parameters&);
  void set_parameters(const std::string str_fft_direction);

  void set_parameter_verboselevel(const Bridge::VerboseLevel vl)
  {
    m_vl = vl;
  }
};
//- #endif of #ifdef USE_FFTWLIB
#endif
#endif
