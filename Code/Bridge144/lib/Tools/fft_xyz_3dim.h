/*!
        @file    $Id:: fft_xyz_3dim.h #$

        @brief

        @author  Yusuke Namekawa (namekawa)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 00:49:46 #$

        @version $LastChangedRevision: 1561 $
*/

#ifndef FFT_XYZ_3DIM_INCLUDED
#define FFT_XYZ_3DIM_INCLUDED

#ifdef USE_FFTWLIB

#include "fft.h"

//! Fast Fourier Transformation in x,y,z directions.

/*!
   This class implements Fast Fourier Transformation
   in x,y,z directions simultaneously.
   NB. FFTW supports parallelization only in z direction.
                                   [06 Jun 2015 Y.Namekawa]
 */

class FFT_xyz_3dim : public FFT
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
  ptrdiff_t m_Nz_p, m_z_start_p;
#endif


 public:
  FFT_xyz_3dim()
  {
    init();
  }

  ~FFT_xyz_3dim()
  {
    tidy_up();
  }

  void fft(Field& field);  // field is overwritten
  void fft(Field& field_out, const Field& field_in);

  void set_parameters(const Parameters&);
  void set_parameters(const std::string str_fft_direction);

  void set_parameter_verboselevel(const Bridge::VerboseLevel vl)
  {
    m_vl = vl;
  }

 private:
  void init();
  void tidy_up();
};
//- #endif of #ifdef USE_FFTWLIB
#endif
#endif
