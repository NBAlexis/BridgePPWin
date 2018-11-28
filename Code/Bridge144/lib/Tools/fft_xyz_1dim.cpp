#include "BridgeLib_Private.h"

/*!
        @file    $Id:: fft_xyz_1dim.cpp #$

        @brief

        @author  Yusuke Namekawa (namekawa)
                 $LastChangedBy: namekawa $

        @date    $LastChangedDate:: 2017-03-31 16:23:13 #$

        @version $LastChangedRevision: 1604 $
*/

#ifdef USE_FFTWLIB

#include "fft_xyz_1dim.h"

#ifdef USE_FACTORY
namespace {
  FFT *create_object()
  {
    return new FFT_xyz_1dim();
  }


  bool init = FFT::Factory::Register("FFT_xyz_1dim", create_object);
}
#endif

const std::string FFT_xyz_1dim::class_name = "FFT_xyz_1dim";

//====================================================================
void FFT_xyz_1dim::set_parameters(const Parameters& params)
{
  const std::string str_vlevel = params.get_string("verbose_level");

  m_vl = vout.set_verbose_level(str_vlevel);

  //- fetch and check input parameters
  std::string str_fft_direction;

  int err = 0;
  err += params.fetch_string("FFT_direction", str_fft_direction);

  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameter not found.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }

  set_parameters(str_fft_direction);
}


//====================================================================
void FFT_xyz_1dim::set_parameters(const std::string str_fft_direction)
{
  //- print input parameters
  vout.general(m_vl, "%s:\n", class_name.c_str());
  vout.general(m_vl, "  FFT_direction = %s\n", str_fft_direction.c_str());

  //- range check

  //- store values
  if (str_fft_direction == "Forward") {
    m_is_forward = true;
  } else if (str_fft_direction == "Backward") {
    m_is_forward = false;
  } else {
    vout.crucial(m_vl, "Error at %s: unsupported FFT direction \"%s\"\n", class_name.c_str(), str_fft_direction.c_str());
    exit(EXIT_FAILURE);
  }
}


//====================================================================
void FFT_xyz_1dim::fft(Field& field)
{
  const int Ndim = CommonParameters::Ndim();

  //- global lattice size
  std::vector<int> Lsize(Ndim - 1);

  Lsize[0] = CommonParameters::Lx();
  Lsize[1] = CommonParameters::Ly();
  Lsize[2] = CommonParameters::Lz();

  const int Lt   = CommonParameters::Lt();
  const int Lxyz = Lsize[0] * Lsize[1] * Lsize[2];

  //- local lattice size
  std::vector<int> Nsize(Ndim - 1);
  Nsize[0] = CommonParameters::Nx();
  Nsize[1] = CommonParameters::Ny();
  Nsize[2] = CommonParameters::Nz();

  const int Nin = field.nin();
  const int Nex = field.nex();

  std::vector<int> NPE(Ndim);
  NPE[0] = CommonParameters::NPEx();
  NPE[1] = CommonParameters::NPEy();
  NPE[2] = CommonParameters::NPEz();
  NPE[3] = CommonParameters::NPEt();

  const int NPE_xyt = NPE[0] * NPE[1] * NPE[3];
  const int NPE_xzt = NPE[0] * NPE[2] * NPE[3];
  const int NPE_yzt = NPE[1] * NPE[2] * NPE[3];

  if ((NPE_xyt != 1) && (NPE_xzt != 1) && (NPE_yzt != 1)) {
    vout.crucial(m_vl, "Error at %s: FFTW supports parallelization only in 1 direction of x,y,z.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

  bool is_allocated     = false;
  int  Lsize_Nsize_prev = 0;

#ifdef USE_OPENMP
  int threads_ok = fftw_init_threads();
#endif
#ifdef USE_MPI
  fftw_mpi_init();
#endif


  //- xyz loop
  for (int i_dim = 0; i_dim < Ndim - 1; ++i_dim) {
    //- allocate m_in,out = m_in,out[Lsize]
    if (Lsize_Nsize_prev != Lsize[i_dim] * Nsize[i_dim]) {
      if (is_allocated) {
        fftw_free(m_in);
        fftw_free(m_out);
        fftw_destroy_plan(m_plan);
      }

      is_allocated = false;
    }
    Lsize_Nsize_prev = Lsize[i_dim] * Nsize[i_dim];

    if (!is_allocated) {
      if (Lsize[i_dim] == Nsize[i_dim]) {
        size_t fftw_size = sizeof(fftw_complex) * Lsize[i_dim];
        m_in  = (fftw_complex *)fftw_malloc(fftw_size);
        m_out = (fftw_complex *)fftw_malloc(fftw_size);

        if (!m_in || !m_out) {
          vout.crucial(m_vl, "Error at %s: failed to allocate memory %d [Byte].\n",
                       class_name.c_str(), (int)fftw_size);
          exit(EXIT_FAILURE);
        }
      } else {
#ifdef USE_MPI
        ptrdiff_t Lsize_p = Lsize[i_dim];
        ptrdiff_t fftw_size_p;

        if (m_is_forward) {
          fftw_size_p = fftw_mpi_local_size_1d(Lsize_p,
                                               Communicator_impl::world(),
                                               FFTW_FORWARD, FFTW_ESTIMATE,
                                               &m_Nsize_in_p, &m_start_in_p,
                                               &m_Nsize_out_p, &m_start_out_p);
        } else {
          fftw_size_p = fftw_mpi_local_size_1d(Lsize_p,
                                               Communicator_impl::world(),
                                               FFTW_BACKWARD, FFTW_ESTIMATE,
                                               &m_Nsize_in_p, &m_start_in_p,
                                               &m_Nsize_out_p, &m_start_out_p);
        }

        m_in  = fftw_alloc_complex(fftw_size_p);
        m_out = fftw_alloc_complex(fftw_size_p);

        if (!m_in || !m_out) {
          vout.crucial(m_vl, "Error at %s: failed to allocate memory %d [Byte].\n",
                       class_name.c_str(), (int)fftw_size_p);
          exit(EXIT_FAILURE);
        }
#endif
      }

      is_allocated = true;
    }


    //- setup FFTW plan
#ifdef USE_OPENMP
    int Nthread = ThreadManager_OpenMP::get_num_threads();
    fftw_plan_with_nthreads(Nthread);
#endif

    if (Lsize[i_dim] == Nsize[i_dim]) {
      if (m_is_forward) {
        m_plan = fftw_plan_dft_1d(Lsize[i_dim],
                                  m_in, m_out,
                                  FFTW_FORWARD, FFTW_ESTIMATE);
      } else {
        m_plan = fftw_plan_dft_1d(Lsize[i_dim],
                                  m_in, m_out,
                                  FFTW_BACKWARD, FFTW_ESTIMATE);
      }
    } else {
#ifdef USE_MPI
      ptrdiff_t Lsize_p = Lsize[i_dim];

      if (m_is_forward) {
        m_plan = fftw_mpi_plan_dft_1d(Lsize_p,
                                      m_in, m_out,
                                      Communicator_impl::world(),
                                      FFTW_FORWARD, FFTW_ESTIMATE);
      } else {
        m_plan = fftw_mpi_plan_dft_1d(Lsize_p,
                                      m_in, m_out,
                                      Communicator_impl::world(),
                                      FFTW_BACKWARD, FFTW_ESTIMATE);
      }
#endif
    }


    // ####  Execution main part  ####
    //- Nin is devided by 2, because of complex(i.e. real and imag)
    for (int in2 = 0; in2 < Nin / 2; ++in2) {
      for (int t_global = 0; t_global < Lt; t_global++) {
        for (int ex = 0; ex < Nex; ++ex) {
          if (i_dim == 0) {
            for (int z = 0; z < Nsize[2]; z++) {
              for (int y = 0; y < Nsize[1]; y++) {
                for (int xyz_local = 0; xyz_local < Nsize[i_dim]; xyz_local++) {
                  int isite  = m_index.site(xyz_local, y, z, t_global);
                  int i_real = 2 * in2;
                  int i_imag = 2 * in2 + 1;

                  m_in[xyz_local][0] = field.cmp(i_real, isite, ex);
                  m_in[xyz_local][1] = field.cmp(i_imag, isite, ex);
                }


                fftw_execute(m_plan);


                for (int xyz_local = 0; xyz_local < Nsize[i_dim]; xyz_local++) {
                  int isite  = m_index.site(xyz_local, y, z, t_global);
                  int i_real = 2 * in2;
                  int i_imag = 2 * in2 + 1;

                  field.set(i_real, isite, ex, m_out[xyz_local][0]);
                  field.set(i_imag, isite, ex, m_out[xyz_local][1]);
                }
              }
            }
          } else if (i_dim == 1) {
            for (int x = 0; x < Nsize[0]; x++) {
              for (int z = 0; z < Nsize[2]; z++) {
                for (int xyz_local = 0; xyz_local < Nsize[i_dim]; xyz_local++) {
                  int isite  = m_index.site(x, xyz_local, z, t_global);
                  int i_real = 2 * in2;
                  int i_imag = 2 * in2 + 1;

                  m_in[xyz_local][0] = field.cmp(i_real, isite, ex);
                  m_in[xyz_local][1] = field.cmp(i_imag, isite, ex);
                }


                fftw_execute(m_plan);


                for (int xyz_local = 0; xyz_local < Nsize[i_dim]; xyz_local++) {
                  int isite  = m_index.site(x, xyz_local, z, t_global);
                  int i_real = 2 * in2;
                  int i_imag = 2 * in2 + 1;

                  field.set(i_real, isite, ex, m_out[xyz_local][0]);
                  field.set(i_imag, isite, ex, m_out[xyz_local][1]);
                }
              }
            }
          } else if (i_dim == 2) {
            for (int y = 0; y < Nsize[1]; y++) {
              for (int x = 0; x < Nsize[0]; x++) {
                for (int xyz_local = 0; xyz_local < Nsize[i_dim]; xyz_local++) {
                  int isite  = m_index.site(x, y, xyz_local, t_global);
                  int i_real = 2 * in2;
                  int i_imag = 2 * in2 + 1;

                  m_in[xyz_local][0] = field.cmp(i_real, isite, ex);
                  m_in[xyz_local][1] = field.cmp(i_imag, isite, ex);
                }


                fftw_execute(m_plan);


                for (int xyz_local = 0; xyz_local < Nsize[i_dim]; xyz_local++) {
                  int isite  = m_index.site(x, y, xyz_local, t_global);
                  int i_real = 2 * in2;
                  int i_imag = 2 * in2 + 1;

                  field.set(i_real, isite, ex, m_out[xyz_local][0]);
                  field.set(i_imag, isite, ex, m_out[xyz_local][1]);
                }
              }
            }
          }
          //- end of if( i_dim == 0 ){
        }
      }
    }
    //- end of outer loops
  }
  //- end of for (int i_dim = 0; i_dim < Ndim - 1; ++i_dim)

  //- normailzation for FFTW_BACKWARD
  if (!m_is_forward) {
    scal(field, 1.0 / Lxyz);
  }

  //- tidy up
  if (is_allocated) {
    fftw_free(m_in);
    fftw_free(m_out);
    fftw_destroy_plan(m_plan);
  }
}


//====================================================================
void FFT_xyz_1dim::fft(Field& field_out, const Field& field_in)
{
  const int Ndim = CommonParameters::Ndim();

  //- global lattice size
  std::vector<int> Lsize(Ndim - 1);

  Lsize[0] = CommonParameters::Lx();
  Lsize[1] = CommonParameters::Ly();
  Lsize[2] = CommonParameters::Lz();

  const int Lt   = CommonParameters::Lt();
  const int Lxyz = Lsize[0] * Lsize[1] * Lsize[2];

  //- local lattice size
  std::vector<int> Nsize(Ndim - 1);
  Nsize[0] = CommonParameters::Nx();
  Nsize[1] = CommonParameters::Ny();
  Nsize[2] = CommonParameters::Nz();

  const int Nin = field_in.nin();
  const int Nex = field_in.nex();

  std::vector<int> NPE(Ndim);
  NPE[0] = CommonParameters::NPEx();
  NPE[1] = CommonParameters::NPEy();
  NPE[2] = CommonParameters::NPEz();
  NPE[3] = CommonParameters::NPEt();

  const int NPE_xyt = NPE[0] * NPE[1] * NPE[3];
  const int NPE_xzt = NPE[0] * NPE[2] * NPE[3];
  const int NPE_yzt = NPE[1] * NPE[2] * NPE[3];

  if ((NPE_xyt != 1) && (NPE_xzt != 1) && (NPE_yzt != 1)) {
    vout.crucial(m_vl, "Error at %s: FFTW supports parallelization only in 1 direction.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

  bool is_allocated     = false;
  int  Lsize_Nsize_prev = 0;

  field_out = field_in;

#ifdef USE_OPENMP
  int threads_ok = fftw_init_threads();
#endif
#ifdef USE_MPI
  fftw_mpi_init();
#endif


  //- xyz loop
  for (int i_dim = 0; i_dim < Ndim - 1; ++i_dim) {
    //- allocate m_in,out = m_in,out[Lsize]
    if (Lsize_Nsize_prev != Lsize[i_dim] * Nsize[i_dim]) {
      if (is_allocated) {
        fftw_free(m_in);
        fftw_free(m_out);
        fftw_destroy_plan(m_plan);
      }

      is_allocated = false;
    }
    Lsize_Nsize_prev = Lsize[i_dim] * Nsize[i_dim];

    if (!is_allocated) {
      if (Lsize[i_dim] == Nsize[i_dim]) {
        size_t fftw_size = sizeof(fftw_complex) * Lsize[i_dim];
        m_in  = (fftw_complex *)fftw_malloc(fftw_size);
        m_out = (fftw_complex *)fftw_malloc(fftw_size);

        if (!m_in || !m_out) {
          vout.crucial(m_vl, "Error at %s: failed to allocate memory %d [Byte].\n",
                       class_name.c_str(), (int)fftw_size);
          exit(EXIT_FAILURE);
        }
      } else {
#ifdef USE_MPI
        ptrdiff_t Lsize_p = Lsize[i_dim];
        ptrdiff_t fftw_size_p;

        if (m_is_forward) {
          fftw_size_p = fftw_mpi_local_size_1d(Lsize_p,
                                               Communicator_impl::world(),
                                               FFTW_FORWARD, FFTW_ESTIMATE,
                                               &m_Nsize_in_p, &m_start_in_p,
                                               &m_Nsize_out_p, &m_start_out_p);
        } else {
          fftw_size_p = fftw_mpi_local_size_1d(Lsize_p,
                                               Communicator_impl::world(),
                                               FFTW_BACKWARD, FFTW_ESTIMATE,
                                               &m_Nsize_in_p, &m_start_in_p,
                                               &m_Nsize_out_p, &m_start_out_p);
        }

        m_in  = fftw_alloc_complex(fftw_size_p);
        m_out = fftw_alloc_complex(fftw_size_p);

        if (!m_in || !m_out) {
          vout.crucial(m_vl, "Error at %s: failed to allocate memory %d [Byte].\n",
                       class_name.c_str(), (int)fftw_size_p);
          exit(EXIT_FAILURE);
        }
#endif
      }

      is_allocated = true;
    }


    //- setup FFTW plan
#ifdef USE_OPENMP
    int Nthread = ThreadManager_OpenMP::get_num_threads();
    fftw_plan_with_nthreads(Nthread);
#endif

    if (Lsize[i_dim] == Nsize[i_dim]) {
      if (m_is_forward) {
        m_plan = fftw_plan_dft_1d(Lsize[i_dim],
                                  m_in, m_out,
                                  FFTW_FORWARD, FFTW_ESTIMATE);
      } else {
        m_plan = fftw_plan_dft_1d(Lsize[i_dim],
                                  m_in, m_out,
                                  FFTW_BACKWARD, FFTW_ESTIMATE);
      }
    } else {
#ifdef USE_MPI
      ptrdiff_t Lsize_p = Lsize[i_dim];

      if (m_is_forward) {
        m_plan = fftw_mpi_plan_dft_1d(Lsize_p,
                                      m_in, m_out,
                                      Communicator_impl::world(),
                                      FFTW_FORWARD, FFTW_ESTIMATE);
      } else {
        m_plan = fftw_mpi_plan_dft_1d(Lsize_p,
                                      m_in, m_out,
                                      Communicator_impl::world(),
                                      FFTW_BACKWARD, FFTW_ESTIMATE);
      }
#endif
    }


    // ####  Execution main part  ####
    //- Nin is devided by 2, because of complex(i.e. real and imag)
    for (int in2 = 0; in2 < Nin / 2; ++in2) {
      for (int t_global = 0; t_global < Lt; t_global++) {
        for (int ex = 0; ex < Nex; ++ex) {

          if (i_dim == 0) {
	    for (int z = 0; z < Nsize[2]; z++) {
              for (int y = 0; y < Nsize[1]; y++) {
                for (int xyz_local = 0; xyz_local < Nsize[i_dim]; xyz_local++) {
                  int isite  = m_index.site(xyz_local, y, z, t_global);
                  int i_real = 2 * in2;
                  int i_imag = 2 * in2 + 1;

                  m_in[xyz_local][0] = field_out.cmp(i_real, isite, ex);
                  m_in[xyz_local][1] = field_out.cmp(i_imag, isite, ex);
                }


                fftw_execute(m_plan);


                for (int xyz_local = 0; xyz_local < Nsize[i_dim]; xyz_local++) {
                  int isite  = m_index.site(xyz_local, y, z, t_global);
                  int i_real = 2 * in2;
                  int i_imag = 2 * in2 + 1;

                  field_out.set(i_real, isite, ex, m_out[xyz_local][0]);
                  field_out.set(i_imag, isite, ex, m_out[xyz_local][1]);
                }
              }
            }

          } else if (i_dim == 1) {
	    for (int x = 0; x < Nsize[0]; x++) {
	      for (int z = 0; z < Nsize[2]; z++) {
                for (int xyz_local = 0; xyz_local < Nsize[i_dim]; xyz_local++) {
                  int isite  = m_index.site(x, xyz_local, z, t_global);
                  int i_real = 2 * in2;
                  int i_imag = 2 * in2 + 1;

                  m_in[xyz_local][0] = field_out.cmp(i_real, isite, ex);
                  m_in[xyz_local][1] = field_out.cmp(i_imag, isite, ex);
                }


                fftw_execute(m_plan);


                for (int xyz_local = 0; xyz_local < Nsize[i_dim]; xyz_local++) {
                  int isite  = m_index.site(x, xyz_local, z, t_global);
                  int i_real = 2 * in2;
                  int i_imag = 2 * in2 + 1;

                  field_out.set(i_real, isite, ex, m_out[xyz_local][0]);
                  field_out.set(i_imag, isite, ex, m_out[xyz_local][1]);
                }
              }
            }

          } else if (i_dim == 2) {
	    for (int y = 0; y < Nsize[1]; y++) {
              for (int x = 0; x < Nsize[0]; x++) {
                for (int xyz_local = 0; xyz_local < Nsize[i_dim]; xyz_local++) {
                  int isite  = m_index.site(x, y, xyz_local, t_global);
                  int i_real = 2 * in2;
                  int i_imag = 2 * in2 + 1;

                  m_in[xyz_local][0] = field_out.cmp(i_real, isite, ex);
                  m_in[xyz_local][1] = field_out.cmp(i_imag, isite, ex);
                }


                fftw_execute(m_plan);


                for (int xyz_local = 0; xyz_local < Nsize[i_dim]; xyz_local++) {
                  int isite  = m_index.site(x, y, xyz_local, t_global);
                  int i_real = 2 * in2;
                  int i_imag = 2 * in2 + 1;

                  field_out.set(i_real, isite, ex, m_out[xyz_local][0]);
                  field_out.set(i_imag, isite, ex, m_out[xyz_local][1]);
                }
              }
            }
          }
          //- end of if( i_dim == 0 ){
        }
      }
    }
    //- end of outer loops
  }
  //- end of for (int i_dim = 0; i_dim < Ndim - 1; ++i_dim)

  //- normailzation for FFTW_BACKWARD
  if (!m_is_forward) {
    scal(field_out, 1.0 / Lxyz);
  }

  //- tidy up
  if (is_allocated) {
    fftw_free(m_in);
    fftw_free(m_out);
    fftw_destroy_plan(m_plan);
  }
}


//==========================================================
//==================================================END=====
#endif
