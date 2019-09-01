/*!
        @file    wilsonLoop.cpp

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/
#include "BridgeLib_Private.h"
#include "wilsonLoop.h"

const std::string WilsonLoop::class_name = "WilsonLoop";

//====================================================================
void WilsonLoop::set_parameters(const Parameters& params)
{
  m_filename_output = params.get_string("filename_output");
  if (m_filename_output.empty()) {
    m_filename_output = "stdout";
  }

  const string str_vlevel = params.get_string("verbose_level");
  m_vl = vout.set_verbose_level(str_vlevel);

  //- fetch and check input parameters
  int Nspc_size, Ntmp_size, Ntype;

  int err = 0;
  err += params.fetch_int("max_spatial_loop_size", Nspc_size);
  err += params.fetch_int("max_temporal_loop_size", Ntmp_size);
  err += params.fetch_int("number_of_loop_type", Ntype);

  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameter not found.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }


  set_parameters(Nspc_size, Ntmp_size, Ntype);
}


//====================================================================
void WilsonLoop::set_parameters(const int Nspc_size, const int Ntmp_size, const int Ntype)
{
  //- print input parameters
  vout.general(m_vl, "Wilson loop measurement:\n");
  vout.general(m_vl, "  Nspc_size = %d\n", Nspc_size);
  vout.general(m_vl, "  Ntmp_size = %d\n", Ntmp_size);
  vout.general(m_vl, "  Ntype     = %d\n", Ntype);

  //- range check
  int err = 0;
  err += ParameterCheck::non_negative(Nspc_size);
  err += ParameterCheck::non_negative(Ntmp_size);
  err += ParameterCheck::non_negative(Ntype);

  //! The following setting explicitly depends on the definition
  //! of unit vectors.
  if (Ntype > 6) ++err;

  if (err) {
    vout.crucial(m_vl, "Error at %s: parameter range check failed.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }

  //- store values
  m_Nspc_size = Nspc_size;
  m_Ntmp_size = Ntmp_size;
  m_Ntype     = Ntype;


  //- post-process

  //- set up internal data members
  m_Nx_ext   = CommonParameters::Nx() + m_Nspc_size + 1;
  m_Ny_ext   = CommonParameters::Ny() + m_Nspc_size + 1;
  m_Nz_ext   = CommonParameters::Nz() + m_Nspc_size + 1;
  m_Nt_ext   = CommonParameters::Nt() + m_Ntmp_size + 1;
  m_Nvol_ext = m_Nx_ext * m_Ny_ext * m_Nz_ext * m_Nt_ext;

  //! The following setting explicitly depends on the definition
  //! of unit vectors.
  m_Nmax[0] = Nspc_size;
  m_Nmax[1] = Nspc_size;
  m_Nmax[2] = Nspc_size / 2;
  m_Nmax[3] = Nspc_size;
  m_Nmax[4] = Nspc_size / 2;
  m_Nmax[5] = Nspc_size / 2;
}


//====================================================================
void WilsonLoop::init()
{
  const int Ndim = CommonParameters::Ndim();

  assert(Ndim == 4);

  m_filename_output = "stdout";

  m_Ntype_max = 6;
  const int Ndim_spc = Ndim - 1;

  m_Nunit.resize(m_Ntype_max);
  m_Nmax.resize(m_Ntype_max);

  for (int i = 0; i < m_Ntype_max; ++i) {
    m_Nunit[i].resize(Ndim_spc);
  }

  // The following setting explicitly depends on the definition
  // of unit vectors.
  assert(m_Ntype_max >= 6);

  m_Nunit[0][0] = 1;
  m_Nunit[0][1] = 0;
  m_Nunit[0][2] = 0;

  m_Nunit[1][0] = 1;
  m_Nunit[1][1] = 1;
  m_Nunit[1][2] = 0;

  m_Nunit[2][0] = 2;
  m_Nunit[2][1] = 1;
  m_Nunit[2][2] = 0;

  m_Nunit[3][0] = 1;
  m_Nunit[3][1] = 1;
  m_Nunit[3][2] = 1;

  m_Nunit[4][0] = 2;
  m_Nunit[4][1] = 1;
  m_Nunit[4][2] = 1;

  m_Nunit[5][0] = 2;
  m_Nunit[5][1] = 2;
  m_Nunit[5][2] = 1;
}


//====================================================================
double WilsonLoop::measure(const Field_G& U)
{
  const int Ndim     = CommonParameters::Ndim();
  const int Ndim_spc = Ndim - 1;
  const int Nc       = CommonParameters::Nc();

  Mat_SU_N Uunit(Nc);

  Uunit.unit();

  //! setting extended config.
  Field_G Uext(m_Nvol_ext, Ndim);
  set_extfield(Uext, U);

  //! fixing extended config to temporal gauge.
  gfix_temporal(Uext);

  vout.paranoiac(m_vl, "%s: measurement start.\n", class_name.c_str());

  std::vector<double> wloop(m_Nspc_size * m_Ntmp_size * m_Ntype);
  for (size_t i = 0, n = wloop.size(); i < n; ++i) {
    wloop[i] = 0.0;
  }

  //! on/off-diagonal direction element_type
  for (int i_type = 0; i_type < m_Ntype; ++i_type) {
    //! permutation of direction
    for (int nu = 0; nu < Ndim_spc; ++nu) {
      std::vector<int> unit_v(Ndim_spc);
      unit_v[0] = m_Nunit[i_type][nu % Ndim_spc];
      unit_v[1] = m_Nunit[i_type][(1 + nu) % Ndim_spc];
      unit_v[2] = m_Nunit[i_type][(2 + nu) % Ndim_spc];

      int unit_v_max = unit_v[0];
      if (unit_v_max < unit_v[1]) unit_v_max = unit_v[1];
      if (unit_v_max < unit_v[2]) unit_v_max = unit_v[2];

      Field_G Uspc(m_Nvol_ext, 1);
      for (int site = 0; site < m_Nvol_ext; ++site) {
        Uspc.set_mat(site, 0, Uunit);
      }

      int Nmax = m_Nmax[i_type];
      for (int j = 0; j < Nmax; ++j) {
        // redef_Uspc(Uspc, Uext, j, unit_v);
        redef_Uspc(Uspc, Uext, j, nu, unit_v);
        //- now Uspc is product of linkv in unit_v*j direction.

        for (int t_sep = 0; t_sep < m_Ntmp_size; ++t_sep) {
          double wloop1 = calc_wloop(Uspc, t_sep + 1);
          vout.detailed(m_vl, "  %d  %d  %d  %d  %f\n",
                        i_type, nu, j + 1, t_sep + 1, wloop1);
          wloop[index_wloop(j, t_sep, i_type)] += wloop1 / 3.0;
        }
      }
    }
  }

  //! output: same format as Fortran code
  std::ostream& log_file_previous = vout.getStream();
  std::ofstream log_file;

  if (m_filename_output != "stdout") {
    log_file.open(m_filename_output.c_str(), std::ios::app);
    vout.init(log_file);
  }

  for (int i_type = 0; i_type < m_Ntype; ++i_type) {
    int Nmax = m_Nmax[i_type];
    for (int x_sep = 0; x_sep < Nmax; ++x_sep) {
      for (int t_sep = 0; t_sep < m_Ntmp_size; ++t_sep) {
        vout.general(m_vl, "  %d  %d  %d  %20.14e\n",
                     i_type + 1, x_sep + 1, t_sep + 1, wloop[index_wloop(x_sep, t_sep, i_type)]);
      }
    }
  }

  if (m_filename_output != "stdout") {
    log_file.close();
    vout.init(log_file_previous);
  }

  vout.paranoiac(m_vl, "%s: measurement finished.\n", class_name.c_str());

  //- return maximum loop with element_type=0.
  return wloop[index_wloop(m_Nmax[0] - 1, m_Ntmp_size - 1, 0)];
}


//====================================================================
double WilsonLoop::calc_wloop(const Field_G& Uspc, const int t_sep)
{
  const int Nx = CommonParameters::Nx();
  const int Ny = CommonParameters::Ny();
  const int Nz = CommonParameters::Nz();
  const int Nt = CommonParameters::Nt();
  const int Nc = CommonParameters::Nc();

  const int Nvol = CommonParameters::Nvol();
  const int NPE  = CommonParameters::NPE();

  Index_lex index_ext(m_Nx_ext, m_Ny_ext, m_Nz_ext, m_Nt_ext);

  double wloop_r = 0.0;
  double wloop_i = 0.0;

  for (int t = 0; t < Nt; ++t) {
    for (int z = 0; z < Nz; ++z) {
      for (int y = 0; y < Ny; ++y) {
        for (int x = 0; x < Nx; ++x) {
          int site1 = index_ext.site(x, y, z, t);
          int site2 = index_ext.site(x, y, z, t + t_sep);

          Mat_SU_N Utmp1(Nc);
          Utmp1 = Uspc.mat(site1, 0);

          Mat_SU_N Utmp2(Nc);
          Utmp2 = Uspc.mat_dag(site2, 0);

          Mat_SU_N Utmp3(Nc);
          Utmp3 = Utmp1 * Utmp2;

          wloop_r += ReTr(Utmp3);
          wloop_i += ImTr(Utmp3);
        }
      }
    }
  }

  wloop_r = Communicator::reduce_sum(wloop_r);
  wloop_i = Communicator::reduce_sum(wloop_i);

  wloop_r = wloop_r / Nc / Nvol / NPE;
  wloop_i = wloop_i / Nc / Nvol / NPE;

  return wloop_r;
}


//====================================================================
void WilsonLoop::redef_Uspc(Field_G& Uspc, const Field_G& Uext,
                            const int j, const int nu, const std::vector<int>& unit_v)
{
  const int Nx = CommonParameters::Nx();
  const int Ny = CommonParameters::Ny();
  const int Nz = CommonParameters::Nz();
  const int Nc = CommonParameters::Nc();

  Index_lex index_ext(m_Nx_ext, m_Ny_ext, m_Nz_ext, m_Nt_ext);

  int unit_v_max = unit_v[0];

  if (unit_v_max < unit_v[1]) unit_v_max = unit_v[1];
  if (unit_v_max < unit_v[2]) unit_v_max = unit_v[2];


  int imx = 0;
  int imy = 0;
  int imz = 0;

  //! this loop counting is for exact comparison to Fortran code.
  for (int k = 0; k < 3 * unit_v_max; ++k) {
    int kmod = (k + 3 - nu) % 3;

    if ((kmod == 0) && (imx < unit_v[0])) {
      for (int t = 0; t < m_Nt_ext; ++t) {
        for (int z = 0; z < Nz; ++z) {
          for (int y = 0; y < Ny; ++y) {
            for (int x = 0; x < Nx; ++x) {
              int x2 = x + unit_v[0] * j + imx;
              int y2 = y + unit_v[1] * j + imy;
              int z2 = z + unit_v[2] * j + imz;

              int site1 = index_ext.site(x, y, z, t);
              int site2 = index_ext.site(x2, y2, z2, t);

              Mat_SU_N Utmp1(Nc);
              Utmp1 = Uspc.mat(site1, 0);

              Mat_SU_N Utmp2(Nc);
              Utmp2 = Uext.mat(site2, 0);

              Mat_SU_N Utmp3(Nc);
              Utmp3 = Utmp1 * Utmp2;

              Uspc.set_mat(site1, 0, Utmp3);
            }
          }
        }
      }
      ++imx;
    }

    if ((kmod == 1) && (imy < unit_v[1])) {
      for (int t = 0; t < m_Nt_ext; ++t) {
        for (int z = 0; z < Nz; ++z) {
          for (int y = 0; y < Ny; ++y) {
            for (int x = 0; x < Nx; ++x) {
              int x2 = x + unit_v[0] * j + imx;
              int y2 = y + unit_v[1] * j + imy;
              int z2 = z + unit_v[2] * j + imz;

              int site1 = index_ext.site(x, y, z, t);
              int site2 = index_ext.site(x2, y2, z2, t);

              Mat_SU_N Utmp1(Nc);
              Utmp1 = Uspc.mat(site1, 0);

              Mat_SU_N Utmp2(Nc);
              Utmp2 = Uext.mat(site2, 1);

              Mat_SU_N Utmp3(Nc);
              Utmp3 = Utmp1 * Utmp2;

              Uspc.set_mat(site1, 0, Utmp3);
            }
          }
        }
      }
      ++imy;
    }

    if ((kmod == 2) && (imz < unit_v[2])) {
      for (int t = 0; t < m_Nt_ext; ++t) {
        for (int z = 0; z < Nz; ++z) {
          for (int y = 0; y < Ny; ++y) {
            for (int x = 0; x < Nx; ++x) {
              int x2 = x + unit_v[0] * j + imx;
              int y2 = y + unit_v[1] * j + imy;
              int z2 = z + unit_v[2] * j + imz;

              int site1 = index_ext.site(x, y, z, t);
              int site2 = index_ext.site(x2, y2, z2, t);

              Mat_SU_N Utmp1(Nc);
              Utmp1 = Uspc.mat(site1, 0);

              Mat_SU_N Utmp2(Nc);
              Utmp2 = Uext.mat(site2, 2);

              Mat_SU_N Utmp3(Nc);
              Utmp3 = Utmp1 * Utmp2;

              Uspc.set_mat(site1, 0, Utmp3);
            }
          }
        }
      }
      ++imz;
    }
  }
}


//====================================================================
void WilsonLoop::set_extfield(Field_G& Uext, const Field_G& Uorg)
{
  const int Ndim = CommonParameters::Ndim();
  const int Nx   = CommonParameters::Nx();
  const int Ny   = CommonParameters::Ny();
  const int Nz   = CommonParameters::Nz();
  const int Nt   = CommonParameters::Nt();
  const int NinG = Uorg.nin();

  Index_lex index_lex;
  Index_lex index_ext(m_Nx_ext, m_Ny_ext, m_Nz_ext, m_Nt_ext);

  Uext.set(0.0);

  //- bulk part of extended field same to the original field.
  for (int it = 0; it < Nt; ++it) {
    for (int iz = 0; iz < Nz; ++iz) {
      for (int iy = 0; iy < Ny; ++iy) {
        for (int ix = 0; ix < Nx; ++ix) {
          int site1 = index_lex.site(ix, iy, iz, it);
          int site2 = index_ext.site(ix, iy, iz, it);

          for (int ex = 0; ex < Ndim; ++ex) {
            Uext.set_mat(site2, ex, Uorg.mat(site1, ex));
          }
        }
      }
    }
  }

  //- setting maximum size of volume for buffer field.
  int Nmin_ext = m_Nx_ext;
  if (m_Ny_ext < Nmin_ext) Nmin_ext = m_Ny_ext;
  if (m_Nz_ext < Nmin_ext) Nmin_ext = m_Nz_ext;
  if (m_Nt_ext < Nmin_ext) Nmin_ext = m_Nt_ext;

  const int Nvol_cp = m_Nvol_ext / Nmin_ext;

  //- buffer field for copy.
  Field_G Ucp1(Nvol_cp, Ndim);
  Field_G Ucp2(Nvol_cp, Ndim);

  const int size_ex = NinG * Nvol_cp * Ndim;


  //- exchange in t-direction
  for (int it_off = 0; it_off < m_Ntmp_size + 1; ++it_off) {
    for (int iz = 0; iz < m_Nz_ext; ++iz) {
      for (int iy = 0; iy < m_Ny_ext; ++iy) {
        for (int ix = 0; ix < m_Nx_ext; ++ix) {
          int site1 = index_ext.site(ix, iy, iz, it_off);
          int site2 = ix + m_Nx_ext * (iy + m_Ny_ext * iz);

          for (int ex = 0; ex < Ndim; ++ex) {
            Ucp1.set_mat(site2, ex, Uext.mat(site1, ex));
          }
        }
      }
    }

    Communicator::exchange(size_ex, Ucp2.ptr(0), Ucp1.ptr(0), 3, 1, 0);

    for (int iz = 0; iz < m_Nz_ext; ++iz) {
      for (int iy = 0; iy < m_Ny_ext; ++iy) {
        for (int ix = 0; ix < m_Nx_ext; ++ix) {
          int site1 = ix + m_Nx_ext * (iy + m_Ny_ext * iz);
          int site2 = index_ext.site(ix, iy, iz, Nt + it_off);

          for (int ex = 0; ex < Ndim; ++ex) {
            Uext.set_mat(site2, ex, Ucp2.mat(site1, ex));
          }
        }
      }
    }
  } // end of it_off loop.

  //- exchange in z-direction
  for (int iz_off = 0; iz_off < m_Nspc_size + 1; ++iz_off) {
    for (int it = 0; it < m_Nt_ext; ++it) {
      for (int iy = 0; iy < m_Ny_ext; ++iy) {
        for (int ix = 0; ix < m_Nx_ext; ++ix) {
          int site1 = index_ext.site(ix, iy, iz_off, it);
          int site2 = ix + m_Nx_ext * (iy + m_Ny_ext * it);

          for (int ex = 0; ex < Ndim; ++ex) {
            Ucp1.set_mat(site2, ex, Uext.mat(site1, ex));
          }
        }
      }
    }

    Communicator::exchange(size_ex, Ucp2.ptr(0), Ucp1.ptr(0), 2, 1, 0);

    for (int it = 0; it < m_Nt_ext; ++it) {
      for (int iy = 0; iy < m_Ny_ext; ++iy) {
        for (int ix = 0; ix < m_Nx_ext; ++ix) {
          int site1 = ix + m_Nx_ext * (iy + m_Ny_ext * it);
          int site2 = index_ext.site(ix, iy, Nz + iz_off, it);

          for (int ex = 0; ex < Ndim; ++ex) {
            Uext.set_mat(site2, ex, Ucp2.mat(site1, ex));
          }
        }
      }
    }
  } // end of iz_off loop.

  //- exchange in y-direction
  for (int iy_off = 0; iy_off < m_Nspc_size + 1; ++iy_off) {
    for (int it = 0; it < m_Nt_ext; ++it) {
      for (int iz = 0; iz < m_Nz_ext; ++iz) {
        for (int ix = 0; ix < m_Nx_ext; ++ix) {
          int site1 = index_ext.site(ix, iy_off, iz, it);
          int site2 = ix + m_Nx_ext * (iz + m_Nz_ext * it);

          for (int ex = 0; ex < Ndim; ++ex) {
            Ucp1.set_mat(site2, ex, Uext.mat(site1, ex));
          }
        }
      }
    }

    Communicator::exchange(size_ex, Ucp2.ptr(0), Ucp1.ptr(0), 1, 1, 0);

    for (int it = 0; it < m_Nt_ext; ++it) {
      for (int iz = 0; iz < m_Nz_ext; ++iz) {
        for (int ix = 0; ix < m_Nx_ext; ++ix) {
          int site1 = ix + m_Nx_ext * (iz + m_Nz_ext * it);
          int site2 = index_ext.site(ix, Ny + iy_off, iz, it);

          for (int ex = 0; ex < Ndim; ++ex) {
            Uext.set_mat(site2, ex, Ucp2.mat(site1, ex));
          }
        }
      }
    }
  } // end of iy_off loop.

  //- exchange in x-direction
  for (int ix_off = 0; ix_off < m_Nspc_size + 1; ++ix_off) {
    for (int it = 0; it < m_Nt_ext; ++it) {
      for (int iz = 0; iz < m_Nz_ext; ++iz) {
        for (int iy = 0; iy < m_Ny_ext; ++iy) {
          int site1 = index_ext.site(ix_off, iy, iz, it);
          int site2 = iy + m_Ny_ext * (iz + m_Nz_ext * it);

          for (int ex = 0; ex < Ndim; ++ex) {
            Ucp1.set_mat(site2, ex, Uext.mat(site1, ex));
          }
        }
      }
    }

    Communicator::exchange(size_ex, Ucp2.ptr(0), Ucp1.ptr(0), 0, 1, 0);

    for (int it = 0; it < m_Nt_ext; ++it) {
      for (int iz = 0; iz < m_Nz_ext; ++iz) {
        for (int iy = 0; iy < m_Ny_ext; ++iy) {
          int site1 = iy + m_Ny_ext * (iz + m_Nz_ext * it);
          int site2 = index_ext.site(Nx + ix_off, iy, iz, it);

          for (int ex = 0; ex < Ndim; ++ex) {
            Uext.set_mat(site2, ex, Ucp2.mat(site1, ex));
          }
        }
      }
    }
  } // end of ix_off loop.
}


//====================================================================
void WilsonLoop::gfix_temporal(Field_G& Uext)
{
  const int Ndim = CommonParameters::Ndim();
  const int Nc   = CommonParameters::Nc();

  const int dir_t = Ndim - 1;

  Index_lex index_ext(m_Nx_ext, m_Ny_ext, m_Nz_ext, m_Nt_ext);

  for (int it = 1; it < m_Nt_ext; ++it) {
    for (int iz = 0; iz < m_Nz_ext; ++iz) {
      for (int iy = 0; iy < m_Ny_ext; ++iy) {
        for (int ix = 0; ix < m_Nx_ext; ++ix) {
          int site0 = index_ext.site(ix, iy, iz, it - 1);

          Mat_SU_N Utrf1(Nc);
          Utrf1 = Uext.mat(site0, dir_t);

          Mat_SU_N Utrf2(Nc);
          Utrf2 = Uext.mat_dag(site0, dir_t);

          Mat_SU_N Utmp2(Nc);
          Utmp2 = Utrf1 * Utrf2;

          Uext.set_mat(site0, 3, Utmp2);

          int site1 = index_ext.site(ix, iy, iz, it);

          for (int ex = 0; ex < Ndim; ++ex) {
            Mat_SU_N Utmp(Nc);
            Utmp  = Uext.mat(site1, ex);
            Utmp2 = Utrf1 * Utmp;
            Uext.set_mat(site1, ex, Utmp2);
          }

          if (ix > 0) {
            int site2 = index_ext.site(ix - 1, iy, iz, it);

            Mat_SU_N Utmp(Nc);
            Utmp  = Uext.mat(site2, 0);
            Utmp2 = Utmp * Utrf2;
            Uext.set_mat(site2, 0, Utmp2);
          }

          if (iy > 0) {
            int site2 = index_ext.site(ix, iy - 1, iz, it);

            Mat_SU_N Utmp(Nc);
            Utmp  = Uext.mat(site2, 1);
            Utmp2 = Utmp * Utrf2;
            Uext.set_mat(site2, 1, Utmp2);
          }

          if (iz > 0) {
            int site2 = index_ext.site(ix, iy, iz - 1, it);

            Mat_SU_N Utmp(Nc);
            Utmp  = Uext.mat(site2, 2);
            Utmp2 = Utmp * Utrf2;
            Uext.set_mat(site2, 2, Utmp2);
          }
        }
      }
    }
  }
}


//====================================================================
//============================================================END=====
