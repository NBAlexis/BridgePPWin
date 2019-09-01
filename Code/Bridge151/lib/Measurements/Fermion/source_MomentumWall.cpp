/*!
        @file    source_MomentumWall.cpp

        @brief   Momentum wall source

        @author  Noriyoshi Ishii (ishii)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 +0900 #$

        @version $LastChangedRevision: 1929 $
*/
#include "BridgeLib_Private.h"
#include "source_MomentumWall.h"

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init = Source_MomentumWall::register_factory();
}
#endif

const std::string Source_MomentumWall::class_name = "Source_MomentumWall";

//====================================================================
void Source_MomentumWall::set_parameters(const Parameters& params)
{
  const string str_vlevel = params.get_string("verbose_level");

  m_vl = vout.set_verbose_level(str_vlevel);

  //- fetch and check input parameters
  std::vector<int> source_position;
  std::vector<int> source_momentum;

  int err = 0;
  err += params.fetch_int_vector("source_position", source_position);
  err += params.fetch_int_vector("source_momentum", source_momentum);

  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameter not found.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }

  set_parameters(source_position, source_momentum);
}


//====================================================================
void Source_MomentumWall::set_parameters(const std::vector<int>& source_position,
                                         const std::vector<int>& source_momentum)
{
  // ####  parameter setup  ####
  const int Ndim  = CommonParameters::Ndim();
  const int t_dir = Ndim - 1;

  //- global lattice size
  std::vector<int> Lsize(Ndim);

  Lsize[0] = CommonParameters::Lx();
  Lsize[1] = CommonParameters::Ly();
  Lsize[2] = CommonParameters::Lz();
  Lsize[3] = CommonParameters::Lt();

  //- local size
  const int Nt = CommonParameters::Nt();

  //- print input parameters
  vout.general(m_vl, "Source for spinor field - momentum wall smeared:\n");
  for (int mu = 0; mu < Ndim; ++mu) {
    vout.general(m_vl, "  source_position[%d] = %d\n",
                 mu, source_position[mu]);
  }
  for (int mu = 0; mu < Ndim - 1; ++mu) {
    vout.general(m_vl, "  source_momentum[%d] = %d\n",
                 mu, source_momentum[mu]);
  }


  //- range check
  int err = 0;
  for (int mu = 0; mu < Ndim; ++mu) {
    // NB. Lsize[mu] > abs(source_position[mu])
    err += ParameterCheck::non_negative(Lsize[mu] - abs(source_position[mu]));
  }

  if (err) {
    vout.crucial(m_vl, "Error at %s: parameter range check failed.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }

  //- store values
  m_source_position.resize(Ndim);
  for (int mu = 0; mu < Ndim; ++mu) {
    m_source_position[mu] = (source_position[mu] + Lsize[mu]) % Lsize[mu];
  }

  m_source_momentum = source_momentum;

  //- post-process

  //- PE location in t-direction.
  const int tpe = m_source_position[t_dir] / Nt;

  m_in_node = false;

  if (tpe == Communicator::ipe(t_dir)) {
    m_in_node = true;
  }
}


//====================================================================
void Source_MomentumWall::set(Field& src, const int idx)
{
  const int Ndim = CommonParameters::Ndim();

  const int x_dir = 0;
  const int y_dir = 1;
  const int z_dir = 2;
  const int t_dir = Ndim - 1;

  //- global lattice size
  const int Lx    = CommonParameters::Lx();
  const int Ly    = CommonParameters::Ly();
  const int Lz    = CommonParameters::Lz();
  const int Lvol3 = Lx * Ly * Lz;

  //- local size
  const int Nx = CommonParameters::Nx();
  const int Ny = CommonParameters::Ny();
  const int Nz = CommonParameters::Nz();
  const int Nt = CommonParameters::Nt();

  static const double PI = 4.0 * atan(1.0);

  std::vector<double> p_unit(Ndim - 1);

  p_unit[x_dir] = (2.0 * PI / Lx) * m_source_momentum[x_dir];
  p_unit[y_dir] = (2.0 * PI / Ly) * m_source_momentum[y_dir];
  p_unit[z_dir] = (2.0 * PI / Lz) * m_source_momentum[z_dir];

  std::vector<int> ipe(Ndim - 1);
  ipe[x_dir] = Communicator::ipe(x_dir);
  ipe[y_dir] = Communicator::ipe(y_dir);
  ipe[z_dir] = Communicator::ipe(z_dir);

  //- clear field
  src.set(0.0);

  if (m_in_node) {
    const int t = m_source_position[t_dir] % Nt;

    for (int z = 0; z < Nz; ++z) {
      for (int y = 0; y < Ny; ++y) {
        for (int x = 0; x < Nx; ++x) {
          int site = m_index.site(x, y, z, t);

          int x_global = x + ipe[x_dir] * Nx;
          int y_global = y + ipe[y_dir] * Ny;
          int z_global = z + ipe[z_dir] * Nz;

          double p_x = p_unit[x_dir] * (x_global - m_source_position[x_dir]);
          double p_y = p_unit[y_dir] * (y_global - m_source_position[y_dir]);
          double p_z = p_unit[z_dir] * (z_global - m_source_position[z_dir]);

          double theta = p_x + p_y + p_z;

          int idx_r = 2 * idx;
          int idx_i = 2 * idx + 1;

          //XXX field layout: complex as two doubles
          src.set(idx_r, site, 0, cos(theta) / Lvol3);
          src.set(idx_i, site, 0, sin(theta) / Lvol3);
        }
      }
    }
  }
}


//====================================================================
void Source_MomentumWall::set(Field& src, const int i_color, const int i_spin)
{
  const int Nc  = CommonParameters::Nc();
  const int idx = i_color + Nc * i_spin;

  set(src, idx);
}


//====================================================================
void Source_MomentumWall::set_all_color(Field& src, const int i_spin)
{
  const int Nc   = CommonParameters::Nc();
  const int Ndim = CommonParameters::Ndim();

  const int x_dir = 0;
  const int y_dir = 1;
  const int z_dir = 2;
  const int t_dir = Ndim - 1;

  //- global lattice size
  const int Lx    = CommonParameters::Lx();
  const int Ly    = CommonParameters::Ly();
  const int Lz    = CommonParameters::Lz();
  const int Lvol3 = Lx * Ly * Lz;

  //- local size
  const int Nx = CommonParameters::Nx();
  const int Ny = CommonParameters::Ny();
  const int Nz = CommonParameters::Nz();
  const int Nt = CommonParameters::Nt();

  static const double PI = 4.0 * atan(1.0);

  std::vector<double> p_unit(Ndim - 1);

  p_unit[x_dir] = (2.0 * PI / Lx) * m_source_momentum[x_dir];
  p_unit[y_dir] = (2.0 * PI / Ly) * m_source_momentum[y_dir];
  p_unit[z_dir] = (2.0 * PI / Lz) * m_source_momentum[z_dir];

  std::vector<int> ipe(Ndim - 1);
  ipe[x_dir] = Communicator::ipe(x_dir);
  ipe[y_dir] = Communicator::ipe(y_dir);
  ipe[z_dir] = Communicator::ipe(z_dir);

  //- clear field
  src.set(0.0);

  if (m_in_node) {
    const int t = m_source_position[t_dir] % Nt;

    for (int z = 0; z < Nz; ++z) {
      for (int y = 0; y < Ny; ++y) {
        for (int x = 0; x < Nx; ++x) {
          int site = m_index.site(x, y, z, t);

          int x_global = x + ipe[x_dir] * Nx;
          int y_global = y + ipe[y_dir] * Ny;
          int z_global = z + ipe[z_dir] * Nz;

          double p_x = p_unit[x_dir] * (x_global - m_source_position[x_dir]);
          double p_y = p_unit[y_dir] * (y_global - m_source_position[y_dir]);
          double p_z = p_unit[z_dir] * (z_global - m_source_position[z_dir]);

          double theta = p_x + p_y + p_z;

          for (int i_color = 0; i_color < Nc; ++i_color) {
            int idx   = i_color + Nc * i_spin;
            int idx_r = 2 * idx;
            int idx_i = 2 * idx + 1;

            //XXX field layout: complex as two doubles
            src.set(idx_r, site, 0, cos(theta) / Lvol3);
            src.set(idx_i, site, 0, sin(theta) / Lvol3);
          }
        }
      }
    }
  }
}


//====================================================================
void Source_MomentumWall::set_all_color_spin(Field& src)
{
  const int Nc   = CommonParameters::Nc();
  const int Nd   = CommonParameters::Nd();
  const int Ndim = CommonParameters::Ndim();

  const int x_dir = 0;
  const int y_dir = 1;
  const int z_dir = 2;
  const int t_dir = Ndim - 1;

  //- global lattice size
  const int Lx    = CommonParameters::Lx();
  const int Ly    = CommonParameters::Ly();
  const int Lz    = CommonParameters::Lz();
  const int Lvol3 = Lx * Ly * Lz;

  //- local size
  const int Nx = CommonParameters::Nx();
  const int Ny = CommonParameters::Ny();
  const int Nz = CommonParameters::Nz();
  const int Nt = CommonParameters::Nt();

  static const double PI = 4.0 * atan(1.0);

  std::vector<double> p_unit(Ndim - 1);

  p_unit[x_dir] = (2.0 * PI / Lx) * m_source_momentum[x_dir];
  p_unit[y_dir] = (2.0 * PI / Ly) * m_source_momentum[y_dir];
  p_unit[z_dir] = (2.0 * PI / Lz) * m_source_momentum[z_dir];

  std::vector<int> ipe(Ndim - 1);
  ipe[x_dir] = Communicator::ipe(x_dir);
  ipe[y_dir] = Communicator::ipe(y_dir);
  ipe[z_dir] = Communicator::ipe(z_dir);

  //- clear field
  src.set(0.0);

  if (m_in_node) {
    const int t = m_source_position[t_dir] % Nt;

    for (int z = 0; z < Nz; ++z) {
      for (int y = 0; y < Ny; ++y) {
        for (int x = 0; x < Nx; ++x) {
          int site = m_index.site(x, y, z, t);

          int x_global = x + ipe[x_dir] * Nx;
          int y_global = y + ipe[y_dir] * Ny;
          int z_global = z + ipe[z_dir] * Nz;

          double p_x = p_unit[x_dir] * (x_global - m_source_position[x_dir]);
          double p_y = p_unit[y_dir] * (y_global - m_source_position[y_dir]);
          double p_z = p_unit[z_dir] * (z_global - m_source_position[z_dir]);

          double theta = p_x + p_y + p_z;

          for (int idx = 0; idx < Nc * Nd; ++idx) {
            int idx_r = 2 * idx;
            int idx_i = 2 * idx + 1;

            //XXX field layout: complex as two doubles
            src.set(idx_r, site, 0, cos(theta) / Lvol3);
            src.set(idx_i, site, 0, sin(theta) / Lvol3);
          }
        }
      }
    }
  }
}


//====================================================================
//============================================================END=====
