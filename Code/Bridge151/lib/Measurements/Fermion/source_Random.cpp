/*!
        @file    source_Random.cpp

        @brief

        @author  Yusuke Namekawa (namekawa)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/
#include "BridgeLib_Private.h"
#include "source_Random.h"

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init = Source_Random::register_factory();
}
#endif

const std::string Source_Random::class_name = "Source_Random";

//====================================================================
void Source_Random::set_parameters(const Parameters& params)
{
  const string str_vlevel = params.get_string("verbose_level");

  m_vl = vout.set_verbose_level(str_vlevel);

  //- fetch and check input parameters
  std::vector<int> source_position;
  std::vector<int> source_momentum;
  string           str_noise_type;

  int err = 0;
  err += params.fetch_int_vector("source_position", source_position);
  err += params.fetch_int_vector("source_momentum", source_momentum);
  err += params.fetch_string("noise_type", str_noise_type);

  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameter not found.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }


  set_parameters(source_position, source_momentum, str_noise_type);
}


//====================================================================
void Source_Random::set_parameters(const std::vector<int>& source_position,
                                   const std::vector<int>& source_momentum,
                                   const string str_noise_type)
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
  vout.general(m_vl, "Source for spinor field - Random:\n");
  for (int mu = 0; mu < Ndim; ++mu) {
    vout.general(m_vl, "  source_position[%d] = %d\n",
                 mu, source_position[mu]);
  }

  for (int mu = 0; mu < Ndim - 1; ++mu) {
    vout.general(m_vl, "  source_momentum[%d] = %d\n",
                 mu, source_momentum[mu]);
  }

  vout.general(m_vl, "  noise_type          = %s\n", str_noise_type.c_str());


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
  m_str_noise_type  = str_noise_type;


  //-- post-process
  // m_rand->reset(static_cast<unsigned long>(i_seed_noise));

  //- PE location in t-direction.
  const int tpe = m_source_position[t_dir] / Nt;

  m_in_node = false;

  if (tpe == Communicator::ipe(t_dir)) {
    m_in_node = true;
  }
}


//====================================================================
void Source_Random::set(Field& src, const int idx)
{
  const int Ndim = CommonParameters::Ndim();

  //- local size
  std::vector<int> Nsize(Ndim);

  Nsize[0] = CommonParameters::Nx();
  Nsize[1] = CommonParameters::Ny();
  Nsize[2] = CommonParameters::Nz();
  Nsize[3] = CommonParameters::Nt();

  Field v(src);

  if (m_str_noise_type == "Gaussian") {
    m_rand->gauss_lex_global(v);
  } else if (m_str_noise_type == "U1") {
    m_rand->U1_lex_global(v);
  } else if (m_str_noise_type == "Z2") {
    m_rand->Z2_lex_global(v);
  } else {
    vout.crucial(m_vl, "Error at %s: Noise type of %s is not supported.\n", class_name.c_str(), m_str_noise_type.c_str());
    exit(EXIT_FAILURE);
  }

  //- clear field
  src.set(0.0);

  const int idx_r = 2 * idx;
  const int idx_i = 2 * idx + 1;

  if (m_in_node) {
    int t = m_source_position[3] % Nsize[3];

    for (int z = 0; z < Nsize[2]; ++z) {
      for (int y = 0; y < Nsize[1]; ++y) {
        for (int x = 0; x < Nsize[0]; ++x) {
          int site = m_index.site(x, y, z, t);

          double rn_r = v.cmp(idx_r, site, 0);
          double rn_i = v.cmp(idx_i, site, 0);

          //XXX field layout: complex as two doubles
          src.set(idx_r, site, 0, rn_r);
          src.set(idx_i, site, 0, rn_i);

          //- for debug
          vout.paranoiac(m_vl, "%s: src(x,y,z) = %e %e (%d %d %d)\n",
                         class_name.c_str(), rn_r, rn_i, x, y, z);
        }
      }
    }
  }

#ifdef DEBUG
  //- global lattice size
  const int Lx    = CommonParameters::Lx();
  const int Ly    = CommonParameters::Ly();
  const int Lz    = CommonParameters::Lz();
  const int Lvol3 = Lx * Ly * Lz;

  const int Nvol = CommonParameters::Nvol();

  double src1_r = 0.0;
  double src1_i = 0.0;
  for (int site = 0; site < Nvol; ++site) {
    src1_r += src.cmp(idx_r, site, 0);
    src1_i += src.cmp(idx_i, site, 0);
  }
  const double src_r = Communicator::reduce_sum(src1_r) / Lvol3;
  const double src_i = Communicator::reduce_sum(src1_i) / Lvol3;
  vout.general(m_vl, "%s: sum_x src = %e %e\n", class_name.c_str(), src_r, src_i);

  const double src2 = src.norm2() / Lvol3;
  vout.general(m_vl, "%s: src.norm2 = %e\n", class_name.c_str(), src2);
#endif
}


//====================================================================
void Source_Random::set(Field& src, const int i_color, const int i_spin)
{
  const int Nc  = CommonParameters::Nc();
  const int idx = i_color + Nc * i_spin;

  set(src, idx);
}


//====================================================================
void Source_Random::set_all_color(Field& src, const int i_spin)
{
  const int Nc   = CommonParameters::Nc();
  const int Ndim = CommonParameters::Ndim();

  //- local size
  std::vector<int> Nsize(Ndim);

  Nsize[0] = CommonParameters::Nx();
  Nsize[1] = CommonParameters::Ny();
  Nsize[2] = CommonParameters::Nz();
  Nsize[3] = CommonParameters::Nt();

  Field v(src);

  if (m_str_noise_type == "Gaussian") {
    m_rand->gauss_lex_global(v);
  } else if (m_str_noise_type == "U1") {
    m_rand->U1_lex_global(v);
  } else if (m_str_noise_type == "Z2") {
    m_rand->Z2_lex_global(v);
  } else {
    vout.crucial(m_vl, "Error at %s: Noise type of %s is not supported.\n", class_name.c_str(), m_str_noise_type.c_str());
    exit(EXIT_FAILURE);
  }

  //- clear field
  src.set(0.0);

  if (m_in_node) {
    int t = m_source_position[3] % Nsize[3];

    for (int z = 0; z < Nsize[2]; ++z) {
      for (int y = 0; y < Nsize[1]; ++y) {
        for (int x = 0; x < Nsize[0]; ++x) {
          int site = m_index.site(x, y, z, t);

          for (int i_color = 0; i_color < Nc; ++i_color) {
            int idx   = i_color + Nc * i_spin;
            int idx_r = 2 * idx;
            int idx_i = 2 * idx + 1;

            double rn_r = v.cmp(idx_r, site, 0);
            double rn_i = v.cmp(idx_i, site, 0);

            //XXX field layout: complex as two doubles
            src.set(idx_r, site, 0, rn_r);
            src.set(idx_i, site, 0, rn_i);

            //- for debug
            vout.paranoiac(m_vl, "%s: src(idx, x,y,z) = %e %e (%d %d %d %d)\n",
                           class_name.c_str(), rn_r, rn_i, idx, x, y, z);
          }
        }
      }
    }
  }

#ifdef DEBUG
  //- global lattice size
  const int Lx    = CommonParameters::Lx();
  const int Ly    = CommonParameters::Ly();
  const int Lz    = CommonParameters::Lz();
  const int Lvol3 = Lx * Ly * Lz;

  const int Nvol = CommonParameters::Nvol();

  const int Nd = CommonParameters::Nd();

  double src1_r = 0.0;
  double src1_i = 0.0;
  for (int site = 0; site < Nvol; ++site) {
    for (int idx = 0; idx < Nc * Nd; ++idx) {
      int idx_r = 2 * idx;
      int idx_i = 2 * idx + 1;

      src1_r += src.cmp(idx_r, site, 0);
      src1_i += src.cmp(idx_i, site, 0);
    }
  }
  const double src_r = Communicator::reduce_sum(src1_r) / (Lvol3 * Nc * Nd);
  const double src_i = Communicator::reduce_sum(src1_i) / (Lvol3 * Nc * Nd);
  vout.general(m_vl, "%s: sum_x src = %e %e\n", class_name.c_str(), src_r, src_i);

  const double src2 = src.norm2() / (Lvol3 * Nc * Nd);
  vout.general(m_vl, "%s: src.norm2 = %e\n", class_name.c_str(), src2);
#endif
}


//====================================================================
void Source_Random::set_all_color_spin(Field& src)
{
  const int Nc   = CommonParameters::Nc();
  const int Nd   = CommonParameters::Nd();
  const int Ndim = CommonParameters::Ndim();

  //- local size
  std::vector<int> Nsize(Ndim);

  Nsize[0] = CommonParameters::Nx();
  Nsize[1] = CommonParameters::Ny();
  Nsize[2] = CommonParameters::Nz();
  Nsize[3] = CommonParameters::Nt();

  Field v(src);

  if (m_str_noise_type == "Gaussian") {
    m_rand->gauss_lex_global(v);
  } else if (m_str_noise_type == "U1") {
    m_rand->U1_lex_global(v);
  } else if (m_str_noise_type == "Z2") {
    m_rand->Z2_lex_global(v);
  } else {
    vout.crucial(m_vl, "Error at %s: Noise type of %s is not supported.\n", class_name.c_str(), m_str_noise_type.c_str());
    exit(EXIT_FAILURE);
  }

  //- clear field
  src.set(0.0);

  if (m_in_node) {
    int t = m_source_position[3] % Nsize[3];

    for (int z = 0; z < Nsize[2]; ++z) {
      for (int y = 0; y < Nsize[1]; ++y) {
        for (int x = 0; x < Nsize[0]; ++x) {
          int site = m_index.site(x, y, z, t);

          for (int idx = 0; idx < Nc * Nd; ++idx) {
            int idx_r = 2 * idx;
            int idx_i = 2 * idx + 1;

            double rn_r = v.cmp(idx_r, site, 0);
            double rn_i = v.cmp(idx_i, site, 0);

            //XXX field layout: complex as two doubles
            src.set(idx_r, site, 0, rn_r);
            src.set(idx_i, site, 0, rn_i);

            //- for debug
            vout.paranoiac(m_vl, "%s: src(idx, x,y,z) = %e %e (%d %d %d %d)\n",
                           class_name.c_str(), rn_r, rn_i, idx, x, y, z);
          }
        }
      }
    }
  }

#ifdef DEBUG
  //- global lattice size
  const int Lx    = CommonParameters::Lx();
  const int Ly    = CommonParameters::Ly();
  const int Lz    = CommonParameters::Lz();
  const int Lvol3 = Lx * Ly * Lz;

  const int Nvol = CommonParameters::Nvol();

  double src1_r = 0.0;
  double src1_i = 0.0;
  for (int site = 0; site < Nvol; ++site) {
    for (int idx = 0; idx < Nc * Nd; ++idx) {
      int idx_r = 2 * idx;
      int idx_i = 2 * idx + 1;

      src1_r += src.cmp(idx_r, site, 0);
      src1_i += src.cmp(idx_i, site, 0);
    }
  }
  const double src_r = Communicator::reduce_sum(src1_r) / (Lvol3 * Nc * Nd);
  const double src_i = Communicator::reduce_sum(src1_i) / (Lvol3 * Nc * Nd);
  vout.general(m_vl, "%s: sum_x src = %e %e\n", class_name.c_str(), src_r, src_i);

  const double src2 = src.norm2() / (Lvol3 * Nc * Nd);
  vout.general(m_vl, "%s: src.norm2 = %e\n", class_name.c_str(), src2);
#endif
}


//====================================================================
void Source_Random::set_all_space_time(Field& src, const int ic)
{
  // This implementation assumes the given field src is complex field.
  int nc   = CommonParameters::Nc();
  int Nex  = src.nex();
  int Nvol = src.nvol();
  int Nin  = src.nin();

  assert(ic < nc);
  assert(Nin >= 2 * nc);
  int nd = Nin / nc / 2;
//  double rn, rn2, rn3;
  double rn;
//  double RF2 = 1.0 / sqrt(2.0);

  Field v(2, Nvol, Nex);
  if (m_str_noise_type == "Gaussian") {
    m_rand->gauss_lex_global(v);
  } else if (m_str_noise_type == "U1") {
    m_rand->U1_lex_global(v);
  } else if (m_str_noise_type == "Z2") {
    m_rand->Z2_lex_global(v);
  } else {
    vout.crucial(m_vl, "Error at %s: Noise type of %s is not supported.\n", class_name.c_str(), m_str_noise_type.c_str());
    exit(EXIT_FAILURE);
  }
  //  m_rand->uniform_lex_global(v);
  src.set(0.0);

  for (int ex = 0; ex < Nex; ++ex) {
    for (int site = 0; site < Nvol; ++site) {
      for (int ri = 0; ri < 2; ++ri) {
        rn = v.cmp(ri, site, ex);
        for (int id = 0; id < nd; ++id) {
          int in = ri + 2 * (ic + nc * id);
//          rn2 = floor(2.0 * rn);
//          rn3 = (2.0 * rn2 - 1.0) * RF2;
//    src.set(in, site, ex, rn3);
          src.set(in, site, ex, rn);
        }
      }
    }
  }
}


//====================================================================
void Source_Random::set_all_space_time(Field& src, const int ic, const int is)
{
  // This implementation assumes the given field src is complex field.
  int nc   = CommonParameters::Nc();
#if DEBUG
  int nd   = CommonParameters::Nd();
#endif
  int Nex  = src.nex();
  int Nvol = src.nvol();
#if DEBUG
  int Nin  = src.nin();
#endif

  assert(ic < nc);
  assert(is < nd);
#if DEBUG
  assert(Nin == 2 * nc * nd);
#endif

  Field v(2, Nvol, Nex);
  if (m_str_noise_type == "Gaussian") {
    m_rand->gauss_lex_global(v);
  } else if (m_str_noise_type == "U1") {
    m_rand->U1_lex_global(v);
  } else if (m_str_noise_type == "Z2") {
    m_rand->Z2_lex_global(v);
  } else {
    vout.crucial(m_vl, "Error at %s: Noise type of %s is not supported.\n", class_name.c_str(), m_str_noise_type.c_str());
    exit(EXIT_FAILURE);
  }
  //  m_rand->uniform_lex_global(v);
  src.set(0.0);

//  double rn, rn2, rn3;
  double rn;
//  double RF2 = 1.0 / sqrt(2.0);

  for (int ex = 0; ex < Nex; ++ex) {
    for (int site = 0; site < Nvol; ++site) {
      for (int ri = 0; ri < 2; ++ri) {
        int in = ri + 2 * (ic + nc * is);
        rn = v.cmp(ri, site, ex);
//        rn2 = floor(2.0 * rn);
//        rn3 = (2.0 * rn2 - 1.0) * RF2;
//        src.set(in, site, ex, rn3);
        src.set(in, site, ex, rn);
      }
    }
  }
}


//====================================================================
//============================================================END=====
