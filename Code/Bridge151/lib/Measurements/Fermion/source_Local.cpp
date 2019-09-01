/*!
        @file    source_Local.cpp

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/
#include "BridgeLib_Private.h"
#include "source_Local.h"

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init = Source_Local::register_factory();
}
#endif

const std::string Source_Local::class_name = "Source_Local";

//====================================================================
void Source_Local::set_parameters(const Parameters& params)
{
  const string str_vlevel = params.get_string("verbose_level");

  m_vl = vout.set_verbose_level(str_vlevel);

  //- fetch and check input parameters
  std::vector<int> source_position;

  int err = 0;
  err += params.fetch_int_vector("source_position", source_position);

  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameter not found.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }

  set_parameters(source_position);
}


//====================================================================
void Source_Local::set_parameters(const std::vector<int>& source_position)
{
  // ####  parameter setup  ####
  int Ndim = CommonParameters::Ndim();

  std::vector<int> Lsize(Ndim);

  Lsize[0] = CommonParameters::Lx();
  Lsize[1] = CommonParameters::Ly();
  Lsize[2] = CommonParameters::Lz();
  Lsize[3] = CommonParameters::Lt();

  //- print input parameters
  vout.general(m_vl, "Source for spinor field - local:\n");
  for (int mu = 0; mu < Ndim; ++mu) {
    vout.general(m_vl, "  source_position[%d] = %d\n",
                 mu, source_position[mu]);
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

  assert(source_position.size() == Ndim);

  //- store values
  m_source_position.resize(Ndim);
  m_source_position = source_position;

  for (int mu = 0; mu < Ndim; ++mu) {
    while (m_source_position[mu] < 0)
    {
      m_source_position[mu] += Lsize[mu];
    }
    m_source_position[mu] %= Lsize[mu];
  }


  //- post-process
  std::vector<int> Nsize(Ndim);
  Nsize[0] = CommonParameters::Nx();
  Nsize[1] = CommonParameters::Ny();
  Nsize[2] = CommonParameters::Nz();
  Nsize[3] = CommonParameters::Nt();

  //- check if source position is in current local node.
  m_in_node = true;
  for (int i = 0; i < Ndim; ++i) {
    int inode = m_source_position[i] / Nsize[i];

    if (inode != Communicator::ipe(i)) {
      m_in_node = false;
      break;
    }
  }

  //- find local coordinate of source position.
  if (m_in_node) {
    for (int i = 0; i < Ndim; ++i) {
      m_source_position[i] %= Nsize[i];
    }
  } else {
    for (int i = 0; i < Ndim; ++i) {
      m_source_position[i] = -1;
    }
  }
}


//====================================================================
void Source_Local::set(Field& src, const int idx)
{
  //- clear field
  src.set(0.0);  // src = 0.0

  if (m_in_node) {
    const int isite = m_index.site(m_source_position[0],
                                   m_source_position[1],
                                   m_source_position[2],
                                   m_source_position[3]);

    //XXX complex as two doubles
    src.set(2 * idx, isite, 0, 1.0);
  }
}


//====================================================================
void Source_Local::set(Field& src, const int i_color, const int i_spin)
{
  const int Nc  = CommonParameters::Nc();
  const int idx = i_color + Nc * i_spin;

  set(src, idx);
}


//====================================================================
void Source_Local::set_all_color(Field& src, const int i_spin)
{
  const int Nc = CommonParameters::Nc();

  //- clear field
  src.set(0.0);

  if (m_in_node) {
    const int isite = m_index.site(m_source_position[0],
                                   m_source_position[1],
                                   m_source_position[2],
                                   m_source_position[3]);

    for (int i_color = 0; i_color < Nc; ++i_color) {
      int idx   = i_color + Nc * i_spin;
      int idx_r = 2 * idx;

      //XXX complex as two doubles
      src.set(idx_r, isite, 0, 1.0);
    }
  }
}


//====================================================================
void Source_Local::set_all_color_spin(Field& src)
{
  const int Nc = CommonParameters::Nc();
  const int Nd = CommonParameters::Nd();

  //- clear field
  src.set(0.0);

  if (m_in_node) {
    const int isite = m_index.site(m_source_position[0],
                                   m_source_position[1],
                                   m_source_position[2],
                                   m_source_position[3]);

    for (int idx = 0; idx < Nc * Nd; ++idx) {
      int idx_r = 2 * idx;

      //XXX complex as two doubles
      src.set(idx_r, isite, 0, 1.0);
    }
  }
}


//====================================================================
//============================================================END=====
