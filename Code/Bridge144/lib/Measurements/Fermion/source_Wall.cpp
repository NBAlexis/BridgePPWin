#include "BridgeLib_Private.h"

/*!
        @file    $Id:: source_Wall.cpp #$

        @brief

        @author  Yusuke Namekawa (namekawa)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 00:49:46 #$

        @version $LastChangedRevision: 1561 $
*/

#include "source_Wall.h"



#ifdef USE_FACTORY
namespace {
  Source *create_object()
  {
    return new Source_Wall();
  }


  bool init = Source::Factory::Register("Wall", create_object);
}
#endif



const std::string Source_Wall::class_name = "Source_Wall";

//====================================================================
void Source_Wall::set_parameters(const Parameters& params)
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
void Source_Wall::set_parameters(const std::vector<int>& source_position)
{
  // ####  parameter setup  ####
  int Ndim = CommonParameters::Ndim();

  //- global lattice size
  std::vector<int> Lsize(Ndim);

  Lsize[0] = CommonParameters::Lx();
  Lsize[1] = CommonParameters::Ly();
  Lsize[2] = CommonParameters::Lz();
  Lsize[3] = CommonParameters::Lt();

  //- local size
  std::vector<int> Nsize(Ndim);
  Nsize[0] = CommonParameters::Nx();
  Nsize[1] = CommonParameters::Ny();
  Nsize[2] = CommonParameters::Nz();
  Nsize[3] = CommonParameters::Nt();

  const int t_dir = Ndim - 1;

  //- print input parameters
  vout.general(m_vl, "Source for spinor field - Wall smeared:\n");
  vout.general(m_vl, "  source_position[t] = %d\n", source_position[t_dir]);

  //- range check
  int err = 0;
  // NB. Lsize[t_dir] > abs(source_position[t_dir])
  err += ParameterCheck::non_negative(Lsize[t_dir] - abs(source_position[t_dir]));

  if (err) {
    vout.crucial(m_vl, "Error at %s: parameter range check failed.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }

  //- store values
  m_source_position.resize(Ndim);
  for (int mu = 0; mu < Ndim; ++mu) {
    m_source_position[mu] = 0;
  }
  m_source_position[t_dir] = (source_position[t_dir] + Lsize[t_dir]) % Lsize[t_dir];


  //- post-process

  //- PE location in t-direction.
  int tpe = m_source_position[t_dir] / Nsize[t_dir];

  m_in_node = false;

  if (tpe == Communicator::ipe(t_dir)) {
    m_in_node = true;
  }
}


//====================================================================
void Source_Wall::set(Field& src, int j)
{
  int Ndim = CommonParameters::Ndim();

  //- global lattice size
  int Lx = CommonParameters::Lx();
  int Ly = CommonParameters::Ly();
  int Lz = CommonParameters::Lz();

  const int Lvol3 = Lx * Ly * Lz;

  //- local size
  std::vector<int> Nsize(Ndim);

  Nsize[0] = CommonParameters::Nx();
  Nsize[1] = CommonParameters::Ny();
  Nsize[2] = CommonParameters::Nz();
  Nsize[3] = CommonParameters::Nt();

  //- clear field
  src.set(0.0);

  if (m_in_node) {
    int t = m_source_position[3] % Nsize[3];

    for (int z = 0; z < Nsize[2]; ++z) {
      for (int y = 0; y < Nsize[1]; ++y) {
        for (int x = 0; x < Nsize[0]; ++x) {
          int isite = m_index.site(x, y, z, t);

          //XXX field layout: complex as two doubles
          src.set(2 * j, isite, 0, 1.0 / Lvol3);
        }
      }
    }
  }
}


//====================================================================
//============================================================END=====
