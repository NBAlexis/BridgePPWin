#include "BridgeLib_Private.h"

/*!
        @file    $Id:: randomNumbers.cpp #$

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 00:49:46 #$

        @version $LastChangedRevision: 1561 $
*/


#include "randomNumbers.h"

#include <fstream>
#include <cassert>

#include "Field/index_lex.h"
#include "Field/index_eo.h"

namespace {
  const double sq2r = 1.0 / sqrt(2.0);
  // const double pi   = 3.141592653589793;
  const double pi  = 4.0 * atan(1.0);
  const double pi2 = 2.0 * pi;
}

const std::string RandomNumbers::class_name = "RandomNumbers";

//====================================================================
void RandomNumbers::gauss(double& rand1, double& rand2)
{
  //   Two Gaussian random number with deviation 1/\sqrt(2).
  double r1 = get();
  double r2 = get();

  double slg1 = sqrt(-log(1.0 - r1) * 2.0) * sq2r;
  double ang1 = pi2 * r2;

  rand1 = slg1 * cos(ang1);
  rand2 = slg1 * sin(ang1);
}


//====================================================================
void RandomNumbers::gauss_lex_global(Field& f)
{
  if (f.nin() % 2 == 0) {
    return generate_global<rand_gauss_even>(f);
  } else {
    return generate_global<rand_gauss_odd>(f);
  }
}


//====================================================================
void RandomNumbers::uniform_lex_global(Field& f)
{
  return generate_global<rand_uniform>(f);
}


//====================================================================
void RandomNumbers::gauss_eo_global(Field& f)
{
  Field g = f.clone();

  gauss_lex_global(g);

  Index_eo idx;
  idx.convertField(f, g);
}


//====================================================================
void RandomNumbers::rand_gauss_even::operator()(const bool do_fill)
{
  if (do_fill) {
    for (size_t i = 0; i < m_block; i += 2) {
      double r1 = m_rand_gauss->get();
      double r2 = m_rand_gauss->get();

      double slg1 = sqrt(-log(1 - r1) * 2.0) * sq2r;
      double ang1 = pi2 * r2;

      *m_ptr++ = slg1 * cos(ang1);
      *m_ptr++ = slg1 * sin(ang1);
    }
  } else {
    // just let rand_gauss state proceed
    for (size_t i = 0; i < m_block; i += 2) {
      m_rand_gauss->get();
      m_rand_gauss->get();
    }
  }
}


size_t RandomNumbers::rand_gauss_even::block_size() const
{
  return m_block;
}


//====================================================================
void RandomNumbers::rand_gauss_odd::operator()(const bool do_fill)
{
  if (do_fill) {
    size_t ngauss = m_block / 2;

    for (size_t i = 0; i < ngauss; ++i) {
      double r1 = m_rand_gauss->get();
      double r2 = m_rand_gauss->get();

      double slg1 = sqrt(-log(1 - r1) * 2.0) * sq2r;
      double ang1 = pi2 * r2;

      *m_ptr++ = slg1 * cos(ang1);
      *m_ptr++ = slg1 * sin(ang1);
    }

    // remaining one
    {
      double r1 = m_rand_gauss->get();
      double r2 = m_rand_gauss->get();

      double slg1 = sqrt(-log(1 - r1) * 2.0) * sq2r;
      double ang1 = pi2 * r2;

      *m_ptr++ = slg1 * cos(ang1);
      // the other one is dropped.
    }
  } else {
    for (size_t i = 0; i < m_block + 1; ++i) {
      m_rand_gauss->get();
    }
  }
}


size_t RandomNumbers::rand_gauss_odd::block_size() const
{
  return m_block + 1;
}


//====================================================================
void RandomNumbers::rand_uniform::operator()(const bool do_fill)
{
  if (do_fill) {
    for (size_t i = 0; i < m_block; ++i) {
      *m_ptr++ = m_rand_gauss->get();
    }
  } else {
    for (size_t i = 0; i < m_block; ++i) {
      m_rand_gauss->get();
    }
  }
}


size_t RandomNumbers::rand_uniform::block_size() const
{
  return m_block;
}


//====================================================================
template<typename InnerGenerator>
void RandomNumbers::generate_global(Field& f)
{
  InnerGenerator fill(f, this);

  //int Nin  = f.nin();
  int Nvol = f.nvol();
  int Nex  = f.nex();

  if (Communicator::size() == 1) {
    vout.detailed(m_vl, "%s: single node. no need to consider division.\n", class_name.c_str());

    for (int j = 0; j < Nex; ++j) {
      for (int i = 0; i < Nvol; ++i) {
        fill(true);
      }
    }
  } else {
    assert(Nvol == CommonParameters::Nvol());

    int Lx = CommonParameters::Lx();
    int Ly = CommonParameters::Ly();
    int Lz = CommonParameters::Lz();
    int Lt = CommonParameters::Lt();

    int Nx = CommonParameters::Nx();
    int Ny = CommonParameters::Ny();
    int Nz = CommonParameters::Nz();
    int Nt = CommonParameters::Nt();

    int gx = Communicator::ipe(0);
    int gy = Communicator::ipe(1);
    int gz = Communicator::ipe(2);
    int gt = Communicator::ipe(3);

    for (int j = 0; j < Nex; ++j) {
      bool in_j = true;

      for (int t = 0; t < Lt; ++t) {
        bool in_t = in_j && (t >= gt * Nt) && (t < (gt + 1) * Nt);

        for (int z = 0; z < Lz; ++z) {
          bool in_z = in_t && (z >= gz * Nz) && (z < (gz + 1) * Nz);

          for (int y = 0; y < Ly; ++y) {
            bool in_y = in_z && (y >= gy * Ny) && (y < (gy + 1) * Ny);

            for (int x = 0; x < Lx; ++x) {
              bool in_x = in_y && (x >= gx * Nx) && (x < (gx + 1) * Nx);

              fill(in_x);
            }
          }
        }
      }
    }
  } // end if communicator::size == 1
}


//====================================================================
//============================================================END=====
