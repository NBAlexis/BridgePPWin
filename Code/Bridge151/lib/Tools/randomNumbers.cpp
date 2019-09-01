/*!
        @file    randomNumbers.cpp

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/
#include "BridgeLib_Private.h"
#include "randomNumbers.h"

#include <fstream>
#include <cassert>

#include "Field/index_eo.h"

#ifdef USE_FACTORY

#ifdef USE_FACTORY_AUTOREGISTER
#else
// setup factories for all subclasses

#include "randomNumbers_Mseries.h"
#include "randomNumbers_MT19937.h"
#ifdef USE_SFMTLIB
#include "randomNumbers_SFMT.h"
#endif

bool RandomNumbers::init_factory()
{
  bool result = true;

  result &= RandomNumbers_Mseries::register_factory();
  result &= RandomNumbers_MT19937::register_factory();
#ifdef USE_SFMTLIB
  result &= RandomNumbers_SFMT::register_factory();
#endif

  return result;
}


#endif /* USE_FACTORY_AUTOREGISTER */

#endif /* USE_FACTORY */


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
void RandomNumbers::lex_global(const std::string& str_rand_type, Field& f)
{
  if (str_rand_type == "Gaussian") {
    gauss_lex_global(f);
  } else if (str_rand_type == "Uniform") {
    uniform_lex_global(f);
  } else if (str_rand_type == "U1") {
    U1_lex_global(f);
  } else if (str_rand_type == "Z2") {
    Z2_lex_global(f);
  } else {
    vout.crucial(m_vl, "Error at %s: unsupported rand_type \"%s\"\n", class_name.c_str(), str_rand_type.c_str());
    exit(EXIT_FAILURE);
  }
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
void RandomNumbers::U1_lex_global(Field& f)
{
  // This implementation assumes the given field f is complex field.
  gauss_lex_global(f);

  const int Nex  = f.nex();
  const int Nvol = f.nvol();
  const int Nin  = f.nin();

  assert(Nin % 2 == 0);

  for (int ex = 0; ex < Nex; ++ex) {
    for (int site = 0; site < Nvol; ++site) {
      for (int in_half = 0; in_half < Nin / 2; ++in_half) {
        int idx_r = 2 * in_half;
        int idx_i = 2 * in_half + 1;

        //- NB. arg = atan2(y,x), not atan2(x,y)
        double arg = atan2(f.cmp(idx_i, site, ex), f.cmp(idx_r, site, ex));

        f.set(idx_r, site, ex, cos(arg));
        f.set(idx_i, site, ex, sin(arg));
      }
    }
  }
}


//====================================================================
void RandomNumbers::Z2_lex_global(Field& f)
{
  // This implementation assumes the given field f is complex field.
  generate_global<rand_uniform>(f);

  const int Nex  = f.nex();
  const int Nvol = f.nvol();
  const int Nin  = f.nin();

  assert(Nin % 2 == 0);

  double RF2 = 1.0 / sqrt(2.0);

  for (int ex = 0; ex < Nex; ++ex) {
    for (int site = 0; site < Nvol; ++site) {
      for (int in = 0; in < Nin; ++in) {
        double rn1 = f.cmp(in, site, ex);
        double rn2 = floor(2.0 * rn1);
        double rn3 = (2.0 * rn2 - 1.0) * RF2;
        f.set(in, site, ex, rn3);
      }
    }
  }
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

  const int Nin  = f.nin();
  const int Nvol = f.nvol();
  const int Nex  = f.nex();

  if (Communicator::size() == 1) {
    vout.detailed(m_vl, "%s: single node. No need to consider division.\n", class_name.c_str());

    for (int j = 0; j < Nex; ++j) {
      for (int i = 0; i < Nvol; ++i) {
        fill(true);
      }
    }
  } else {
    assert(Nvol == CommonParameters::Nvol());

    const int Lx = CommonParameters::Lx();
    const int Ly = CommonParameters::Ly();
    const int Lz = CommonParameters::Lz();
    const int Lt = CommonParameters::Lt();

    const int Nx = CommonParameters::Nx();
    const int Ny = CommonParameters::Ny();
    const int Nz = CommonParameters::Nz();
    const int Nt = CommonParameters::Nt();

    const int ipe_x = Communicator::ipe(0);
    const int ipe_y = Communicator::ipe(1);
    const int ipe_z = Communicator::ipe(2);
    const int ipe_t = Communicator::ipe(3);

    for (int j = 0; j < Nex; ++j) {
      bool in_j = true;

      for (int t = 0; t < Lt; ++t) {
        bool in_t = in_j && (t >= ipe_t * Nt) && (t < (ipe_t + 1) * Nt);

        for (int z = 0; z < Lz; ++z) {
          bool in_z = in_t && (z >= ipe_z * Nz) && (z < (ipe_z + 1) * Nz);

          for (int y = 0; y < Ly; ++y) {
            bool in_y = in_z && (y >= ipe_y * Ny) && (y < (ipe_y + 1) * Ny);

            for (int x = 0; x < Lx; ++x) {
              bool in_x = in_y && (x >= ipe_x * Nx) && (x < (ipe_x + 1) * Nx);

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
