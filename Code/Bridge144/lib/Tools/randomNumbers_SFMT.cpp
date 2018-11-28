#include "BridgeLib_Private.h"

/*!
        @file    $Id:: randomNumbers_SFMT.cpp #$

        @brief

        @author  T.Aoyama (aoym)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 18:35:38 #$

        @version $LastChangedRevision: 1571 $
*/

#ifdef USE_SFMTLIB

#include "randomNumbers_SFMT.h"

#include "timer.h"

#include <iomanip>

const std::string RandomNumbers_SFMT::class_name = "RandomNumbers:SFMT";

#ifdef ENABLE_SFMT_JUMP
static const NTL::GF2X ch_poly = sfmt_characteristic_polynomial();
#endif

//====================================================================
RandomNumbers_SFMT::RandomNumbers_SFMT(const int s)
{
  // initialize with seed s
  sfmt_init_gen_rand(&m_state, s);
}


//====================================================================
void RandomNumbers_SFMT::reset(unsigned long seed)
{
  // initialize with seed
  sfmt_init_gen_rand(&m_state, seed);
}


//====================================================================
double RandomNumbers_SFMT::get()
{
  return sfmt_genrand_res53(&m_state);
}


//====================================================================
void RandomNumbers_SFMT::get_block(double *v, const size_t n)
{
  for (size_t i = 0; i < n; ++i) {
    v[i] = get();
  }
}


//====================================================================
#ifdef ENABLE_SFMT_JUMP

void RandomNumbers_SFMT::uniform_lex_global(Field& f)
{
  return generate_global_jump<RandomNumbers::rand_uniform>(f);
}


//====================================================================
void RandomNumbers_SFMT::gauss_lex_global(Field& f)
{
  if (f.nin() % 2 == 0) {
    return generate_global_jump<RandomNumbers::rand_gauss_even>(f);
  } else {
    return generate_global_jump<RandomNumbers::rand_gauss_odd>(f);
  }
}


//====================================================================
template<typename InnerGenerator>
void RandomNumbers_SFMT::generate_global_jump(Field& f)
{
  InnerGenerator fill(f, this);

  const int Nin  = f.nin();
  const int Nvol = f.nvol();
  const int Nex  = f.nex();

  if (Communicator::size() == 1) {
    vout.detailed(m_vl, "%s: single node. no need to consider division.\n", class_name.c_str());

    for (int j = 0; j < Nex; ++j) {
      for (int isite = 0; isite < Nvol; ++isite) {
        fill(true);
      }
    }
  } else {
    if (f.nvol() != CommonParameters::Nvol()) {
      vout.crucial(m_vl, "Error at %s::%s: Nvol mismatch.\n", class_name.c_str(), __func__);
      exit(EXIT_FAILURE);
    }

    // duplicate mt state for backup.
    sfmt_t rng_state = m_state;

    const int Lx   = CommonParameters::Lx();
    const int Ly   = CommonParameters::Ly();
    const int Lz   = CommonParameters::Lz();
    const int Lt   = CommonParameters::Lt();
    const int Lvol = CommonParameters::Lvol();

    const int Nx   = CommonParameters::Nx();
    const int Ny   = CommonParameters::Ny();
    const int Nz   = CommonParameters::Nz();
    const int Nt   = CommonParameters::Nt();
    const int Nvol = CommonParameters::Nvol();

    const int Gx  = Communicator::npe(0);
    const int Gy  = Communicator::npe(1);
    const int Gz  = Communicator::npe(2);
    const int Gt  = Communicator::npe(3);
    const int NPE = Communicator::size();

    vout.detailed(m_vl, "L = { %2d, %2d, %2d, %2d }, Lvol = %d\n", Lx, Ly, Lz, Lt, Lvol);
    vout.detailed(m_vl, "N = { %2d, %2d, %2d, %2d }, Nvol = %d\n", Nx, Ny, Nz, Nt, Nvol);
    vout.detailed(m_vl, "G = { %2d, %2d, %2d, %2d }, NPE  = %d\n", Gx, Gy, Gz, Gt, NPE);

    // calculate jump polynomials
    sfmt_jump_t jump_y;
    sfmt_jump_t jump_z;
    sfmt_jump_t jump_t;
    // for external degree of freedom
    sfmt_jump_t jump_w;

    size_t block_size = fill.block_size();

    if (Lx - Nx != 0)
      sfmt_calculate_jump_polynomial(jump_y, block_size * (Lx - Nx), ch_poly);
    if (Ly - Ny != 0)
      sfmt_calculate_jump_polynomial(jump_z, block_size * (Lx * (Ly - Ny)), ch_poly);
    if (Lz - Nz != 0)
      sfmt_calculate_jump_polynomial(jump_t, block_size * (Lx * Ly * (Lz - Nz)), ch_poly);
    if (Lt - Nt != 0)
      sfmt_calculate_jump_polynomial(jump_w, block_size * (Lx * Ly * Lz * (Lt - Nt)), ch_poly);

#define calc_global_index(x, y, z, t) \
  (x) + Lx * ((y) + Ly * ((z) + Lz * (t)))

#define calc_local_index(x, y, z, t) \
  (x) + Nx * ((y) + Ny * ((z) + Nz * (t)))

    int ipe = Communicator::self();

    int gx = Communicator::ipe(0);
    int gy = Communicator::ipe(1);
    int gz = Communicator::ipe(2);
    int gt = Communicator::ipe(3);

//    vout.detailed(m_vl, "ipe = %d: (%d, %d, %d, %d)\n", ipe, gx, gy, gz, gt);

    int global_index = calc_global_index(Nx * gx, Ny * gy, Nz * gz, Nt * gt);

    // initial shift
    sfmt_jump(&m_state, block_size * global_index, ch_poly);

    // generate random numbers
    double *p = f.ptr(0);

    for (int j = 0; j < Nex; ++j) {
      for (int t = 0; t < Nt; ++t) {
        for (int z = 0; z < Nz; ++z) {
          for (int y = 0; y < Ny; ++y) {
            for (int x = 0; x < Nx; ++x) {
              fill(true);
            }

            // jump to next row
            if (Lx - Nx != 0)
              sfmt_jump_by_polynomial(&m_state, jump_y);
          }

          if (Ly - Ny != 0)
            sfmt_jump_by_polynomial(&m_state, jump_z);
        }

        if (Lz - Nz != 0)
          sfmt_jump_by_polynomial(&m_state, jump_t);
      }

      if (Lt - Nt != 0)
        sfmt_jump_by_polynomial(&m_state, jump_w);
    }

    // restore original state and jump by global size.
    m_state = rng_state;
    sfmt_jump(&m_state, block_size * Lvol * Nex, ch_poly);


#undef calc_global_index
#undef calc_local_index
  }  // if communicator::size == 1
}
#endif /* ENABLE_SFMT_JUMP */

//====================================================================

// /**
//  * SFMT internal state
//  */
// struct SFMT_T {
//     /** the 128-bit internal state array */
//     w128_t state[SFMT_N];
//     /** index counter to the 32-bit internal state array */
//     int idx;
// };
//
// typedef struct SFMT_T sfmt_t;

void RandomNumbers_SFMT::write_file(const std::string& filename)
{
  vout.detailed(m_vl, "%s: save random number state to file %s\n", class_name.c_str(), filename.c_str());

  if (Communicator::is_primary()) {
    std::ofstream out_file(filename.c_str());

    if (!out_file) {
      vout.crucial(m_vl, "Error at %s: unable to open output file.\n", class_name.c_str());
      exit(EXIT_FAILURE);
    }

    for (int i = 0; i < SFMT_N; ++i) {
      out_file << std::setw(8) << std::setfill('0') << std::hex << m_state.state[i].u[0] << " ";
      out_file << std::setw(8) << std::setfill('0') << std::hex << m_state.state[i].u[1] << " ";
      out_file << std::setw(8) << std::setfill('0') << std::hex << m_state.state[i].u[2] << " ";
      out_file << std::setw(8) << std::setfill('0') << std::hex << m_state.state[i].u[3] << " ";
      out_file << std::endl;
    }

    out_file << std::dec << m_state.idx << std::endl;

    out_file.close();
  }
}


//====================================================================
void RandomNumbers_SFMT::read_file(const std::string& filename)
{
  vout.detailed(m_vl, "%s: read random number state from file %s\n", class_name.c_str(), filename.c_str());

  if (Communicator::is_primary()) {
    std::ifstream in_file(filename.c_str());

    if (!in_file) {
      vout.crucial(m_vl, "Error at %s: unable to open output file.\n", class_name.c_str());
      exit(EXIT_FAILURE);
    }

    for (int i = 0; i < SFMT_N; ++i) {
      in_file >> std::hex >> m_state.state[i].u[0];
      in_file >> std::hex >> m_state.state[i].u[1];
      in_file >> std::hex >> m_state.state[i].u[2];
      in_file >> std::hex >> m_state.state[i].u[3];
    }

    in_file >> std::dec >> m_state.idx;

    in_file.close();
  }

  Communicator::Base::broadcast(sizeof(sfmt_t), &m_state, 0);
}


//====================================================================
//============================================================END=====
#endif  // USE_SFMTLIB
