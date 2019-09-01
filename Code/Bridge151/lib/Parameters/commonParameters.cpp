/*!
        @file    commonParameters.cpp

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/
#include "BridgeLib_Private.h"
//#include "commonParameters.h"
//
//#include "Communicator/communicator.h"
//
//#include "IO/bridgeIO.h"

using Bridge::vout;

// initialize static parameters
// color, spinor and space-time dimension
int       CommonParameters::m_Nc   = 3; // `const' removed. [07 May 2014]
const int CommonParameters::m_Nd   = 4;
const int CommonParameters::m_Ndim = 4;

const double CommonParameters::m_epsilon_criterion = 1.0e-16;

int CommonParameters::m_Lx = 0;
int CommonParameters::m_Ly = 0;
int CommonParameters::m_Lz = 0;
int CommonParameters::m_Lt = 0;

long_t CommonParameters::m_Lvol = 0;

int CommonParameters::m_NPEx = 0;
int CommonParameters::m_NPEy = 0;
int CommonParameters::m_NPEz = 0;
int CommonParameters::m_NPEt = 0;
int CommonParameters::m_NPE  = 0;

// note: obsolete. moved to MPI/Communicator_impl.
char CommonParameters::m_map_grid[] = "xyzt";

int CommonParameters::m_Nx   = 0;
int CommonParameters::m_Ny   = 0;
int CommonParameters::m_Nz   = 0;
int CommonParameters::m_Nt   = 0;
int CommonParameters::m_Nvol = 0;

//CommonParameters* CommonParameters::m_instance = 0;
bool CommonParameters::m_initialized = false;

Bridge::VerboseLevel CommonParameters::m_vlevel = Bridge::GENERAL;

const std::string CommonParameters::class_name = "CommonParameters";

//====================================================================
void CommonParameters::init(const std::vector<int>& lattice_size,
                            const std::vector<int>& grid_size)
{
  init(lattice_size, grid_size, 3);
}


//====================================================================
void CommonParameters::init(const std::vector<int>& lattice_size,
                            const std::vector<int>& grid_size,
                            const int Nc)
{
  if (m_initialized == true) {
    vout.crucial("Error at %s: CommonParameters already initialized.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }

  m_Nc = Nc;

  m_Lx = lattice_size[0];
  m_Ly = lattice_size[1];
  m_Lz = lattice_size[2];
  m_Lt = lattice_size[3];

  m_NPEx = grid_size[0];
  m_NPEy = grid_size[1];
  m_NPEz = grid_size[2];
  m_NPEt = grid_size[3];

  if (check_parameters() == false) {
    vout.crucial("Error at %s: check_parameters failed.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }

  print_parameters();

  m_initialized = true;
}


//====================================================================
void CommonParameters::init_Vlevel(Bridge::VerboseLevel vlevel)
{
  m_vlevel = vlevel;
}


//====================================================================
int CommonParameters::Lsize(const int dir)
{
  switch (dir)
  {
  case 0:
    return m_Lx;

  case 1:
    return m_Ly;

  case 2:
    return m_Lz;

  case 3:
    return m_Lt;

  default:
    vout.crucial("Error at %s: invalid argument value.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
    //return -1;
  }
}


//====================================================================
int CommonParameters::Nsize(const int dir)
{
  switch (dir)
  {
  case 0:
    return m_Nx;

  case 1:
    return m_Ny;

  case 2:
    return m_Nz;

  case 3:
    return m_Nt;

  default:
    vout.crucial("Error at %s: invalid argument value.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
    //return -1;
  }
}


//====================================================================
int CommonParameters::NPEsize(const int dir)
{
  switch (dir)
  {
  case 0:
    return m_NPEx;

  case 1:
    return m_NPEy;

  case 2:
    return m_NPEz;

  case 3:
    return m_NPEt;

  default:
    vout.crucial("Error at %s: invalid argument value.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
    //return -1;
  }
}


//====================================================================
bool CommonParameters::check_parameters()
{
  //if (m_Ndim == 0) return false;

  if (m_Nc == 0) return false;

  //if (m_Nd == 0) return false;

  if (m_Lx == 0) return false;

  if (m_Ly == 0) return false;

  if (m_Lz == 0) return false;

  if (m_Lt == 0) return false;

  m_Lvol = m_Lx * m_Ly * m_Lz * m_Lt;

  if (m_NPE == 1) {
    // if (m_NPEx != 1 || m_NPEy != 1 || m_NPEz != 1 || m_NPEt != 1) {
    //   vout.general(m_vl,stderr, "assume NPE_mu = 1\n");
    // }
    m_NPEx = m_NPEy = m_NPEz = m_NPEt = 1;
  }

  if ((m_Nx != 0) && (m_NPEx != 0) && (m_Lx != m_Nx * m_NPEx)) return false;

  if ((m_Ny != 0) && (m_NPEy != 0) && (m_Ly != m_Ny * m_NPEy)) return false;

  if ((m_Nz != 0) && (m_NPEz != 0) && (m_Lz != m_Nz * m_NPEz)) return false;

  if ((m_Nt != 0) && (m_NPEt != 0) && (m_Lt != m_Nt * m_NPEt)) return false;

  if ((m_Nx == 0) && (m_NPEx != 0)) {
    if (m_Lx % m_NPEx != 0) return false; else m_Nx = m_Lx / m_NPEx;
  }
  if ((m_Ny == 0) && (m_NPEy != 0)) {
    if (m_Ly % m_NPEy != 0) return false; else m_Ny = m_Ly / m_NPEy;
  }
  if ((m_Nz == 0) && (m_NPEz != 0)) {
    if (m_Lz % m_NPEz != 0) return false; else m_Nz = m_Lz / m_NPEz;
  }
  if ((m_Nt == 0) && (m_NPEt != 0)) {
    if (m_Lt % m_NPEt != 0) return false; else m_Nt = m_Lt / m_NPEt;
  }

  if ((m_Nx != 0) && (m_NPEx == 0)) {
    if (m_Lx % m_Nx != 0) return false; else m_NPEx = m_Lx / m_Nx;
  }
  if ((m_Ny != 0) && (m_NPEy == 0)) {
    if (m_Ly % m_Ny != 0) return false; else m_NPEy = m_Ly / m_Ny;
  }
  if ((m_Nz != 0) && (m_NPEz == 0)) {
    if (m_Lz % m_Nz != 0) return false; else m_NPEz = m_Lz / m_Nz;
  }
  if ((m_Nt != 0) && (m_NPEt == 0)) {
    if (m_Lt % m_Nt != 0) return false; else m_NPEt = m_Lt / m_Nt;
  }

  if ((m_NPEx != 0) && (m_NPEy != 0) && (m_NPEz != 0) && (m_NPEt != 0)) {
    if ((m_NPE != 0) && (m_NPE != m_NPEx * m_NPEy * m_NPEz * m_NPEt)) return false;

    if (m_NPE == 0) m_NPE = m_NPEx * m_NPEy * m_NPEz * m_NPEt;
  }

  if ((m_Nx != 0) && (m_Ny != 0) && (m_Nz != 0) && (m_Nt != 0)) {
    if ((m_Nvol != 0) && (m_Nvol != m_Nx * m_Ny * m_Nz * m_Nt)) return false;

    if (m_Nvol == 0) m_Nvol = m_Nx * m_Ny * m_Nz * m_Nt;
  }

  return true;
}


//====================================================================
void CommonParameters::print_parameters()
{
  vout.general("\n");
  vout.general("Lattice parameters:\n");
  vout.general("  Lx = %d,  Ly = %d,  Lz = %d,  Lt = %d,\n",
               m_Lx, m_Ly, m_Lz, m_Lt);
  vout.general("  Nx = %d,  Ny = %d,  Nz = %d,  Nt = %d,\n",
               m_Nx, m_Ny, m_Nz, m_Nt);
  vout.general("  NPEx = %d,  NPEy  = %d,  NPEz  = %d,  NPEt  = %d,\n",
               m_NPEx, m_NPEy, m_NPEz, m_NPEt);
  vout.general("  Lvol = %ld,  Nvol = %d,  NPE = %d,\n",
               m_Lvol, m_Nvol, m_NPE);
  vout.general("  Ndim = %d\n", m_Ndim);
  vout.general("  Nc   = %d\n", m_Nc);
  vout.general("  Nd   = %d\n", m_Nd);
  vout.general("\n");
}


//====================================================================
//============================================================END=====
