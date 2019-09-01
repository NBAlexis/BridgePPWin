/*!
        @file    layout.cpp

        @brief

        @author  Hideo Matsufuru  (matsufuru)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/
#include "BridgeLib_Private.h"
//#include "layout.h"

// prototype declaration
namespace {
  static int pe_logical_layout(const int ndim, const int *dims, int nproc, int *npe);
  //static int find_primes(const int n, int *p);
}

// static members
int Communicator_impl::Layout::m_ndim  = 0;
int *Communicator_impl::Layout::m_dims = 0;  //< lattice extent

// logical layout
int *Communicator_impl::Layout::m_grid_dims  = 0;
int *Communicator_impl::Layout::m_grid_coord = 0;

int *Communicator_impl::Layout::m_ipe_up = 0;
int *Communicator_impl::Layout::m_ipe_dn = 0;

// physical mapping
char Communicator_impl::Layout::m_map_grid[16];
int  *Communicator_impl::Layout::m_physical_to_logical = 0; //< map between physical and logical grid

// subdimensional-slices
MPI_Comm *Communicator_impl::Layout::m_sub_comm = 0;


//====================================================================
// communicator_impl methods
int Communicator_impl::Layout::ipe(const int idir)
{
  return m_grid_coord[idir];
}


//====================================================================
int Communicator_impl::Layout::npe(const int idir)
{
  return m_grid_dims[idir];
}


//====================================================================
int Communicator_impl::Layout::grid_dims(int *dims)
{
  for (int i = 0; i < m_ndim; ++i) {
    dims[i] = m_grid_dims[i];
  }
  return EXIT_SUCCESS;
}


//====================================================================
// tag() : identify transfer channel by source rank and direction (incl fw/bw)
int Communicator_impl::Layout::tag(int rank, int idir, int ipm)
{
  return 2 * m_ndim * rank + idir + ((ipm > 0) ? 0 : m_ndim);
}


//====================================================================
//! layout_setup()  -- setup logical layout

/**
   layout_setup()

   find dimension, global lattice size, local lattice size, grid division
   from CommonParameters.
   if grid division or local size is specified, check consistency.
   otherwise, try to divide so that global lattice size is divisable and
   the total number of rank is equal to number of processor elements.

   N.B. the result should be stored back to communicator for further use
   from other components.

   N.B. initial parameter settings are taken from CommonParameters; should
   rather be fetched from parameter file or something else?

*/
int Communicator_impl::Layout::layout_setup()
{
  // fetch parameter values from CommonParameters

  m_ndim = CommonParameters::Ndim();

  // lattice size
  m_dims = new int [m_ndim];

  m_dims[0] = CommonParameters::Lx();
  m_dims[1] = CommonParameters::Ly();
  m_dims[2] = CommonParameters::Lz();
  m_dims[3] = CommonParameters::Lt();

  // suggested logical layout
  m_grid_dims = new int [m_ndim];

  m_grid_dims[0] = CommonParameters::NPEx();
  m_grid_dims[1] = CommonParameters::NPEy();
  m_grid_dims[2] = CommonParameters::NPEz();
  m_grid_dims[3] = CommonParameters::NPEt();

  // suggested local volume
  int *local_dims = new int [m_ndim];

  local_dims[0] = CommonParameters::Nx();
  local_dims[1] = CommonParameters::Ny();
  local_dims[2] = CommonParameters::Nz();
  local_dims[3] = CommonParameters::Nt();

  for (int i = 0; i < m_ndim; ++i) {
    if (local_dims[i] != 0) {  // there is a suggestion of local extent.
      if (m_grid_dims[i] != 0) {
        if (m_grid_dims[i] * local_dims[i] != m_dims[i]) {
          fprintf(stderr, "layout mismatch.\n");
          exit(EXIT_FAILURE);
        } else {
          // ok. m_grid_dims[i] is as specified.
        }
      } else {  // use specified local size and find number of divisions.
        if (m_dims[i] % local_dims[i] != 0) {
          fprintf(stderr, "layout mismatch. lattice undivisable by local volume.\n");
          exit(EXIT_FAILURE);
        } else {
          m_grid_dims[i] = m_dims[i] / local_dims[i];
        }
      }
    }
  }

  delete [] local_dims;


  std::vector<int> lattice_size(m_dims, m_dims + 4);
  std::vector<int> grid_size(m_grid_dims, m_grid_dims + 4);

  // continue setup
  return layout_setup(lattice_size, grid_size);
}


int Communicator_impl::Layout::layout_setup(
  const std::vector<int>& lattice_size,
  std::vector<int>& grid_size)
{
  m_ndim = static_cast<int>(lattice_size.size());

  // lattice size
  m_dims = new int [m_ndim];

  for (int i = 0; i < m_ndim; ++i) {
    m_dims[i] = lattice_size[i];
  }

  // suggested logical layout
  m_grid_dims = new int [m_ndim];

  for (int i = 0; i < m_ndim; ++i) {
    m_grid_dims[i] = 0;
  }

  for (int i = 0; i < grid_size.size(); ++i) {
    m_grid_dims[i] = grid_size[i];
  }


  // find logical layout
  int retv = pe_logical_layout(m_ndim, m_dims, m_grid_size, m_grid_dims);

  // note: Communicator_impl::m_grid_size = size of local communicator m_comm.

  if (retv != EXIT_SUCCESS) {
    fprintf(stderr, "layout failed.\n");
    exit(EXIT_FAILURE);
  }


  // store back to grid_size
  grid_size.resize(m_ndim);
  for (int i = 0; i < m_ndim; ++i) {
    grid_size[i] = m_grid_dims[i];
  }


  // physical to logical map
  physical_map_setup();

  // find logical coordinate of my rank
  m_grid_coord = new int [m_ndim];
  grid_coord(m_grid_coord, m_grid_rank);

  // find neighbour rank
  m_ipe_up = new int [m_ndim];
  m_ipe_dn = new int [m_ndim];

  int *coord = new int [m_ndim];
  for (int i = 0; i < m_ndim; ++i) {
    for (int j = 0; j < m_ndim; ++j) {
      coord[j] = m_grid_coord[j];
    }

    // upward
    coord[i] = (m_grid_coord[i] + 1) % m_grid_dims[i];
    grid_rank(&m_ipe_up[i], coord);

    // downward
    coord[i] = (m_grid_coord[i] - 1 + m_grid_dims[i]) % m_grid_dims[i];
    grid_rank(&m_ipe_dn[i], coord);
  }
  delete [] coord;


#ifdef DEBUG
  printf("rank %d: up=(%d,%d,%d,%d), dn=(%d,%d,%d,%d)\n",
         Communicator::self(),
         m_ipe_up[0], m_ipe_up[1], m_ipe_up[2], m_ipe_up[3],
         m_ipe_dn[0], m_ipe_dn[1], m_ipe_dn[2], m_ipe_dn[3]);
#endif

  // subgrid of reduced dimension.
  subgrid_setup();

  return EXIT_SUCCESS;
}


//====================================================================
int Communicator_impl::Layout::layout_delete()
{
  subgrid_delete();

  physical_map_delete();

  delete [] m_dims;
  delete [] m_grid_dims;
  delete [] m_grid_coord;
  delete [] m_ipe_up;
  delete [] m_ipe_dn;

  return EXIT_SUCCESS;
}


//====================================================================
// subgrid for communication in reduced dimensional subsets
int Communicator_impl::Layout::subgrid_setup()
{
  //  int Nmask = (1 << m_ndim);
  int Nmask = 1;

  // temporary modified to avoid error on BG/Q. [23 May 2014 H.M.]

  m_sub_comm = new MPI_Comm [Nmask];

  for (int imask = 0; imask < Nmask; ++imask) 
  {
      int* coord = (int*)alloca(sizeof(int) * m_ndim);
    //int coord[m_ndim];

    for (int i = 0; i < m_ndim; ++i) 
    {
      coord[i] = m_grid_coord[i];
    }

    for (int i = 0; i < m_ndim; ++i) {
      bool mask = ((imask >> i) & 1) == 1 ? true : false;
      if (!mask) coord[i] = 0;
    }

    int rank = 0;
    grid_rank(&rank, coord);

/*
  printf("split %2d: rank %2d: (%d,%d,%d,%d) -> (%d,%d,%d,%d), %2d\n",
  imask,
  m_grid_rank,
  m_grid_coord[0], m_grid_coord[1], m_grid_coord[2], m_grid_coord[3],
  coord[0], coord[1], coord[2], coord[3],
  rank);
*/
    MPI_Comm_split(m_comm, rank, 0 /* key */, &m_sub_comm[imask]);
  }

  return EXIT_SUCCESS;
}


//====================================================================
int Communicator_impl::Layout::subgrid_delete()
{
  delete [] m_sub_comm;

  return EXIT_SUCCESS;
}


//====================================================================
namespace { // anonymous namespace
// layout
  static const int prime_table[] =
  {
    2,     3,   5,   7,  11,  13,  17,  19,
    23,   29,  31,  37,  41,  43,  47,  53,
    59,   61,  67,  71,  73,  79,  83,  89,
    97,  101, 103, 107, 109, 113, 127, 131,
    137, 139, 149, 151, 157, 163, 167, 173,
    179, 181, 191, 193, 197, 199, 211, 223,
    227, 229, 233, 239, 241, 251, 257, 263,
    269, 271, 277, 281, 283, 293, 307, 311,
  };

  static const int nprime = sizeof(prime_table) / sizeof(int);

//* find logical layout
//*   divide dims_mu into npe_mu x local_mu where prod npe_mu = nproc.
//*   npe_mu != 0 is kept intact.
  static int pe_logical_layout(const int ndim, const int *dims, int nproc, int *npe)
  {
    int retv = EXIT_SUCCESS;

    int nfreedim  = 0;
    int nfreeproc = nproc;

    for (int i = 0; i < ndim; ++i) {
      if (npe[i] == 0) {
        ++nfreedim;
      } else {
        if (npe[i] < 0) {
          fprintf(stderr, "illegal value: npe[%d]=%d.\n", i, npe[i]);
          return EXIT_FAILURE;
        } else if (nproc % npe[i] != 0) {
          fprintf(stderr, "illegal value: npe[%d]=%d does not divide NPE=%d.\n", i, npe[i], nproc);
          return EXIT_FAILURE;
        } else if (nfreeproc % npe[i] != 0) {
          fprintf(stderr, "illegal value: NPE=%d is not divisable by %d.\n", nproc, nproc / nfreeproc * npe[i]);
          return EXIT_FAILURE;
        } else if (dims[i] % npe[i] != 0) {
          fprintf(stderr, "illegal value: npe[%d]=%d does not divide L=%d.\n", i, npe[i], dims[i]);
          return EXIT_FAILURE;
        } else {
          nfreeproc /= npe[i];
        }
      }
    }

    if (nfreeproc < 1) {
      fprintf(stderr, "impossible layout.\n");
      return EXIT_FAILURE;
    } else if (nfreeproc == 1) {
      for (int i = 0; i < ndim; ++i) {
        if (npe[i] == 0) npe[i] = 1;
      }
      return EXIT_SUCCESS;  // complete.
    } else {
      if (nfreedim == 0) {
        fprintf(stderr, "impossible layout. no room to divide.\n");
        return EXIT_FAILURE;
      }
    }

    // divide nfreeproc to nfreedim npe's.
    int np = nfreeproc;

    int *subdims = new int [ndim];
    int nf       = 0;
    for (int i = 0; i < ndim; ++i) {
      if (npe[i] == 0) subdims[nf++] = dims[i]; }

    int *count = new int [nprime];
    for (int i = 0; i < nprime; ++i) {
      count[i] = 0;
    }

    for (int i = 0; i < nprime; ++i) {
      int p = prime_table[i];
      while (np > 1 && np % p == 0)
      {
        ++count[i];
        np /= p;
      }
      if (np == 1) break;
    }

//   printf("%d = (", nfreeproc);
//   for (int i=0; i<nprime; ++i) {
//     printf("%d, ", count[i]);
//   }
//   printf(")\n");

    if (np > 1) {
      fprintf(stderr, "insufficient prime table.\n");
      retv = EXIT_FAILURE;
      goto finalize;
    }

//   printf("subdims=(");
//   for (int i=0; i<nf; ++i) {
//     printf("%d, ", subdims[i]);
//   }
//   printf(")\n");

    for (int i = nprime - 1; i >= 0; --i) {
      if (count[i] == 0) continue;

      int p = prime_table[i];

      for (int j = 0; j < count[i]; ++j) {
        int maxsubdim = 0;
        int maxk      = -1;
        for (int k = 0; k < nf; ++k) {
          if ((subdims[k] >= maxsubdim) && (subdims[k] % p == 0)) {
            maxsubdim = subdims[k];
            maxk      = k;
          }
        }
//      printf("prime=%d, maxk=%d, maxsubdim=%d\n", p, maxk, maxsubdim);

        if (maxk == -1) {
          fprintf(stderr, "not divisable. %d\n", p);
          retv = EXIT_FAILURE;
          goto finalize;
        }

        subdims[maxk] /= p;
      }
    }

/*
  printf("subdims=(");
  for (int i=0; i<nf; ++i) {
  printf("%d, ", subdims[i]);
  }
  printf(")\n");
*/
    // store
    for (int i = 0, k = 0; i < ndim; ++i) {
      if (npe[i] == 0) {
        npe[i] = dims[i] / subdims[k];
        ++k;
      }
    }

finalize:
    delete [] subdims;
    delete [] count;

    return retv;
  }


  //* find n primes.
  //static int find_primes(const int n, int *p)
  //{
  //  if (n < 1) return EXIT_FAILURE;

  //  int i = 0;
  //  int k = 2;

  //  p[i++] = k++;

  //  while (1)
  //  {
  //    int j;
  //    for (j = 0; j < i; ++j) {
  //      if (k % p[j] == 0) break;
  //    }

  //    if (j >= i) {  // not divisable by p[0]..p[i-1]. found new prime.
  //      p[i++] = k;
  //      if (i >= n) return EXIT_FAILURE;
  //    }

  //    ++k;
  //  }

  //  return EXIT_SUCCESS;
  //}
} // anonymous namespace

//====================================================================
//============================================================END=====
