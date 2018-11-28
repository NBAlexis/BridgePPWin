#include "BridgeLib_Private.h"

#if COMMUNICATOR_USE_BGNET

/*!
        @file    $Id: physical_map.cpp #$

        @brief

        @author  Hideo Matsufuru  (matsufuru)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 00:49:46 #$

        @version $LastChangedRevision: 1561 $
*/

#include "layout.h"

#include <cstdarg>
#include <cstring>

#include "Parameters/commonParameters.h"

static int find_coord_map(const char map_expr[], int ndim, int *map);

// physical rank is assigned to physical grid coord in column major order.
//   phys_coord[0] x phys_coord[1] x ... (fastest to slowest)
// other assignment may be chosen.

//====================================================================
// logical grid coord to physical rank number
int Communicator_impl::Layout::grid_rank(int *rank, const int *gcoord)
{
  int r = 0;

  for (int i = m_ndim - 1; i >= 0; --i) {
    int k = m_physical_to_logical[i];
    r = r * m_grid_dims[k] + gcoord[k];
  }

  *rank = r;

  return EXIT_SUCCESS;
}


//====================================================================
//* physical rank to logical grid coord
int Communicator_impl::Layout::grid_coord(int *gcoord, const int rank)
{
  int r = rank;

  for (int i = 0; i < m_ndim; ++i) {
    int k = m_physical_to_logical[i];

    gcoord[k] = r % m_grid_dims[k];
    r        /= m_grid_dims[k];
  }

  return EXIT_SUCCESS;
}


//====================================================================
int Communicator_impl::Layout::physical_map_setup()
{
  strncpy(m_map_grid, CommonParameters::Grid_map(), sizeof(m_map_grid) / sizeof(char));

  // physical to logical map
  m_physical_to_logical = new int [m_ndim];

  find_coord_map(m_map_grid, m_ndim, m_physical_to_logical);

  return EXIT_SUCCESS;
}


//====================================================================
int Communicator_impl::Layout::physical_map_delete()
{
  delete [] m_physical_to_logical;

  return EXIT_SUCCESS;
}


//====================================================================
//* translate {x, y, z, t, (w)} into {0, 1, 2, 3, (4)}
static int find_coord_map(const char map_expr[], int ndim, int *map)
{
  int len = strlen(map_expr);

  if (len != ndim) {
    fprintf(stderr, "ERROR: find_map: length mismatch: %d, expected %d.\n", len, ndim);
    return EXIT_FAILURE;
  }

  for (int i = 0; i < len; ++i) {
    char c   = map_expr[i];
    int  idx = -1;

    switch (c)
    {
    case 'x':
      idx = 0;
      break;

    case 'y':
      idx = 1;
      break;

    case 'z':
      idx = 2;
      break;

    case 't':
      idx = 3;
      break;

    case 'w':
      idx = 4;
      break;                  // extra 5th dimension

    default:
      fprintf(stderr, "ERROR: find_map: unknown symbol %c.\n", c);
      return EXIT_FAILURE;
    }

    map[i] = idx;
  }

  return EXIT_SUCCESS;
}


//====================================================================
//============================================================END=====

#endif //COMMUNICATOR_USE_BGNET
