/*!
        @file    $Id: layout.h #$

        @brief

        @author  Hideo Matsufuru  (matsufuru)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 18:35:38 #$

        @version $LastChangedRevision: 1571 $
*/

#ifndef LAYOUT_INCLUDED
#define LAYOUT_INCLUDED

#include "communicator_mpi.h"
#include "Parameters/commonParameters.h"

//! Layout class for logical organisation of parallel nodes

/**
   Layout class inside MPI communicator class describes
   logical organisation of parallel nodes.
   not intended for public use.
 */

class Communicator_impl::Layout {
 public:

  static int grid_rank(int *rank, const int *gcoord); //!< find rank from grid coordinate.
  static int grid_coord(int *gcoord, const int rank); //!< find grid coordinate corresponding to rank.
  static int grid_dims(int *gdims);                   //!< find grid size.

  static int ipe(const int idir);                     //!< grid coordinate along idir direction.
  static int npe(const int idir);                     //!< grid size in idir direction.

  static int tag(int rank, int idir, int ipm);        //!< generate tag for communication.

  static int layout_setup();                          //!< initialise layout.
  static int layout_delete();                         //!< destroy layout and clean up.

  static int m_ndim;                                  //!< number of dimensions.
  static int *m_dims;                                 //!< lattice extent (Lx, Ly, Lz, Lt)

  // logical layout
  static int *m_grid_dims;  //!< grid dimensions in directions.
  static int *m_grid_coord; //!< grid coordinate.

  static int *m_ipe_up;     //!< rank of upward neighbour in directions.
  static int *m_ipe_dn;     //!< rank of downward neighbour in directions.

  // physical mapping
  static char m_map_grid[16];
  static int  *m_physical_to_logical; //!< map between physical and logical grid

  static int physical_map_setup();
  static int physical_map_delete();

  // subdimensional-slices
  static MPI_Comm *m_sub_comm;  //!< subgrid

  static int subgrid_setup();
  static int subgrid_delete();

 private:
  Layout() {}                       //!< constructor is hidden as private; no instantiation.
  Layout(const Layout&) {} //!< copy constructor is hidden as private; no instantiation.
  Layout& operator=(const Layout&); //!< assignment is hidden as private.
};

/*
  struct LogicalLayout {
    int ndim;
    int *dims; // global lattice size

    int *grid_dims;  // number of division along axis
    int *grid_coord;  // logical grid coord

    int *ipe_up;  // rank of upward neighbour along ith axis
    int *ipe_dn;  // rank of downward neighbour along ith axis

    MPI_Comm *sub_comm;  // communicator of reduced dimensions
  };

  static LogicalLayout m_layout;

  static int layout_setup ();
  static int layout_delete ();

  static int subgrid_setup ();
  static int subgrid_delete ();

  static int physical_map_setup ();
  static int physical_map_delete ();
*/
#endif /* m__COMMUNICATOR_MPI_LAYOUT_H */
