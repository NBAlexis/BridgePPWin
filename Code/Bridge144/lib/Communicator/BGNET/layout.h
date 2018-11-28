/*!
        @file    $Id: layout.h #$

        @brief

        @author  Hideo Matsufuru  (matsufuru)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 18:35:38 #$

        @version $LastChangedRevision: 1571 $
*/

#ifndef COMMUNICATOR_MPI_LAYOUT_INCLUDED
#define COMMUNICATOR_MPI_LAYOUT_INCLUDED

#include "communicator_bgnet.h"

class Communicator_impl::Layout {
 public:

  static int grid_rank(int *rank, const int *gcoord);
  static int grid_coord(int *gcoord, const int rank);
  static int grid_dims(int *gdims);

  static int ipe(const int idir);
  static int npe(const int idir);

  static int tag(int rank, int idir, int ipm);

  static int layout_setup();
  static int layout_delete();

  static int m_ndim;
  static int *m_dims;  //< lattice extent

  // logical layout
  static int *m_grid_dims;
  static int *m_grid_coord;

  static int *m_ipe_up;
  static int *m_ipe_dn;

  // physical mapping
  static char m_map_grid[16];
  static int  *m_physical_to_logical; //< map between physical and logical grid

  static int physical_map_setup();
  static int physical_map_delete();

  // subdimensional-slices
  //  static MPI_Comm *m_sub_comm;
  static int *m_sub_comm;

  static int subgrid_setup();
  static int subgrid_delete();

 private:
  Layout() {}
  Layout(const Layout&) {}
  Layout& operator=(const Layout&);
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
#endif /* COMMUNICATOR_MPI_LAYOUT_INCLUDED */
