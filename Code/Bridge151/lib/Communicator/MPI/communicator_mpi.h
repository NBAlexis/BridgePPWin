/*!
        @file    communicator_mpi.h

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/

#ifndef COMMUNICATOR_MPI_INCLUDED
#define COMMUNICATOR_MPI_INCLUDED

#include <cstdarg>
#include <cstring>
#include <exception>

#include "Communicator/communicator.h"
#include "channel.h"

//! MPI-realisation of communicator class BAPI implementation

/**
   Communicator_impl class BAPI provides implementation of communication
   between parallel processes using MPI.
   Interface is basically the same as the abstract communicator class BAPI.
 */

class BAPI Communicator_impl {
 public:
  static int init(int *pargc, char ***pargv);
  static int finalize();
  static void abort();

  static int setup(int ninstance = 1);

  static int setup(const std::vector<int>& lattice_size,
                   std::vector<int>& grid_size,
                   int ninstance = 1);

  // info about rank
  static bool is_primary();

#ifdef ENABLE_MULTI_INSTANCE
  static bool is_primary_master();
#endif

  static int self();   //< rank within small world.
  static int size();   //< size of small world.

#ifdef ENABLE_MULTI_INSTANCE
  static int self_global();
  static int world_id();
#endif

  // world
  static MPI_Comm& world() { return m_comm; }  //!< retrieves current communicator.

  // synchronize
  static int sync();   //< synchronize within small world.

#ifdef ENABLE_MULTI_INSTANCE
  static int sync_global();   //< synchronize all processes.
#endif

  // info
  static double get_time();

  // debug
  static int status();

  // base case
  class BAPI Base {
   public:
    static int reduce(int count, void *recv_buf, void *send_buf, MPI_Datatype type, MPI_Op op, int pattern);
    static int broadcast(size_t size, void *data, int sender);
    static int exchange(size_t size, void *recv_buf, void *send_buf, int idir, int ipm, int tag);

    static int send_1to1(size_t size, void *recv_buf, void *send_buf, int send_to, int recv_from, int tag);
  };

  // for specific datatypes
  static int broadcast_string(int count, string& data, int sender);

  // async communication
  static Channel *send_init(int count, int idir, int ipm);
  static Channel *recv_init(int count, int idir, int ipm);

  // logical and physical layout
  class BAPI Layout;

 private:
  Communicator_impl() {}
  Communicator_impl(const Communicator_impl&) {}
  Communicator_impl& operator=(const Communicator_impl&);

  ~Communicator_impl() {}

#ifdef ENABLE_MULTI_INSTANCE
  static int m_n_instance;  // number of instances
  static int m_instance_id; // id of present instance

  static int m_global_rank;
  static int m_global_size;
#endif

  static int m_grid_rank;
  static int m_grid_size;

  static MPI_Comm m_comm;

  static char default_grid_map[16];
};
#endif /* COMMUNICATOR_MPI_INCLUDED */
