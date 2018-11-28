/*!
        @file    $Id: communicator_bgnet.h #$

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 00:49:46 #$

        @version $LastChangedRevision: 1561 $
*/

#ifndef COMMUNICATOR_BGNET_INCLUDED
#define COMMUNICATOR_BGNET_INCLUDED

#include <exception>
#include <cassert>

#include <bgnet.h>

#include "channel.h"

//! Implementation of Communicator with BGNET library.

/*!
  This class is an implementation of Communicator using the
  low latency library BGNET for Blue Gene/Q.
  MPI is not used (while possible to initialize for the purpose
  of using Parallel I/O).
  Multi-instance option is not available.
  If setting ENABLE_MULTI_INSTANCE and putting number of
  instance other than 1, setup() aborts.
                                   [29 Jun 2013 H.Matsufuru]
  Simultaneous use of BGNET and MPI libraries become possible
  thanks to Takumi Doi's contribution. [21 Aug 2014 H.Matsufuru]
 */
class Communicator_impl {
 public:
  static int init(int *pargc, char ***pargv);
  static int finalize();
  static void abort();

  static int setup(int ninstance = 1);

  //! info about rank
  static bool is_primary();

#ifdef ENABLE_MULTI_INSTANCE
  static bool is_primary_master();
#endif

  static int self();  //< rank within small world.
  static int size();  //< size of small world.

#ifdef ENABLE_MULTI_INSTANCE
  static int self_global();
  static int world_id();
#endif


  // static MPI_Comm& world() { return m_comm; }
  static int world() { return m_comm; }

  // synchronize
  static int sync();   //!< synchronize within small world.

#ifdef ENABLE_MULTI_INSTANCE
  static int sync_global();  //!< synchronize all processes.
#endif

  //! for getting time interval using clock count.
  static double get_time();

  //! for debug
  static int status();

  //! base class
  class Base {
   public:
    static int reduce(int count, void *recv_buf, void *send_buf,
                      int type, int op, int pattern);

    static int broadcast(size_t size, void *data, int sender);

    static int exchange(size_t size, void *recv_buf, void *send_buf,
                        int idir, int ipm, int tag);

    static int send_1to1(size_t size, void *recv_buf, void *send_buf,
                         int send_to, int recv_from, int tag);
  };

  //! for specific datatypes
  static int broadcast_string(int count, string& data, int sender);
  static int broadcast_double(int count, double *data, int sender);
  static int broadcast_int(int count, int *data, int sender);

  //! async communication
  static Channel *send_init(int count, int idir, int ipm);
  static Channel *recv_init(int count, int idir, int ipm);

  //! logical and physical layout
  class Layout;

 private:
  Communicator_impl() {}
  Communicator_impl(const Communicator_impl&) {}
  Communicator_impl& operator=(const Communicator_impl&);

  ~Communicator_impl() {}

#ifdef ENABLE_MULTI_INSTANCE
  static int m_n_instance;   //!< number of instances
  static int m_instance_id;  //!< id of present instance

  static int m_global_rank;
  static int m_global_size;
#endif

  static int m_grid_rank;
  static int m_grid_size;

  static int m_comm;      //!< instead of MPI_Comm m_comm;
};
#endif /* #ifndef COMMUNICATOR_BGNET_INCLUDED */
