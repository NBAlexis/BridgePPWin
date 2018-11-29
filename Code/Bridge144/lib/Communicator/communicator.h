/*!
        @file    $Id:: communicator.h #$

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2017-11-13 18:52:57 #$

        @version $LastChangedRevision: 1679 $
*/

#ifndef COMMUNICATOR_INCLUDED
#define COMMUNICATOR_INCLUDED

#include "configure.h"
#include "defs.h"

#include <cstdio>
#include <cstdlib>
#include <cstddef>
#include <string>
using std::string;

class BAPI Channel;

//! Communication library which wraps MPI.

/**
  This class BAPI provides communication library which wraps
  MPI (message passing interface) if the implementation
  file communicator_mpi.cpp is bound.
  For single processor environment with no MPI library,
  communicator_dummy.cpp should be bound instead.
  [28 Dec 2011 H.Matsufuru]

  This class BAPI defines interface of inter-node communication routines.
  All methods are static, i.e. class-methods (like global).
  The explicit definitions are hidden in the implementation classes.
*/

// forward declaration
//class BAPI Channel;

class BAPI Communicator {
 public:
  //! initialize communicator

  /**
     initialize communicator.
     @param pargc pointer to argc passed from main function.
     @param pargv pointer to argv passed from main function.

     they would further be passed to communication library
     for hints to machine-dependent configurations.
  */
  static int init(int *pargc, char ***pargv);

  //! finalize communicator

  /**
     finalize communicator.
     terminates communication environment.
  */
  static int finalize();

  //! terminate communicator

  /**
     terminate communicator immediately
     at some erroneous situations.
  */
  static void abort();

  //! setup communicator

  /**
     setup communicator environment such as logical layout.
     called after the parameters are obtained.
     @param ninstance specifies multiplicity of trivial parallelism.
     (expects 1 at present).
   */
  static int setup(int ninstance = 1);

// info about rank
  static bool is_primary();        //!< check if the present node is primary in small communicator.
  static bool is_primary_master(); //!< check if the present node is primary in global communicator.

  static int self();               //!< rank within small world.

  static int nodeid() { return self(); }  //!< alternative name for self().
  static int size();   //!< size of small world.

#ifdef ENABLE_MULTI_INSTANCE
  static int self_global(); //!< rank within global communicator.
  static int world_id();    //!< id for the current small world.
#endif

// layout
  static int ipe(const int dir);                          //!< logical coordinate of current proc.
  static int npe(const int dir);                          //!< logical grid extent

  static int grid_rank(int *rank, const int *grid_coord); //!< find rank number from grid coordinate.
  static int grid_coord(int *grid_coord, const int rank); //!< find grid coordinate from rank number.
  static int grid_dims(int *grid_dims);                   //!< find grid dimensions.

// synchronize
  static int sync();   //!< synchronize within small world.

#ifdef ENABLE_MULTI_INSTANCE
  static int sync_global();   //!< synchronize all processes.
#endif

// data transfer
  static int broadcast(int count, double *data, int sender);                                          //!< broadcast array of double from sender.
  static int broadcast(int count, int *data, int sender);                                             //!< broadcast array of integer from sender.
  static int broadcast(int count, string& data, int sender);                                          //!< broadcast a string from sender. count is insignificant.

  static int exchange(int count, double *recv_buf, double *send_buf, int idir, int ipm, int tag);     //!< receive array of double from upstream specified by idir and ipm, and send array to downstream.
  static int exchange(int count, int *recv_buf, int *send_buf, int idir, int ipm, int tag);           //!< receive array of int from upstream specified by idir and ipm, and send array to downstream.

  static int send_1to1(int count, double *recv_buf, double *send_buf, int p_to, int p_from, int tag); //!< send array of double from rank p_from to rank p_to. communication distinguished by tag.
  static int send_1to1(int count, int *recv_buf, int *send_buf, int p_to, int p_from, int tag);       //!< send array of int from rank p_from to rank p_to. communication distinguished by tag.

  static int reduce_sum(int count, double *recv_buf, double *send_buf, int pattern = 0);              //!< make a global sum of an array of double over the communicator. pattern specifies the dimensions to be reduced.
  static int reduce_sum(int count, int *recv_buf, int *send_buf, int pattern = 0);                    //!< make a global sum of an array of int over the communicator. pattern specifies the dimensions to be reduced.

  static int reduce_max(int count, double *recv_buf, double *send_buf, int pattern = 0);              //!< find a global maximum of an array of double over the communicator. pattern specifies the dimensions to be reduced.
  static int reduce_max(int count, int *recv_buf, int *send_buf, int pattern = 0);                    //!< find a global maximum of an array of int over the communicator. pattern specifies the dimensions to be reduced.

  static int reduce_min(int count, double *recv_buf, double *send_buf, int pattern = 0);              //!< find a global minimum of an array of double over the communicator. pattern specifies the dimensions to be reduced.
  static int reduce_min(int count, int *recv_buf, int *send_buf, int pattern = 0);                    //!< find a global minimum of an array of int over the communicator. pattern specifies the dimensions to be reduced.

  static double reduce_sum(double);                                                                   //!< alternative interface to reduce_sum(). returns the global sum of a double over the whole communicator.
  static double reduce_max(double);                                                                   //!< alternative interface to reduce_max(). returns the global sum of a double over the whole communicator.
  static double reduce_min(double);                                                                   //!< alternative interface to reduce_min(). returns the global sum of a double over the whole communicator.

  static double get_time();                                                                           //!< obtain a wall-clock time.

  // async communication
  static Channel *send_init(int count, int idir, int ipm);
  static Channel *recv_init(int count, int idir, int ipm);

  // async communication with given buffer [2017.09.02 H.Matsufuru]
  static Channel *send_init(int count, int idir, int ipm, void* buf);
  static Channel *recv_init(int count, int idir, int ipm, void* buf);

// debug
  static int status();

//! base case

  /**
     this class BAPI defines base case for data exchange that are specified
     by plain streams of bytes of size.
   */
  class BAPI Base {
   public:
    static int broadcast(size_t size, void *data, int sender);
    static int exchange(size_t size, void *recv_buf, void *send_buf, int idir, int ipm, int tag);

    static int send_1to1(size_t size, void *recv_buf, void *send_buf, int send_to, int recv_from, int tag);
  };

 private:

//! no instance at all

/**
   communicator class BAPI is not indented to be instantiated.
   constructor, copy constroctor, assignment operator, and destructor
   are defined as private.
*/
  Communicator() {}
  Communicator(const Communicator&) {}
  Communicator& operator=(const Communicator&);

  ~Communicator() {}
};
#endif /* COMMUNICATOR_INCLUDED */
