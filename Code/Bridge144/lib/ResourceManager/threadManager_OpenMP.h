/*!
        @file    $Id:: threadManager_OpenMP.h #$

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 18:35:38 #$

        @version $LastChangedRevision: 1571 $
*/

#ifndef THREADMANAGER_OPENMP_INCLUDED
#define THREADMANAGER_OPENMP_INCLUDED

#include <string>

//#include "configure.h"
//#include "defs.h"
#include "Parameters/commonParameters.h"
#include "IO/bridgeIO.h"


//! Thread manager with OpenMP.

/*!
   This class wraps OpenMP API functions.
   The number of available threads is specified as an argument
   of init() function, which is to be called at the begining of
   main, soon after Communicator and CommonParameter are
   initialized.
   If the specified number of threads (nthread) is not the same
   as that of OMP_NUM_THREADS environment variable, the number
   of threads is set to the specified one.
   If Nthread = 0 is specified, the number of threads is set
   to OMP_NUM_THREADS or maximum number of platform.
                                      [29 Aug 2013 H.Matsufuru]
 */
class ThreadManager_OpenMP {
 private:
  static int m_Nthread;              //!< number of threads.
  static Bridge::VerboseLevel m_vl;  //!< verbose level.
  static std::vector<double>  m_darray_reduction;

 private:
  // non-copyable
  ThreadManager_OpenMP(const ThreadManager_OpenMP&);
  ThreadManager_OpenMP& operator=(const ThreadManager_OpenMP&);

 public:
  static const std::string class_name;

  //! setup: called in main only once.
  static void init(int Nthread);

  //! finalization.
  static void finalize();

  //! returns number of threads (works outside of parallel region).
  static int get_num_threads_available() { return m_Nthread; }

  //! returns available number of threads.
  static int get_num_threads();

  //! returns thread id.
  static int get_thread_id();

  //! barrier among threads inside a node.
  static void wait();

  //! barrier among threads inside a node.
  static void barrier(const int Nthread);

  //! barrier among all the threads and nodes.
  static void sync_barrier_all();

  //! global reduction with summation: value is assumed thread local.
  static void reduce_sum_global(double& value,
                                const int i_thread, const int Nthread);

  //! assert currently running on single thread.
  static void assert_single_thread(const std::string& class_name);
};
#endif //THREADMANAGER_OPENMP_INCLUDED
