#include "BridgeLib_Private.h"

#if COMMUNICATOR_USE_MPI

/*!
        @file    $Id:: threadManager_OpenMP.cpp #$

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 18:35:38 #$

        @version $LastChangedRevision: 1571 $
*/

#include "ResourceManager/threadManager_OpenMP.h"

#include <omp.h>

#include "Communicator/communicator.h"
#include "IO/bridgeIO.h"
using Bridge::vout;

//====================================================================
// initialization of static member variables.

int ThreadManager_OpenMP::m_Nthread             = 0;
Bridge::VerboseLevel ThreadManager_OpenMP::m_vl = Bridge::CRUCIAL;
std::vector<double>  ThreadManager_OpenMP::m_darray_reduction(0);

const std::string ThreadManager_OpenMP::class_name = "ThreadManager_OpenMP";

//====================================================================
void ThreadManager_OpenMP::init(int Nthread)
{
  m_vl = CommonParameters::Vlevel();

  vout.general(m_vl, "%s: initialization\n", class_name.c_str());

  int Nthread_env = 0;

#pragma omp parallel
  {
    if (omp_get_thread_num() == 0) {
      Nthread_env = omp_get_num_threads();
    }
  }


  if ((Nthread == Nthread_env) || (Nthread == 0)) {
    m_Nthread = Nthread_env;
  } else {
    vout.general(m_vl, "  Number of threads(env)   = %d\n", Nthread_env);
    vout.general(m_vl, "  Number of threads(input) = %d\n", Nthread);
    vout.general(m_vl, "  resetting Number of threads.\n");
    omp_set_num_threads(Nthread);
    m_Nthread = Nthread;
  }

  vout.general(m_vl, "  Number of thread = %d\n", m_Nthread);

  m_darray_reduction.resize(m_Nthread);
}


//====================================================================
void ThreadManager_OpenMP::finalize()
{
  vout.paranoiac(m_vl, "%s: finalize.\n", class_name.c_str());
}


//====================================================================
int ThreadManager_OpenMP::get_num_threads()
{
  return omp_get_num_threads();
}


//====================================================================
int ThreadManager_OpenMP::get_thread_id()
{
  return omp_get_thread_num();
}


//====================================================================
void ThreadManager_OpenMP::wait()
{
  int Nthread = get_num_threads();

  barrier(Nthread);
}


//====================================================================
void ThreadManager_OpenMP::barrier(int )
{
#pragma omp barrier
}


//====================================================================
void ThreadManager_OpenMP::sync_barrier_all()
{
#pragma omp barrier
#pragma omp master
  {
    Communicator::sync();
  }
#pragma omp barrier
}


//====================================================================
void ThreadManager_OpenMP::reduce_sum_global(double& a,
                                             const int i_thread, const int Nthread)
{
  m_darray_reduction[i_thread] = a;
  barrier(Nthread);

#pragma omp barrier
#pragma omp master
  {
    double b = 0.0;

#ifdef NECSX
#pragma omp flush
#else
#pragma omp flush (m_darray_reduction)
#endif

    for (int i = 0; i < Nthread; ++i) {
      b += m_darray_reduction[i];
    }
    b = Communicator::reduce_sum(b);

#ifdef NECSX
#pragma omp flush
#else
#pragma omp flush (m_darray_reduction)
#endif

    m_darray_reduction[0] = b;
  }

#pragma omp barrier

  a = m_darray_reduction[0];

#pragma omp barrier
}


//====================================================================
void ThreadManager_OpenMP::assert_single_thread(const std::string& name)
{
  int Nthread = get_num_threads();

  if (Nthread != 1) {
    vout.crucial(m_vl, "\n");
    vout.crucial(m_vl, "##### Caution #####\n");
    vout.crucial(m_vl, "Single-thread %s is called in parallel region.\n", name.c_str());
    vout.crucial(m_vl, "Current number of thread = %d.\n", Nthread);

    exit(EXIT_FAILURE);
  }
}


//====================================================================
//============================================================END=====
#endif //COMMUNICATOR_USE_MPI
