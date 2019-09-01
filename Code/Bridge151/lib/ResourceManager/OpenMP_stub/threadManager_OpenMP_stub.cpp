/*!
        @file    threadManager_OpenMP.cpp

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/
#include "BridgeLib_Private.h"

#if COMMUNICATOR_USE_MPI

#include "ResourceManager/threadManager_OpenMP.h"


#include "Communicator/communicator.h"
#include "IO/bridgeIO.h"
using Bridge::vout;

// These are stub routine when OpenMP is not used.
//                                          [11 Mar 2014 H.Matsufuru]
//====================================================================
// initialization of static member variables.

int ThreadManager_OpenMP::m_Nthread             = 0;
Bridge::VerboseLevel ThreadManager_OpenMP::m_vl = Bridge::CRUCIAL;
std::vector<double>  ThreadManager_OpenMP::m_darray_reduction(0);
std::vector<float>   ThreadManager_OpenMP::m_darray_reductionF(0);

//====================================================================
void ThreadManager_OpenMP::init(int Nthread)
{
  m_vl = CommonParameters::Vlevel();

  m_Nthread = 1;

  vout.general(m_vl, "ThreadManager_OpenMP being setup.\n");
  vout.general(m_vl, "  Number of thread = %d\n", m_Nthread);
  vout.general(m_vl, "  OpenMP is not used: stub implementation.\n");
}


//====================================================================
void ThreadManager_OpenMP::finalize()
{
  vout.paranoiac(m_vl, "Thread manager says good-bye.\n");
}


//====================================================================
int ThreadManager_OpenMP::get_num_threads()
{
  return 1;
}


//====================================================================
int ThreadManager_OpenMP::get_thread_id()
{
  return 0;
}


//====================================================================
void ThreadManager_OpenMP::wait()
{
}


//====================================================================
void ThreadManager_OpenMP::barrier(int Nthread)
{
}


//====================================================================
void ThreadManager_OpenMP::sync_barrier_all()
{
  Communicator::sync();
}


//====================================================================
void ThreadManager_OpenMP::reduce_sum_global(double& a,
                                             const int i_thread, const int Nthread)
{
  a = Communicator::reduce_sum(a);
}


//====================================================================
void ThreadManager_OpenMP::reduce_sum_global(float& a,
                                             const int i_thread, const int Nthread)
{
  a = Communicator::reduce_sum(a);
}


//====================================================================
void ThreadManager_OpenMP::assert_single_thread(
  const std::string& p_class_name)
{ // do nothing.
}


//====================================================================
//============================================================END=====

#endif
