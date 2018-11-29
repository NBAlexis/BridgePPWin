#include "BridgeLib_Private.h"

#if COMMUNICATOR_USE_BGNET

/*!
        @file    $Id: communicator_bgnet.cpp #$

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate: 2015-10-20 19:54:08 #$

        @version $LastChangedRevision: 1571 $
*/

#include "communicator_bgnet.h"

#include <cstdarg>
#include <cstring>

#ifdef USE_MPI
#include <mpi.h>
#endif

// BGNET: for Blue Gene/Q
#include <hwi/include/bqc/A2_inlines.h>
#define  CLOCKRATE    1.6e9;
// This is temporary setting and valid only for BG/Q.
// If the system provide corresponding constant, to be replaced.

#include "layout.h"

#include "Parameters/commonParameters.h"
#include "IO/bridgeIO.h"
using Bridge::vout;

// exception handler for uncaught throw.
static std::terminate_handler default_handler = std::set_terminate(Communicator::abort);

// define static members
#ifdef ENABLE_MULTI_INSTANCE
int Communicator_impl::m_n_instance  = 1; // number of instances
int Communicator_impl::m_instance_id = 0; // id of present instance
#endif

int Communicator_impl::m_global_rank = 0;
int Communicator_impl::m_global_size = 1;

int Communicator_impl::m_grid_rank = 0;
int Communicator_impl::m_grid_size = 1;

int Communicator_impl::m_comm;

//----------------------------------------------------------------

// class BAPI methods
//====================================================================
int Communicator_impl::init(int *pargc, char ***pargv)
{
  LOG;

#ifdef USE_MPI
  int required = MPI_THREAD_FUNNELED;
  int provided;
  MPI_Init_thread(pargc, pargv, required, &provided);
  if (provided < required) {
    fprintf(stderr,
            "MPI implementation provides insufficient threading support.\n"
            "required = %d, provided = %d : (SINGLE, FUNNELED, SERIALIZED, MULT\
IPLE) = (%d, %d, %d, %d)\n",
            required, provided,
            MPI_THREAD_SINGLE,
            MPI_THREAD_FUNNELED,
            MPI_THREAD_SERIALIZED,
            MPI_THREAD_MULTIPLE);
    abort();
  }
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  // BGNET
  BGNET_Init();

  m_global_size = BGNET_Procs();
  m_global_rank = BGNET_Rank();

  // grid size and rank equal to global ones for a moment until layout is set.
  m_grid_size = m_global_size;
  m_grid_rank = m_global_rank;

  return EXIT_SUCCESS;
}


//====================================================================
int Communicator_impl::finalize()
{
  LOG;

#ifdef USE_MPI
  MPI_Finalize();
#endif

  return EXIT_SUCCESS;
}


//====================================================================
int Communicator_impl::setup(int ninstance)
{
  LOG;

#ifdef ENABLE_MULTI_INSTANCE
  if ((ninstance == 0) || (m_global_size % ninstance != 0)) {
    printf("%s: invalid number of instance: %d\n", "Communicator::init", ninstance);
    abort();
  }

  if (ninstance != 1) {
    printf("%s: sorry, multi-instance is disabled.\n", "Communicator::init");
    abort();
  }

  m_n_instance = ninstance;

  int gsize = m_global_size / ninstance;
  m_instance_id = m_global_rank / gsize;
#else
  m_n_instance  = 1;
  m_instance_id = 0;
#endif

  m_grid_size = m_global_size;
  m_grid_rank = m_global_rank;

  Communicator_impl::Layout::layout_setup();

  status();

  // BGNET parameters
  int N_send_buf = BGNET_GetNumInjBuffer();
  int N_fifo     = BGNET_GetNumInjFIFO();
  int N_recv_buf = BGNET_GetNumRecBuffer();
  int N_group    = BGNET_GetNumGroup();
  int N_counter  = BGNET_GetNumRecCounter();

  Bridge::VerboseLevel vl = CommonParameters::Vlevel();
  vout.general(vl, "BGNET parameters:\n");
  vout.general(vl, "  Number of send buffer: %d\n", N_send_buf);
  vout.general(vl, "  Number of recv buffer: %d\n", N_recv_buf);
  vout.general(vl, "  Number of fifo:        %d\n", N_fifo);
  vout.general(vl, "  Number of group:       %d\n", N_group);
  vout.general(vl, "  Number of counter:     %d\n", N_counter);

  return EXIT_SUCCESS;
}


//====================================================================
void Communicator_impl::abort()
{
  LOG;

#ifdef USE_MPI
  MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
#endif

  abort();

  // unreached.
}


// information
//====================================================================
bool Communicator_impl::is_primary()
{
  return m_grid_rank == 0;
}


//====================================================================
int Communicator_impl::self()
{
  return m_grid_rank;
}


//====================================================================
int Communicator_impl::size()
{
  return m_grid_size;
}


//====================================================================
bool Communicator_impl::is_primary_master()
{
  return m_global_rank == 0;
}


#ifdef ENABLE_MULTI_INSTANCE
//====================================================================
int Communicator_impl::self_global()
{
  return m_global_rank;
}


//====================================================================
int Communicator_impl::world_id()
{
  return m_instance_id;
}
#endif

// synchronize
//====================================================================
int Communicator_impl::sync()
{
  LOG;
  BGNET_GlobalBarrier();
  return 0;
}


#ifdef ENABLE_MULTI_INSTANCE
//====================================================================
int Communicator_impl::sync_global()
{
  LOG;
  BGNET_GlobalBarrier();

  return 0;
}
#endif

// data transfer: base cases
//====================================================================
int Communicator_impl::Base::broadcast(size_t count, void *data, int sender)
{
  LOG;
  BGNET_BCast(data, count, BGNET_COLLECTIVE_BYTE, sender,
              BGNET_COMM_WORLD);

  return 0;
}


//====================================================================
int Communicator_impl::Base::exchange(size_t size,
                                      void *recv_buf, void *send_buf,
                                      int idir, int ipm, int itag)
{
  LOG;

  //  MPI_Status status;
  int p_send, p_recv;

  assert(ipm == 1 || ipm == -1);

  if (Layout::m_grid_dims[idir] == 1) {  // no need to transfer
    memcpy(recv_buf, send_buf, size);
    return EXIT_SUCCESS;
  }

  if (ipm == 1) {  // downward shift
    p_send = Layout::m_ipe_dn[idir];
    p_recv = Layout::m_ipe_up[idir];
  } else {  // upward shift
    p_send = Layout::m_ipe_up[idir];
    p_recv = Layout::m_ipe_dn[idir];
  }

  int tag_send = Layout::tag(self(), idir, -ipm);
  int tag_recv = Layout::tag(p_recv, idir, -ipm);

  return BGNET_Sendrecv(0, send_buf, size, p_send, recv_buf, size, p_recv);
  //  return BGNET_Sendrecv(itag,send_buf,size,p_send,recv_buf,size,p_recv);
}


//====================================================================
int Communicator_impl::Base::send_1to1(size_t size,
                                       void *recv_buf, void *send_buf,
                                       int send_to, int recv_from, int tag)
{
  LOG;

  /*
  MPI_Status status;

  if (send_to == recv_from) {

    memcpy(recv_buf, send_buf, size);

  } else {

    if (self() == recv_from)
      MPI_Send(send_buf, size, MPI_BYTE, send_to, tag, MPI_COMM_WORLD);

    if (self() == send_to)
      MPI_Recv(recv_buf, size, MPI_BYTE, recv_from, tag, MPI_COMM_WORLD, &status);

  };
  */

  if (send_to == recv_from) {
    memcpy(recv_buf, send_buf, size);
  } else {
    if (self() == recv_from)
      BGNET_Send(0, send_buf, size, send_to);

    if (self() == send_to)
      BGNET_Recv(0, recv_buf, size, recv_from);
  }


  // sync should be taken outside.

  return EXIT_SUCCESS;
}


//====================================================================
int Communicator_impl::Base::reduce(int count,
                                    void *recv_buf, void *send_buf,
                                    int type, int op, int pComm)
{
  LOG;

  BGNET_AllReduce(send_buf, recv_buf, count, type, op, BGNET_COMM_WORLD);

  return EXIT_SUCCESS;
}


// data transfer for specific datatypes
//====================================================================
int Communicator_impl::broadcast_double(int count, double *data,
                                        int sender)
{
  LOG;

  BGNET_BCast(data, count, BGNET_COLLECTIVE_DOUBLE, sender,
              BGNET_COMM_WORLD);

  return EXIT_SUCCESS;
}


//====================================================================
int Communicator_impl::broadcast_int(int count, int *data, int sender)
{
  LOG;

  BGNET_BCast(data, count, BGNET_COLLECTIVE_INT32, sender,
              BGNET_COMM_WORLD);

  return EXIT_SUCCESS;
}


//====================================================================
int Communicator_impl::broadcast_string(int count, string& data,
                                        int sender)
{
  LOG;
  assert(count == 1);

  size_t size = 0;

  // broadcast string length.
  if (Communicator::self() == sender) {
    size = data.size();
  }

  if (sizeof(size_t) == 4) {
    BGNET_BCast(&size, 1, BGNET_COLLECTIVE_INT32, sender,
                BGNET_COMM_WORLD);
  } else if (sizeof(size_t) == 8) {
    BGNET_BCast(&size, 1, BGNET_COLLECTIVE_INT64, sender,
                BGNET_COMM_WORLD);
  } else {
    abort();
  }

  // allocate buffer. pack data at sender.

  /*
  char *buf = new char[size+1];
  memset(buf, '\0', size+1);
  */
  int size2 = (size + 1) / 4;
  if ((size + 1) % 4 != 0) size2 += 1;
  size2 *= 4;
  char *buf = new char[size2];
  memset(buf, '\0', size2);


  if (Communicator::self() == sender) {
    data.copy(buf, size, 0);
  }

  // do broadcast.
  //  int retv = MPI_Bcast((void*)buf, size, MPI_BYTE, sender, m_comm);
  BGNET_BCast((int *)buf, size2 / 4, BGNET_COLLECTIVE_INT32, sender,
              BGNET_COMM_WORLD);

  if (Communicator::self() != sender) {
    //data = string(buf);
    data.assign(buf, size);
  }

  delete [] buf;

  //  return retv;
  return EXIT_SUCCESS;
}


// info
//====================================================================
double Communicator_impl::get_time()
{
  uint64_t clock = GetTimeBase();
  double   dtime = clock / CLOCKRATE;

  return dtime;
}


// debug
//====================================================================
int Communicator_impl::status()
{
#ifdef DEBUG
#ifdef ENABLE_MULTI_INSTANCE
  printf("global_rank=%2d/%2d: ngrid=%d, grid_id=%d: grid_rank=%2d/%2d\n",
         m_global_rank, m_global_size,
         m_n_instance, m_instance_id,
         m_grid_rank, m_grid_size);
#else
  printf("grid_rank=%2d/%2d\n",
         m_grid_rank, m_grid_size);
#endif
#endif

  return EXIT_SUCCESS;
}


//====================================================================
//============================================================END=====

#endif //COMMUNICATOR_USE_BGNET
