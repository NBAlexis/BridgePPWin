#include "BridgeLib_Private.h"

#if COMMUNICATOR_USE_MPI

/*!
        @file    $Id: communicator_mpi.cpp #$

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate: 2012-11-15 16:42:21 #$

        @version $LastChangedRevision: 1571 $
*/

#include "layout.h"

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

MPI_Comm Communicator_impl::m_comm;

//====================================================================
// class methods
int Communicator_impl::init(int *pargc, char ***pargv)
{
  LOG;

  int is_initialized = 0;
  MPI_Initialized(&is_initialized);

  if (is_initialized) {
    fprintf(stderr, "Communicator: MPI already initialized. skip.\n");
  } else {
#ifdef USE_OPENMP
    int required = MPI_THREAD_FUNNELED;
    int provided;

    MPI_Init_thread(pargc, pargv, required, &provided);

    if (provided < required) {
      fprintf(stderr, "Communicator: MPI not supporting sufficient thread level. exiting.\n");
      exit(EXIT_FAILURE);
    }
#else
    MPI_Init(pargc, pargv);
#endif
  }

  MPI_Comm_size(MPI_COMM_WORLD, &m_global_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &m_global_rank);

  // grid size and rank equal to global ones for a moment until layout is set.
  m_grid_size = m_global_size;
  m_grid_rank = m_global_rank;

  //- initialize m_comm, thanks to Aoyama-san.
  MPI_Comm_dup(MPI_COMM_WORLD, &m_comm);

  return EXIT_SUCCESS;
}


//====================================================================
int Communicator_impl::finalize()
{
  LOG;
  return MPI_Finalize();
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

  m_n_instance = ninstance;

  int gsize = m_global_size / ninstance;
  m_instance_id = m_global_rank / gsize;

  MPI_Comm_split(MPI_COMM_WORLD, m_instance_id, 0 /* key */, &m_comm);
#else
//  m_n_instance = 1;
//  m_instance_id = 0;

  MPI_Comm_dup(MPI_COMM_WORLD, &m_comm);
#endif

  MPI_Comm_size(m_comm, &m_grid_size);
  MPI_Comm_rank(m_comm, &m_grid_rank);

  Communicator_impl::Layout::layout_setup();

  status();

  return EXIT_SUCCESS;
}


//====================================================================
void Communicator_impl::abort()
{
  LOG;
  MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  // unreached.
}


//====================================================================
// information
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


//====================================================================
#ifdef ENABLE_MULTI_INSTANCE

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

//====================================================================
// synchronize
int Communicator_impl::sync()
{
  LOG;
  return MPI_Barrier(m_comm);
}


//====================================================================
#ifdef ENABLE_MULTI_INSTANCE

int Communicator_impl::sync_global()
{
  LOG;
  return MPI_Barrier(MPI_COMM_WORLD);
}
#endif

//====================================================================
// data transfer: base cases
int Communicator_impl::Base::broadcast(size_t size, void *data, int sender)
{
  LOG;
  return MPI_Bcast(data, size, MPI_BYTE, sender, m_comm);
}


//====================================================================
int Communicator_impl::Base::exchange(size_t size, void *recv_buf, void *send_buf, int idir, int ipm, int itag)
{
  LOG;

  MPI_Status status;
  int        p_send, p_recv;

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

  return MPI_Sendrecv(
    send_buf, size, MPI_BYTE, p_send, tag_send,
    recv_buf, size, MPI_BYTE, p_recv, tag_recv,
    m_comm, &status);
}


//====================================================================
int Communicator_impl::Base::send_1to1(size_t size, void *recv_buf, void *send_buf, int send_to, int recv_from, int tag)
{
  LOG;

  MPI_Status status;

  if (send_to == recv_from) {
    memcpy(recv_buf, send_buf, size);
  } else {
    if (self() == recv_from)
      MPI_Send(send_buf, size, MPI_BYTE, send_to, tag, m_comm);

    if (self() == send_to)
      MPI_Recv(recv_buf, size, MPI_BYTE, recv_from, tag, m_comm, &status);
  }

  // sync should be taken outside.

  return EXIT_SUCCESS;
}


//====================================================================
int Communicator_impl::Base::reduce(int count, void *recv_buf, void *send_buf, MPI_Datatype type, MPI_Op op, int pattern)
{
  LOG;
  return MPI_Allreduce((void *)send_buf, (void *)recv_buf, count, type, op, Layout::m_sub_comm[pattern]);
}


//====================================================================
// data transfer for specific datatypes
int Communicator_impl::broadcast_string(int count, string& data, int sender)
{
  LOG;
  assert(count == 1);

  size_t size = 0;

  // broadcast string length.
  if (Communicator::self() == sender) {
    size = data.size();
  }
  MPI_Bcast((void *)&size, sizeof(size_t), MPI_BYTE, sender, m_comm);

  // allocate buffer. pack data at sender.
  char *buf = new char[size + 1];
  memset(buf, '\0', size + 1);

  if (Communicator::self() == sender) {
    data.copy(buf, size, 0);
  }

  // do broadcast.
  int retv = MPI_Bcast((void *)buf, size, MPI_BYTE, sender, m_comm);

  if (Communicator::self() != sender) {
    data = string(buf);
  }

  delete [] buf;

  return retv;
}


//====================================================================
// info
double Communicator_impl::get_time()
{
  return MPI_Wtime();
}


//====================================================================
// debug
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

#endif //COMMUNICATOR_USE_MPI
