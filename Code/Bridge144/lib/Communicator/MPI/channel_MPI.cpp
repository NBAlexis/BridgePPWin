#include "BridgeLib_Private.h"

#if COMMUNICATOR_USE_MPI

/*!
        @file    $Id: channel.cpp #$

        @brief

        @author  Hideo Matsufuru  (matsufuru)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate: 2012-11-15 16:42:21 #$

        @version $LastChangedRevision: 1678 $
*/

#include "channel.h"

#include "communicator_mpi.h"
#include "layout.h"

//====================================================================

/**
   Communicator::send_init() method
   creates a channel instance that wraps persistent communication object
   to send data to neighbour node located at upward or downward (specified
   by ipm) in the direction idir.
   count represents size of buffer by number of elements (bytes).
 */
Channel *Communicator_impl::send_init(int count, int idir, int ipm)
{
  LOG;

  assert(ipm == 1 || ipm == -1);
  assert(Layout::m_ndim <= max_dimension);

  Channel *ch = new Channel(count);
  if (!ch) {
#ifdef DEBUG
    printf("%s: allocate channel failed.\n", __func__);
#endif
    exit(EXIT_FAILURE);
  }

  int dest = (ipm == Forward) ? Layout::m_ipe_up[idir] : Layout::m_ipe_dn[idir];
  int tag  = idir + max_dimension * (((ipm == Forward) ? 0 : 1) + 2 * m_grid_rank);

  int retv = MPI_Send_init((void *)&ch->m_buf[0], sizeof(Channel::element_type) * count, MPI_BYTE, dest, tag, m_comm, &ch->m_request);

  if (retv != MPI_SUCCESS) {
    return (Channel *)0;
  }

  return ch;
}


//====================================================================

/**
   Communicator::send_init() method
   creates a channel instance that wraps persistent communication object
   to receive data from neighbour node located at upward or downward (specified
   by ipm) in the direction idir.
   count represents size of buffer by number of elements (bytes).
 */
Channel *Communicator_impl::recv_init(int count, int idir, int ipm)
{
  LOG;

  assert(ipm == 1 || ipm == -1);
  assert(Layout::m_ndim <= max_dimension);

  Channel *ch = new Channel(count);
  if (!ch) {
#ifdef DEBUG
    printf("%s: allocate channel failed.\n", __func__);
#endif
    exit(EXIT_FAILURE);
  }

  int src = (ipm == Forward) ? Layout::m_ipe_up[idir] : Layout::m_ipe_dn[idir];
  int tag = idir + max_dimension * (((ipm == Forward) ? 1 : 0) + 2 * src);

  int retv = MPI_Recv_init((void *)&ch->m_buf[0], sizeof(Channel::element_type) * count, MPI_BYTE, src, tag, m_comm, &ch->m_request);

  if (retv != MPI_SUCCESS) {
    return (Channel *)0;
  }

  return ch;
}


//====================================================================
// class Channel
Channel::Channel() : m_buf(0)
{
  LOG;
}


Channel::Channel(const int count) : m_buf(count)
{
  LOG;
}


Channel::~Channel()
{
  LOG;
}


//====================================================================
int Channel::start()
{
  LOG;
  return MPI_Start(&m_request);
}


//====================================================================
int Channel::wait()
{
  LOG;
  return MPI_Wait(&m_request, &m_status);
}


//====================================================================
// class ChannelSet
ChannelSet::ChannelSet(int count)
 : m_array(count), m_status(count), m_nreq(0)
{
  LOG;
}


//====================================================================
int ChannelSet::append(Channel *p)
{
  LOG;

  if (m_nreq >= m_array.size()) {
    return MPI_ERR_BUFFER;
  }
  m_array[m_nreq++] = p->m_request;  //* need grant access to private data of channel class.

  return MPI_SUCCESS;
}


//====================================================================
int ChannelSet::start()
{
  LOG;
  return MPI_Startall(m_nreq, &m_array[0]);
}


//====================================================================
int ChannelSet::wait()
{
  LOG;
  //  return MPI_Waitall(m_nreq, &m_array[0], (MPI_Status *)0);
  return MPI_Waitall(m_nreq, &m_array[0], &m_status[0]);
}


//====================================================================
//============================================================END=====

#endif //COMMUNICATOR_USE_MPI
