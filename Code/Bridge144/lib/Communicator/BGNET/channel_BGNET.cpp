#include "BridgeLib_Private.h"

#if COMMUNICATOR_USE_BGNET

/*!
        @file    $Id: channel.cpp #$

        @brief

        @author  Hideo Matsufuru  (matsufuru)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 00:49:46 #$

        @version $LastChangedRevision: 1561 $
*/

#include "channel.h"

#include "Parameters/commonParameters.h"
#include "IO/bridgeIO.h"
using Bridge::vout;

// BGNET
// #include "bgnet.h"
#include "communicator_bgnet.h"

#include "layout.h"


// class BAPI Channel
// class BAPI Communicator
#define MAX_BUFFER_DIM    8 // this value should be twice dimension.
//====================================================================
// define static members
int Channel::m_set_id = 0;

//std::vector<Channel::container_type> Channel::m_buf_send(MAX_BUFFER_DIM);
//std::vector<Channel::container_type> Channel::m_buf_recv(MAX_BUFFER_DIM);
std::vector<Channel::container_type> Channel::m_buf(2 * MAX_BUFFER_DIM);
//====================================================================
void Channel::set_thread(int Ntask, const std::vector<int>& destid,
                         const std::vector<int>& offset,
                         const std::vector<int>& datasize)
{
  m_bgnet_ids_thread.resize(Ntask);
  assert(destid.size() == Ntask);
  assert(offset.size() == Ntask);
  for (int itask = 0; itask < Ntask; itask++) {
    m_bgnet_ids_thread[itask].fifoID     = itask;
    m_bgnet_ids_thread[itask].sendBufID  = m_bgnet_ids.sendBufID;
    m_bgnet_ids_thread[itask].sendOffset = sizeof(element_type) * offset[itask];
    m_bgnet_ids_thread[itask].size       = datasize[itask];
    m_bgnet_ids_thread[itask].destRank   = m_bgnet_ids.destRank;
    m_bgnet_ids_thread[itask].groupID    = destid[itask];
    m_bgnet_ids_thread[itask].recvBufID  = m_bgnet_ids.recvBufID;
    m_bgnet_ids_thread[itask].recvOffset = sizeof(element_type) * offset[itask];
    m_bgnet_ids_thread[itask].counterID  = m_bgnet_ids.counterID;
  }
}


// class BAPI Communicator
//====================================================================
Channel *Communicator_impl::send_init(int count, int idir, int ipm)
{
  LOG;

  assert(ipm == 1 || ipm == -1);
  assert(Layout::m_ndim <= max_dimension);

  int dest = (ipm == Forward) ? Layout::m_ipe_up[idir] : Layout::m_ipe_dn[idir];
  int tag  = idir + max_dimension * ((ipm == Forward) ? 0 : 1);


  Channel *ch = new Channel(count, tag, Channel::SEND);

  if (!ch) {
#ifdef DEBUG
    printf("%s: allocate channel failed.\n", __func__);
#endif
    exit(EXIT_FAILURE);
  }

  ch->m_bgnet_ids.fifoID     = 0;
  ch->m_bgnet_ids.sendBufID  = Channel::m_set_id + tag;
  ch->m_bgnet_ids.sendOffset = 0;
  ch->m_bgnet_ids.size       = sizeof(Channel::element_type) * count;
  ch->m_bgnet_ids.destRank   = dest;
  ch->m_bgnet_ids.groupID    = 0;
  ch->m_bgnet_ids.recvBufID  = tag;
  ch->m_bgnet_ids.recvOffset = 0;
  ch->m_bgnet_ids.counterID  = tag;

  Communicator::sync();

  return ch;
}


//====================================================================
Channel *Communicator_impl::recv_init(int count, int idir, int ipm)
{
  LOG;

  assert(ipm == 1 || ipm == -1);
  assert(Layout::m_ndim <= max_dimension);

  int src = (ipm == Forward) ? Layout::m_ipe_up[idir] : Layout::m_ipe_dn[idir];
  int tag = idir + max_dimension * ((ipm == Forward) ? 1 : 0);

  Channel *ch = new Channel(count, tag, Channel::RECV);

  if (!ch) {
#ifdef DEBUG
    printf("%s: allocate channel failed.\n", __func__);
#endif
    exit(EXIT_FAILURE);
  }

  ch->m_bgnet_ids.fifoID     = 0;
  ch->m_bgnet_ids.sendBufID  = Channel::m_set_id + tag;
  ch->m_bgnet_ids.sendOffset = 0;
  ch->m_bgnet_ids.size       = sizeof(Channel::element_type) * count;
  ch->m_bgnet_ids.destRank   = -1;
  ch->m_bgnet_ids.groupID    = 0;
  ch->m_bgnet_ids.recvBufID  = tag;
  ch->m_bgnet_ids.recvOffset = 0;
  ch->m_bgnet_ids.counterID  = tag;

  Communicator::sync();

  return ch;
}


// class BAPI Channel
//====================================================================
Channel::Channel(const int count, const int tag, const channel_mode mode)
{
  LOG;
  m_count = count;

  Bridge::VerboseLevel vl = CommonParameters::Vlevel();
  //  Bridge::VerboseLevel vl = Bridge::PARANOIAC;

  if (tag >= MAX_BUFFER_DIM) {
    vout.crucial(vl, "%s: too large value of tag.\n", __func__);
    exit(EXIT_FAILURE);
  }

  if (mode == SEND) {
    m_ibuf = tag;

    if (count > m_buf[m_ibuf].size()) {
      vout.general(vl, "%s: buffer size reset: ibuf = %d  count = %d.\n",
                   __func__, tag, count);
      m_buf[m_ibuf].resize(count);
      //      m_buf = &(m_buf_send[tag][0]);
      BGNET_SetSendBuffer((void *)&m_buf[m_ibuf][0], tag,
                          sizeof(element_type) * count);
      vout.paranoiac(vl, "pointer to m_buf[%d] = %x.\n", m_ibuf, &m_buf[m_ibuf]);
    } else {
      vout.paranoiac(vl, "%s: buffer size unchanged: ibuf = %d  count = %d.\n",
                     __func__, m_ibuf, count);
      vout.paranoiac(vl, "pointer to m_buf[%d] = %x.\n", m_ibuf, &m_buf[m_ibuf]);
    }
  } else if (mode == RECV) {
    m_ibuf = tag + MAX_BUFFER_DIM;

    if (count > m_buf[m_ibuf].size()) {
      vout.general(vl, "%s: buffer size reset: ibuf = %d  count = %d.\n",
                   __func__, tag, count);
      m_buf[m_ibuf].resize(count);
      BGNET_SetRecvBuffer((void *)&m_buf[m_ibuf][0], tag,
                          sizeof(element_type) * count);
      vout.paranoiac(vl, "pointer to m_buf[%d] = %x.\n", m_ibuf, &m_buf[m_ibuf]);
    } else {
      vout.paranoiac(vl, "%s: buffer size unchanged: ibuf = %d  count = %d.\n",
                     __func__, m_ibuf, count);
      vout.paranoiac(vl, "pointer to m_buf[%d] = %x.\n", m_ibuf, &m_buf[m_ibuf]);
    }
  } else {
    vout.crucial(vl, "%s: irrelevant input mode.\n", __func__);
    exit(EXIT_FAILURE);
  }
}


//====================================================================
Channel::~Channel()
{
  LOG;
}


//====================================================================
int Channel::start()
{
  LOG;

  int ret;
  if (m_bgnet_ids.destRank != -1) {
    ret = BGNET_Put(m_bgnet_ids.fifoID,
                    m_bgnet_ids.sendBufID,
                    m_bgnet_ids.sendOffset,
                    m_bgnet_ids.size,
                    m_bgnet_ids.destRank,
                    m_bgnet_ids.groupID,
                    m_bgnet_ids.recvBufID,
                    m_bgnet_ids.recvOffset,
                    m_bgnet_ids.counterID);
  } else {
    ret = 0;
  }

  return ret;
}


//====================================================================
int Channel::wait()
{
  LOG;

  int ret;
  if (m_bgnet_ids.destRank == -1) {
    ret = BGNET_WaitForRecv(m_bgnet_ids.groupID,
                            m_bgnet_ids.counterID,
                            m_bgnet_ids.size);
  } else {
    ret = 0;
  }

  return ret;
}


//====================================================================
int Channel::start_thread(int itask)
{
  LOG;

  int ret;
  if ((m_bgnet_ids_thread[itask].destRank != -1) &&
      (m_bgnet_ids_thread[itask].size != 0)) {
    ret = BGNET_Put(m_bgnet_ids_thread[itask].fifoID,
                    m_bgnet_ids_thread[itask].sendBufID,
                    m_bgnet_ids_thread[itask].sendOffset,
                    m_bgnet_ids_thread[itask].size,
                    m_bgnet_ids_thread[itask].destRank,
                    m_bgnet_ids_thread[itask].groupID,
                    m_bgnet_ids_thread[itask].recvBufID,
                    m_bgnet_ids_thread[itask].recvOffset,
                    m_bgnet_ids_thread[itask].counterID);
  } else {
    ret = 0;
  }

  return ret;
}


//====================================================================
int Channel::wait_thread(int itask)
{
  LOG;

  int ret;
  if ((m_bgnet_ids_thread[itask].destRank == -1) &&
      (m_bgnet_ids_thread[itask].size != 0)) {
    ret = BGNET_WaitForRecv(itask,
                            m_bgnet_ids_thread[itask].counterID,
                            m_bgnet_ids_thread[itask].size);
  } else {
    ret = 0;
  }

  return ret;
}


//====================================================================
//============================================================END=====

#endif //COMMUNICATOR_USE_BGNET
