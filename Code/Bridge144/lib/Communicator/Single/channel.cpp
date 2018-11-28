#include "BridgeLib_Private.h"

#ifdef COMMUNICATOR_USE_SINGLE

/*!
        @file    $Id: channel.cpp #$

        @brief

        @author  Hideo Matsufuru  (matsufuru)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate: 2012-11-15 16:42:21 #$

        @version $LastChangedRevision: 1679 $
*/

#include "channel.h"

//====================================================================
Channel *Communicator::send_init(int count, int idir, int ipm)
{
  LOG;

  assert(ipm == 1 || ipm == -1);

  Channel *ch = new Channel(count);
  if(!ch){
    printf("%s: allocate channel failed.\n", __func__);
    exit(EXIT_FAILURE);
  }

  return ch;

}

//====================================================================
Channel *Communicator::recv_init(int count, int idir, int ipm)
{
  LOG;

  assert(ipm == 1 || ipm == -1);

  Channel *ch = new Channel(count);
  if(!ch){
    printf("%s: allocate channel failed.\n", __func__);
    exit(EXIT_FAILURE);
  }

  return ch;

}

//====================================================================
Channel *Communicator::send_init(int count, int idir, int ipm, void* buf)
{
  LOG;

  assert(ipm == 1 || ipm == -1);

  Channel *ch = new Channel(buf);
  if(!ch){
    printf("%s: allocate channel failed.\n", __func__);
    exit(EXIT_FAILURE);
  }

  return ch;

}

//====================================================================
Channel *Communicator::recv_init(int count, int idir, int ipm, void* buf)
{
  LOG;

  assert(ipm == 1 || ipm == -1);

  Channel *ch = new Channel(buf);
  if(!ch){
    printf("%s: allocate channel failed.\n", __func__);
    exit(EXIT_FAILURE);
  }

  return ch;

}

//====================================================================
// class Channel
Channel::Channel()
  : m_buf(NULL),
    m_my_buf(false)
{
}

Channel::Channel(void* buf)
  : m_buf(reinterpret_cast<element_type*>(buf)),
    m_my_buf(false)
{
}

Channel::Channel(const int count)
{
  m_buf = new element_type[count];
  m_my_buf = true;
}

Channel::~Channel()
{
  if (m_my_buf) delete [] m_buf;

  m_buf = NULL;
  m_my_buf = false;
}

//====================================================================
int Channel::start()
{
  return 0;
}

//====================================================================
int Channel::wait()
{
  return 0;
}

//====================================================================
// class ChannelSet
ChannelSet::ChannelSet(int count)
{
  LOG;
}


//====================================================================
int ChannelSet::append(Channel *p)
{
  return 0;
}


//====================================================================
int ChannelSet::start()
{
  return 0;
}

//====================================================================
int ChannelSet::wait()
{
  return 0;
}

//====================================================================
//============================================================END=====

#endif //COMMUNICATOR_USE_SINGLE
